#!/usr/bin/env bash
#
# aws_run_gpu.sh — run a single ANUGA GPU job on AWS with no standing infra.
#
# Launches ONE GPU EC2 instance that pulls the anuga:gpu container, runs your
# command with S3 input/output, uploads results, and TERMINATES itself. You pay
# only for the instance-hours the job uses (it self-destructs when done).
#
# Prereqs: awscli configured (`aws configure`); an S3 bucket you can read/write;
# a one-time GPU vCPU service-quota increase for the instance family (new
# accounts start at 0 for G/P instances).
#
# Typical use — upload your project (run script + data), then launch:
#   ./aws_run_gpu.sh \
#       --upload   ./my_tohoku_project \
#       --input    s3://my-bucket/anuga/in \
#       --output   s3://my-bucket/anuga/out \
#       --command  "python run_Tohoku.py -alg DE0"
#
# Results (and anuga-run.log) appear under --output. The instance is gone.
#
# Add --dry-run to validate creds / GPU quota / AMI and print the launch plan
# without creating or charging anything.
#
# Image: defaults to the public GHCR image. Pass --ecr to use the private
# in-region ECR image <acct>.dkr.ecr.<region>.amazonaws.com/anuga:gpu-slim
# instead (~459 MB slim GPU build — far faster cold-start than the ~16.7 GB
# public one). The instance logs in to ECR with its role, and the auto-created
# role is granted ECR read. Any --image whose URI is an ECR registry is
# auto-detected and gets the same login. (Push the image first — see docker/README.md.)
#
# Robustness: creates a default VPC if the account has none (or pass --subnet-id);
# on InsufficientInstanceCapacity it falls back across instance types and AZs
# (--instance takes a single type or a comma list). The run log records timing
# (image-pull vs container vs total).
#
# NOTE: authored but not yet exercised against a live account — test with your
# own creds and a tiny run first.
set -euo pipefail

# ---- defaults ---------------------------------------------------------------
IMAGE="ghcr.io/anuga-community/anuga:develop-gpu"   # pre-release (built from develop)
USE_ECR=0                          # --ecr: use the private in-region ECR image
ECR_REPO="anuga:gpu-slim"          # repo:tag used to build the ECR URI for --ecr
# Comma-separated capacity-fallback chain: the first type with capacity (in any
# AZ) wins. g5.2xlarge (A10G/cc86, 8 vCPU) -> g5.xlarge -> g4dn.xlarge (T4/cc75,
# cheap + widely available). Override with --instance (single or comma list).
INSTANCE_TYPES="g5.2xlarge,g5.xlarge,g4dn.xlarge"
REGION="${AWS_REGION:-${AWS_DEFAULT_REGION:-}}"
DISK_GB=100                       # root volume; DEMs + SWW output can be large
COMPUTE_MODE="unified"            # GPU offload for mode-2 domains
INPUT_S3="" ; OUTPUT_S3="" ; COMMAND="" ; UPLOAD_DIR="" ; OUTPUT_DIR=""
INSTANCE_PROFILE="" ; AMI="" ; SSH_KEY="" ; SPOT=0 ; KEEP=0 ; DRY=0 ; SUBNET_ID=""

usage() {
  sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//'
  exit "${1:-0}"
}

# ---- args -------------------------------------------------------------------
while [ $# -gt 0 ]; do
  case "$1" in
    --image)            IMAGE="$2"; shift 2 ;;
    --ecr)              USE_ECR=1; shift ;;              # use ECR image built from account+region
    --ecr-repo)         ECR_REPO="$2"; USE_ECR=1; shift 2 ;;  # override repo:tag (implies --ecr)
    --instance)         INSTANCE_TYPES="$2"; shift 2 ;;   # single or comma list
    --subnet-id)        SUBNET_ID="$2"; shift 2 ;;        # else default-VPC subnets
    --region)           REGION="$2"; shift 2 ;;
    --disk)             DISK_GB="$2"; shift 2 ;;
    --input)            INPUT_S3="$2"; shift 2 ;;
    --output)           OUTPUT_S3="$2"; shift 2 ;;
    --command)          COMMAND="$2"; shift 2 ;;
    --upload)           UPLOAD_DIR="$2"; shift 2 ;;
    --output-dir)       OUTPUT_DIR="$2"; shift 2 ;;   # ANUGA_OUTPUT_DIR in container
    --compute-mode)     COMPUTE_MODE="$2"; shift 2 ;;
    --instance-profile) INSTANCE_PROFILE="$2"; shift 2 ;;
    --ami)              AMI="$2"; shift 2 ;;
    --ssh-key)          SSH_KEY="$2"; shift 2 ;;
    --spot)             SPOT=1; shift ;;
    --keep)             KEEP=1; shift ;;              # don't self-terminate (debug)
    --dry-run)          DRY=1; shift ;;              # validate + print plan; launch nothing
    -h|--help)          usage 0 ;;
    *) echo "unknown arg: $1" >&2; usage 1 ;;
  esac
done

die() { echo "ERROR: $*" >&2; exit 1; }
[ -n "$OUTPUT_S3" ] || die "--output s3://... is required"
[ -n "$COMMAND" ]   || die "--command '...' is required"
[ -n "$REGION" ]    || die "no region (pass --region or set AWS_REGION)"
[ -n "$UPLOAD_DIR" ] && [ -z "$INPUT_S3" ] && die "--upload needs --input to upload to"
AWS=(aws --region "$REGION")

# ---- resolve ECR image (--ecr) and detect ECR registries --------------------
# --ecr builds <account>.dkr.ecr.<region>.amazonaws.com/<repo:tag> from the
# caller's account. Any ECR image (via --ecr or an ECR --image) needs the
# instance to `docker login` first (private registry) and the role to have ECR
# read — both handled below.
if [ "$USE_ECR" = "1" ]; then
  ACCOUNT=$("${AWS[@]}" sts get-caller-identity --query Account --output text 2>/dev/null) \
    || die "--ecr: could not resolve AWS account id (is the CLI configured?)"
  IMAGE="${ACCOUNT}.dkr.ecr.${REGION}.amazonaws.com/${ECR_REPO}"
  echo ">> --ecr: using image $IMAGE"
fi
case "$IMAGE" in
  *.dkr.ecr.*.amazonaws.com/*) IMAGE_IS_ECR=1 ; ECR_REGISTRY="${IMAGE%%/*}" ;;
  *)                           IMAGE_IS_ECR=0 ; ECR_REGISTRY="" ;;
esac

bucket_of() { echo "$1" | sed -E 's#^s3://([^/]+).*#\1#'; }

# ---- optional: upload a local project dir to the input prefix ---------------
if [ -n "$UPLOAD_DIR" ]; then
  if [ "$DRY" = "1" ]; then
    echo "[dry-run] would upload ${UPLOAD_DIR} -> ${INPUT_S3}"
  else
    echo ">> uploading ${UPLOAD_DIR} -> ${INPUT_S3}"
    aws s3 sync "$UPLOAD_DIR" "$INPUT_S3"
  fi
fi

# ---- resolve a GPU-ready AMI (Docker + NVIDIA toolkit + driver preinstalled) -
if [ -z "$AMI" ]; then
  echo ">> resolving Deep Learning Base GPU AMI via SSM"
  AMI=$("${AWS[@]}" ssm get-parameter \
    --name /aws/service/deeplearning/ami/x86_64/base-oss-nvidia-driver-gpu-ubuntu-22.04/latest/ami-id \
    --query 'Parameter.Value' --output text 2>/dev/null) \
    || die "could not resolve the DLAMI AMI id via SSM; pass --ami ami-xxxxx (a GPU AMI with Docker + the NVIDIA Container Toolkit)"
fi
echo ">> AMI: $AMI   instances: $INSTANCE_TYPES   region: $REGION"

# ---- ensure an instance profile with S3 access ------------------------------
if [ -z "$INSTANCE_PROFILE" ]; then
  INSTANCE_PROFILE="anuga-gpu-runner"
  ROLE="$INSTANCE_PROFILE"
  if aws iam get-instance-profile --instance-profile-name "$INSTANCE_PROFILE" >/dev/null 2>&1; then
    : # already exists
  elif [ "$DRY" = "1" ]; then
    echo "[dry-run] would create IAM instance profile '$INSTANCE_PROFILE' (S3 access to the input/output buckets)"
  else
    echo ">> creating IAM instance profile '$INSTANCE_PROFILE' (S3 access)"
    IN_B=$(bucket_of "${INPUT_S3:-$OUTPUT_S3}") ; OUT_B=$(bucket_of "$OUTPUT_S3")
    aws iam create-role --role-name "$ROLE" \
      --assume-role-policy-document '{"Version":"2012-10-17","Statement":[{"Effect":"Allow","Principal":{"Service":"ec2.amazonaws.com"},"Action":"sts:AssumeRole"}]}' >/dev/null
    aws iam put-role-policy --role-name "$ROLE" --policy-name s3-access \
      --policy-document "{\"Version\":\"2012-10-17\",\"Statement\":[{\"Effect\":\"Allow\",\"Action\":[\"s3:GetObject\",\"s3:PutObject\",\"s3:DeleteObject\",\"s3:ListBucket\"],\"Resource\":[\"arn:aws:s3:::${IN_B}\",\"arn:aws:s3:::${IN_B}/*\",\"arn:aws:s3:::${OUT_B}\",\"arn:aws:s3:::${OUT_B}/*\"]}]}" >/dev/null
    aws iam create-instance-profile --instance-profile-name "$INSTANCE_PROFILE" >/dev/null
    aws iam add-role-to-instance-profile --instance-profile-name "$INSTANCE_PROFILE" --role-name "$ROLE" >/dev/null
    echo ">> waiting for IAM propagation..."; sleep 15
  fi
  # ECR pulls need the role to read the registry. attach-role-policy is
  # idempotent, so this also fixes a role first created for an S3-only run.
  if [ "$IMAGE_IS_ECR" = "1" ]; then
    if [ "$DRY" = "1" ]; then
      echo "[dry-run] would attach AmazonEC2ContainerRegistryReadOnly to role '$ROLE'"
    else
      echo ">> granting role '$ROLE' ECR read (AmazonEC2ContainerRegistryReadOnly)"
      aws iam attach-role-policy --role-name "$ROLE" \
        --policy-arn arn:aws:iam::aws:policy/AmazonEC2ContainerRegistryReadOnly >/dev/null 2>&1 || true
    fi
  fi
fi

# ---- resolve subnet(s) to launch into ---------------------------------------
# EC2 places the instance in the default VPC when no subnet is given; if the
# account has no default VPC (a common surprise), run-instances fails with
# VPCIdNotSpecified. Resolve one subnet per AZ so the capacity loop below can
# also try alternate AZs. --subnet-id overrides all of this.
SUBNETS=()
if [ -n "$SUBNET_ID" ]; then
  SUBNETS=("$SUBNET_ID")
else
  VPC=$("${AWS[@]}" ec2 describe-vpcs --filters Name=isDefault,Values=true \
        --query 'Vpcs[0].VpcId' --output text 2>/dev/null || echo None)
  if [ "$VPC" = "None" ] || [ -z "$VPC" ]; then
    if [ "$DRY" = "1" ]; then
      echo "[dry-run] no default VPC in $REGION — would run 'aws ec2 create-default-vpc'"
    else
      echo ">> no default VPC in $REGION — creating one"
      "${AWS[@]}" ec2 create-default-vpc >/dev/null \
        || die "no default VPC and create-default-vpc failed; pass --subnet-id of an existing subnet"
      VPC=$("${AWS[@]}" ec2 describe-vpcs --filters Name=isDefault,Values=true \
            --query 'Vpcs[0].VpcId' --output text)
    fi
  fi
  if [ -n "${VPC:-}" ] && [ "$VPC" != "None" ]; then
    mapfile -t SUBNETS < <("${AWS[@]}" ec2 describe-subnets \
      --filters Name=vpc-id,Values="$VPC" Name=default-for-az,Values=true \
      --query 'Subnets[].SubnetId' --output text | tr '\t' '\n' | grep -v '^$')
  fi
fi
[ ${#SUBNETS[@]} -eq 0 ] && SUBNETS=("<none>")   # dry-run w/o default VPC

# ---- build user-data (runs at boot on the GPU instance) ---------------------
TERMINATE="shutdown -h now"
[ "$KEEP" = "1" ] && TERMINATE="echo '[--keep] leaving instance up; terminate it manually when done.'"
USER_DATA=$(cat <<EOF
#!/bin/bash
set -x
exec > /var/log/anuga-run.log 2>&1
T0=\$(date +%s)
${ECR_REGISTRY:+aws ecr get-login-password --region ${REGION} | docker login --username AWS --password-stdin ${ECR_REGISTRY}}
docker pull "${IMAGE}" || true
TR=\$(date +%s)
docker run --rm --gpus all \
  -e ANUGA_DEFAULT_COMPUTE_MODE='${COMPUTE_MODE}' \
  ${INPUT_S3:+-e INPUT_S3='${INPUT_S3}'} \
  -e OUTPUT_S3='${OUTPUT_S3}' \
  ${OUTPUT_DIR:+-e ANUGA_OUTPUT_DIR='${OUTPUT_DIR}'} \
  "${IMAGE}" sh -c '${COMMAND}'
rc=\$?
TE=\$(date +%s)
echo "anuga timing: image-pull \$((TR-T0))s | container \$((TE-TR))s | total(boot-script) \$((TE-T0))s | rc=\$rc"
aws s3 cp /var/log/anuga-run.log "${OUTPUT_S3%/}/anuga-run.log" || true
${TERMINATE}
EOF
)

# ---- launch (capacity fallback across instance types x AZs) -----------------
# Static args; --instance-type and --subnet-id are added per attempt below.
RUN_ARGS=(
  --image-id "$AMI"
  --count 1
  --iam-instance-profile "Name=$INSTANCE_PROFILE"
  --instance-initiated-shutdown-behavior terminate
  --metadata-options "HttpEndpoint=enabled,HttpTokens=required,HttpPutResponseHopLimit=2"
  --block-device-mappings "DeviceName=/dev/sda1,Ebs={VolumeSize=${DISK_GB},VolumeType=gp3}"
  --user-data "$USER_DATA"
  --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=anuga-gpu-run}]"
)
[ -n "$SSH_KEY" ] && RUN_ARGS+=(--key-name "$SSH_KEY")
[ "$SPOT" = "1" ] && RUN_ARGS+=(--instance-market-options "MarketType=spot")
# hop limit 2 lets the container reach IMDSv2 for the instance-profile creds the
# entrypoint's `aws s3 sync` uses.

if [ "$DRY" = "1" ]; then
  echo
  echo "================= DRY RUN — nothing launched, nothing charged ================="
  echo "-- identity --"
  aws sts get-caller-identity --output text 2>&1 | sed 's/^/  /' || echo "  (could not verify creds — is the CLI configured?)"
  echo "-- GPU vCPU quotas (each needs >= this instance's vCPUs) --"
  for qq in "L-DB2E81BA:G/VT on-demand" "L-417A185B:P on-demand" \
            "L-3819A6DF:G/VT spot" "L-7212CCBC:P spot"; do
    qcode=${qq%%:*}; qname=${qq##*:}
    q=$("${AWS[@]}" service-quotas get-service-quota --service-code ec2 --quota-code "$qcode" \
          --query 'Quota.Value' --output text 2>/dev/null) \
      && printf '  %-16s %s vCPUs\n' "$qname" "$q" \
      || printf '  %-16s (could not read — see AWS_SETUP.md step 6)\n' "$qname"
  done
  [ "$SPOT" = "1" ] && echo "  (--spot: tries the spot quota first, then falls back to on-demand)"
  echo "-- plan --"
  printf '  %-16s %s\n' image "$IMAGE" instances "$INSTANCE_TYPES (spot=$SPOT, disk=${DISK_GB}GB)" \
    region "$REGION" ami "$AMI" profile "$INSTANCE_PROFILE" subnets "${SUBNETS[*]}" \
    input "${INPUT_S3:-<none>}" output "$OUTPUT_S3" command "$COMMAND" compute-mode "$COMPUTE_MODE"
  echo "-- ec2 run-instances permission check --"
  if aws iam get-instance-profile --instance-profile-name "$INSTANCE_PROFILE" >/dev/null 2>&1; then
    dt="${INSTANCE_TYPES%%,*}"; ds="${SUBNETS[0]}"
    DA=( "${RUN_ARGS[@]}" --instance-type "$dt" )
    [ "$ds" != "<none>" ] && DA+=( --subnet-id "$ds" )
    "${AWS[@]}" ec2 run-instances "${DA[@]}" --dry-run 2>&1 | sed 's/^/  /' | head -3
  else
    echo "  (skipped — instance profile not created yet; it would be created on a real run)"
  fi
  echo "==============================================================================="
  echo "Re-run without --dry-run to launch."
  exit 0
fi

# With --spot we try the spot market first and fall back to on-demand, rather
# than failing outright: spot is much cheaper but its capacity pools are
# separate and often empty, and for these short jobs an on-demand run is far
# better than no run. Without --spot only on-demand is attempted.
MARKETS=("on-demand")
[ "$SPOT" = "1" ] && MARKETS=("spot" "on-demand")

echo ">> launching (types: $INSTANCE_TYPES; ${#SUBNETS[@]} AZ subnet(s); markets: ${MARKETS[*]})..."
IID="" ; USED_TYPE="" ; USED_SUBNET="" ; USED_MARKET=""
for market in "${MARKETS[@]}"; do
  MARKET_ARGS=()
  [ "$market" = "spot" ] && MARKET_ARGS=(--instance-market-options "MarketType=spot")
  for itype in ${INSTANCE_TYPES//,/ }; do
    for subnet in "${SUBNETS[@]}"; do
      ATTEMPT=( "${RUN_ARGS[@]}" "${MARKET_ARGS[@]}" --instance-type "$itype" )
      [ "$subnet" != "<none>" ] && ATTEMPT+=( --subnet-id "$subnet" )
      if out=$("${AWS[@]}" ec2 run-instances "${ATTEMPT[@]}" \
                 --query 'Instances[0].InstanceId' --output text 2>&1); then
        IID="$out" ; USED_TYPE="$itype" ; USED_SUBNET="$subnet" ; USED_MARKET="$market" ; break 3
      elif echo "$out" | grep -qE "InsufficientInstanceCapacity|Unsupported|InvalidParameterValue|MaxSpotInstanceCountExceeded"; then
        echo "   $itype${subnet:+ / $subnet} [$market]: unavailable — trying next"
        continue
      else
        die "run-instances failed: $out"
      fi
    done
  done
  [ "$market" = "spot" ] && echo ">> no spot capacity/quota — falling back to on-demand"
done
[ -z "$IID" ] && die "no capacity for any of [$INSTANCE_TYPES] across ${#SUBNETS[@]} AZ(s) in ${MARKETS[*]}; try --instance <other type>, another --region, or retry later"

# Report the market we ACTUALLY got, not the one requested: with --spot the
# launch may have fallen back to on-demand, and the price difference matters.
SPOT_LABEL="" ; [ "$USED_MARKET" = "spot" ] && SPOT_LABEL=", spot"
echo ">> instance $IID launched ($USED_TYPE${SPOT_LABEL}${USED_SUBNET:+, $USED_SUBNET})."
echo ">> it will run: ${COMMAND}"
echo ">> results ->   ${OUTPUT_S3}    (+ anuga-run.log with timing)"
if [ "$KEEP" = "1" ]; then
  echo ">> --keep set: instance will NOT self-terminate. Terminate with:"
  echo "     aws --region $REGION ec2 terminate-instances --instance-ids $IID"
else
  echo ">> instance self-terminates when the job finishes. Watch it:"
  echo "     aws --region $REGION ec2 describe-instances --instance-ids $IID --query 'Reservations[0].Instances[0].State.Name' --output text"
fi
