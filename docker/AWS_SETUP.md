# Setting up AWS to run ANUGA GPU jobs

A practical, ordered walkthrough from "no account" to running `aws_run_gpu.sh`.
Two GPU-specific gotchas are called out — **a service-quota increase** and
**billing guardrails** — do those early.

> **Cost note:** GPU instances are **not** free-tier. You pay from the first run
> (~US$1.2/hr for a `g5.2xlarge` on-demand, ~$0.4 spot). Set up the budget in
> step 5 before running anything. Each user runs in **their own AWS account and
> pays for their own usage** — there is no shared/central bill.

## 1. Create the account
- aws.amazon.com → **Create an AWS Account**. Needs an email, a **credit card**,
  and phone verification.
- The signup email/password is the **root user** (all-powerful). Use it only for
  account/billing setup, then stop using it day-to-day.

## 2. Lock down the root user
- Sign in as root → **IAM → enable MFA** on the root user (authenticator app).
  This is the single most important security step.
- Do **not** create access keys for the root user.

## 3. Create your everyday identity
- **Simplest (solo account):** IAM → **Users → Create user** → attach
  **AdministratorAccess** (fine in your own account) → create an **access key**
  (type: CLI). Save the Key ID + Secret. Enable MFA on this user too.
- **More secure/modern:** IAM **Identity Center** → user + admin permission set,
  then sign in with `aws configure sso`. Either works with the launch script.

## 4. Install and configure the AWS CLI
```bash
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o awscliv2.zip
unzip awscliv2.zip && sudo ./aws/install
aws configure                 # paste Key ID + Secret, set region + output=json
aws sts get-caller-identity   # confirms creds work
```

## 5. Set a budget + billing alert — **before running jobs**
- Root → **Account → "IAM user and role access to billing"** = Activate.
- **Billing → Budgets → Create budget** → monthly cost budget (e.g. $50) with an
  email alert at 80% / 100%. The first two budgets are free. This is your safety
  net against a forgotten/runaway instance.

## 6. Request the GPU quota increase — **do this early (can take hours–days)**
New accounts have a **0-vCPU quota** for GPU instances, so the first GPU launch
*fails* until this is approved.
- **Service Quotas → Amazon EC2 →** search **"Running On-Demand G and P
  instances"** → request an increase to at least **8 vCPUs** (a `g5.2xlarge` =
  8 vCPU; ask for 16–32 for headroom).
- Using Spot too? Also raise **"All G and P Spot Instance Requests"**.

## 7. Pick a region + make an S3 bucket
- Choose a region near you with GPU availability. In Australia:
  **`ap-southeast-2` (Sydney)** — has the `g5`/`g6` families. Set it as your
  default in `aws configure`.
- Create a bucket (globally-unique name), same region:
  ```bash
  aws s3 mb s3://my-anuga-bucket --region ap-southeast-2
  ```

## 8. Dry-run, then run a GPU job
Check everything (creds, GPU quota, AMI, the launch plan) **without spending**:
```bash
docker/aws_run_gpu.sh --region ap-southeast-2 \
  --output s3://my-anuga-bucket/anuga/out \
  --command "python run.py" --dry-run
```
Then launch for real:
```bash
docker/aws_run_gpu.sh --region ap-southeast-2 \
  --upload  ./my_project \
  --input   s3://my-anuga-bucket/anuga/in \
  --output  s3://my-anuga-bucket/anuga/out \
  --command "python run_Tohoku.py -alg DE0"
```
It launches one `g5.2xlarge` (A10G / cc86, 8 vCPU, 32 GB), runs the container,
uploads results (+ `anuga-run.log`) to `--output`, and **self-terminates**.
A ~1M-triangle single-GPU run fits comfortably (a few GB of GPU memory).

## 9. Cost hygiene
- The instance self-terminates, but verify occasionally:
  ```bash
  aws ec2 describe-instances --region ap-southeast-2 \
    --filters Name=tag:Name,Values=anuga-gpu-run \
              Name=instance-state-name,Values=running \
    --query 'Reservations[].Instances[].InstanceId' --output text
  ```
- Delete outputs you don't need (`aws s3 rm --recursive s3://.../out`) or set an
  S3 lifecycle rule. Your Budget alert (step 5) is the backstop.

---

**Do first / in parallel:** step 5 (budget) and step 6 (quota) — the quota
approval is the long pole, and the budget protects you from surprises.
