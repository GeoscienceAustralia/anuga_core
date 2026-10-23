# Validation reports for each release

Decided 2026-10-01. Every release ships its validation reports. They are built
from the tagged commit for DE0, DE1 and DE_ader2, attached to the GitHub
release, and, for final releases, archived on Zenodo as versions of one
record, *ANUGA validation reports*, filed in the ANUGA Zenodo community.

Zenodo accounts are personal: there is no organisation login. A maintainer
logs in to Zenodo with their own GitHub account, and the API tokens come from
that account. The Zenodo **community** is the shared home: its owners and
managers are the ANUGA maintainers, so the collection does not depend on one
person. Share each record with at least one other maintainer, with manage
rights, so someone else can publish a version if needed.

Workflow: `.github/workflows/validation-report.yml`.
Zenodo upload: `.github/scripts/zenodo_validation_reports.py`.

## What runs

Publishing a GitHub release (the same event that triggers PyPI) starts it:

1. **report**: one job per algorithm (DE0, DE1, DE_ader2), run in parallel.
   Each runs `validations_produce_results.py -alg <alg> -l` (every case's
   `produce_results.py`, then the LaTeX), and keeps the PDF and logs as an
   artifact for 90 days. It is a gate: the script exits non-zero if any case
   fails or the report does not typeset, and then nothing below runs.
2. **combine**: builds `anuga-<tag>-validation-combined.pdf` from the
   figures each report job collected (`validations_combine_algs.py`): every
   case once, a summary table, and each figure that differs between the
   algorithms as one panel per algorithm. It needs only LaTeX.
3. **attach**: uploads the combined PDF and `anuga-<tag>-validation-<alg>.pdf`
   to the release. This runs for release candidates too.
4. **zenodo**: final releases only. Creates a **draft** new version of the
   Zenodo record with the four PDFs and the metadata: version = tag, creators
   from `CITATION.cff`, CC-BY 4.0, *isSupplementTo* the GitHub release. It
   never publishes, because a published Zenodo record cannot be deleted.

## One-off setup

- [x] Create the ANUGA Zenodo community, identifier `anuga`
      (https://zenodo.org/communities/anuga), and store the identifier as the
      repository variable `ZENODO_COMMUNITY` (both done 2026-10-01). Drafts
      are then filed into it, and publishing one submits it to the
      community, where an owner accepts it.
- [ ] Add the other maintainers to the community as managers.
- [ ] On zenodo.org, logged in with your own GitHub account:
      Applications → Personal access tokens → new token, scopes
      `deposit:write` and `deposit:actions`. Store it as the repository
      secret `ZENODO_TOKEN` (Settings → Secrets and variables → Actions).
- [ ] The same on sandbox.zenodo.org (a separate site with separate
      accounts and tokens), stored as `ZENODO_SANDBOX_TOKEN`. To test the
      community step too, create a community there and store its identifier
      as `ZENODO_SANDBOX_COMMUNITY`; leave it unset to skip that step.
- [ ] Optional but recommended: in Zenodo's GitHub settings (your account;
      owning anuga-community on GitHub is what lets you enable the repo),
      switch on `anuga-community/anuga_core`. Each release is then also
      archived as **software** with its own DOI, a different record from the
      reports. Once it exists, add the concept DOI to `docs/source/citing.rst`
      and `CITATION.cff` (`doi:`).
- [x] Dry run against the sandbox (2026-10-02). It was run locally with
      `zenodo_validation_reports.py --sandbox`, because "Run workflow" only
      offers workflows that are on the default branch (`main`), which this
      one will not be until the next release. Both paths worked against the
      real API: the first run created and published record 612560 (concept
      612559, metadata and files checked), and the second, with
      `--concept-recid 612559`, made a new-version draft, removed the
      carried-over files and uploaded the new ones. `ZENODO_SANDBOX_CONCEPT_RECID`
      is set to 612559. Once the workflow is on `main`, a by-hand run with
      zenodo = sandbox tests the workflow wiring itself (secrets, artifacts).
      Publishing on the sandbox is harmless.

The first production run has no `ZENODO_CONCEPT_RECID` variable, so it
creates the record. Publish that draft, then store its concept record id as
the variable `ZENODO_CONCEPT_RECID`. Every later release then becomes a new
version of the same record, under the same concept DOI.

## Each release

- [ ] After publishing the GitHub release, watch *Validation reports* in the
      Actions tab (about 2-3 h; the three algorithms run in parallel).
- [ ] If a report job fails, its artifact holds `produce_<alg>.log` and the
      LaTeX log. Fix the case, then re-run by hand with tag = the release.
- [ ] Final releases: open the Zenodo draft linked from the zenodo job's
      summary, check the four PDFs open and the metadata is right, and press
      **Publish**.
- [ ] Add the new version's DOI to the GitHub release notes.

## Re-running by hand

Actions → Validation reports → Run workflow takes the tag, whether to attach
to the GitHub release (`--clobber` replaces existing assets), and `none`,
`sandbox` or `production` for Zenodo.

If a draft new version is already open on Zenodo, a new version cannot be
created: the zenodo job fails with Zenodo's message. Publish or discard the
existing draft first.
