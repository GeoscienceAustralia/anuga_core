"""Upload the validation reports of an ANUGA release to Zenodo as a DRAFT.

The reports live in their own Zenodo record ("ANUGA validation reports"), one
version per ANUGA release, linked to the release it validates. This script
creates the next version of that record (or, given no concept record id, the
record itself), replaces its files with the given PDFs, sets the metadata,
and stops: it never publishes. A published Zenodo record cannot be deleted,
so a person looks at the draft and presses Publish (see
claude/RELEASE_VALIDATION_REPORTS.md).

    python zenodo_validation_reports.py --tag 4.1.0 \
        --release-url https://github.com/anuga-community/anuga_core/releases/tag/4.1.0 \
        [--concept-recid 1234567] [--community anuga] [--sandbox] [--dry-run] \
        report1.pdf report2.pdf ...

The token comes from the ZENODO_TOKEN environment variable (scope
deposit:write); --sandbox talks to sandbox.zenodo.org, which needs its own
token. Creators are read from CITATION.cff.
"""
import argparse
import datetime
import os
import sys

import requests
import yaml

PRODUCTION = 'https://zenodo.org/api'
SANDBOX = 'https://sandbox.zenodo.org/api'


def creators_from_citation(path):
    with open(path) as f:
        cff = yaml.safe_load(f)
    creators = []
    for a in cff.get('authors', []):
        c = {'name': '%s, %s' % (a['family-names'], a['given-names'])}
        if a.get('affiliation'):
            c['affiliation'] = a['affiliation']
        if a.get('orcid'):
            c['orcid'] = a['orcid'].rsplit('/', 1)[-1]
        creators.append(c)
    return creators


def metadata(tag, release_url, creators, community=''):
    meta = {
        'title': 'ANUGA validation reports',
        'upload_type': 'publication',
        'publication_type': 'report',
        'version': tag,
        'publication_date': datetime.date.today().isoformat(),
        'description': (
            '<p>Validation reports for ANUGA %s, built from the tagged release '
            'by its release workflow. They cover the analytical, '
            'experimental, behaviour and sediment cases of the ANUGA '
            'validation suite; the large case studies are not included.</p>'
            '<p>The combined report (file name ending -combined.pdf) runs every '
            'case with the flow algorithms DE0, DE1 and DE_ader2 and shows '
            'their results side by side, with a summary table at the front. '
            'The other files are the full report for each algorithm on its '
            'own, with full-size figures.</p>'
            '<p>The cases, their reference solutions and the scripts that '
            'produce these reports are in the validation_tests directory of '
            'the ANUGA repository.</p>' % tag),
        'creators': creators,
        'license': 'cc-by-4.0',
        'keywords': ['ANUGA', 'shallow water equations', 'validation',
                     'finite volume method', 'flooding', 'tsunami'],
        'related_identifiers': [{
            'identifier': release_url,
            'relation': 'isSupplementTo',
            'resource_type': 'software',
        }],
    }
    if community:
        # Publishing the draft then also submits it to the community
        meta['communities'] = [{'identifier': community}]
    return meta


class Zenodo:
    def __init__(self, base, token, dry_run=False):
        self.base = base
        self.dry_run = dry_run
        self.session = requests.Session()
        self.session.headers['Authorization'] = 'Bearer %s' % token

    def call(self, method, url, **kw):
        if not url.startswith('http'):
            url = self.base + url
        r = self.session.request(method, url, timeout=300, **kw)
        if r.status_code >= 400:
            sys.exit('%s %s failed (%d): %s' % (method, url, r.status_code, r.text[:1000]))
        return r.json() if r.content else {}

    def new_draft(self, concept_recid):
        """A new draft: the next version of concept_recid, or a new record."""
        if not concept_recid:
            print('No concept record id: creating a new record')
            return self.call('POST', '/deposit/depositions', json={})
        # The concept id resolves to the latest published version
        latest = self.call('GET', '/records/%s' % concept_recid)
        print('Latest published version: record %s (version %s)'
              % (latest['id'], latest.get('metadata', {}).get('version')))
        version = self.call('POST', '/deposit/depositions/%s/actions/newversion' % latest['id'])
        draft = self.call('GET', version['links']['latest_draft'])
        # A new version starts with the previous version's files
        for f in draft.get('files', []):
            print('Removing %s carried over from the previous version' % f['filename'])
            self.call('DELETE', '/deposit/depositions/%s/files/%s' % (draft['id'], f['id']))
        return draft

    def upload(self, draft, path):
        name = os.path.basename(path)
        print('Uploading %s (%.1f MB)' % (name, os.path.getsize(path)/1e6))
        with open(path, 'rb') as f:
            self.call('PUT', '%s/%s' % (draft['links']['bucket'], name), data=f)

    def set_metadata(self, draft, meta):
        return self.call('PUT', '/deposit/depositions/%s' % draft['id'], json={'metadata': meta})


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('files', nargs='+')
    p.add_argument('--tag', required=True)
    p.add_argument('--release-url', required=True)
    p.add_argument('--concept-recid', default='')
    p.add_argument('--community', default='',
                   help='Zenodo community identifier to file the record in')
    p.add_argument('--sandbox', action='store_true')
    p.add_argument('--citation', default='CITATION.cff')
    p.add_argument('--api', help=argparse.SUPPRESS)   # testing only
    p.add_argument('--dry-run', action='store_true',
                   help='print the metadata and the files, contact nobody')
    a = p.parse_args()

    for path in a.files:
        if not os.path.isfile(path):
            sys.exit('No such file: %s' % path)
    meta = metadata(a.tag, a.release_url, creators_from_citation(a.citation), a.community)

    if a.dry_run:
        import json
        print('Would create a draft on %s%s with metadata:'
              % (SANDBOX if a.sandbox else PRODUCTION,
                 ' as a new version of %s' % a.concept_recid if a.concept_recid else ''))
        print(json.dumps(meta, indent=2))
        print('and files: %s' % ', '.join(os.path.basename(f) for f in a.files))
        return

    token = os.environ.get('ZENODO_TOKEN')
    if not token:
        sys.exit('ZENODO_TOKEN is not set')
    z = Zenodo(a.api or (SANDBOX if a.sandbox else PRODUCTION), token)
    draft = z.new_draft(a.concept_recid)
    for path in a.files:
        z.upload(draft, path)
    draft = z.set_metadata(draft, meta)

    html = draft['links'].get('html', '')
    print()
    print('Draft ready (NOT published): %s' % html)
    print('Concept record id: %s' % draft.get('conceptrecid', '?'))
    summary = os.environ.get('GITHUB_STEP_SUMMARY')
    if summary:
        with open(summary, 'a') as f:
            f.write('## Zenodo draft (not published)\n\n')
            f.write('- Draft: %s\n- Concept record id: `%s`\n' % (html, draft.get('conceptrecid', '?')))
            if not a.concept_recid:
                f.write('\nThis was a new record. Store its concept record id in the '
                        'repository variable so later releases become new versions of it.\n')


if __name__ == '__main__':
    main()
