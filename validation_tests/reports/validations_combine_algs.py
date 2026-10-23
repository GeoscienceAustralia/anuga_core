"""
Build the combined validation report: every case once, with each of its
figures shown for DE0, DE1 and DE_ader2 side by side.

Two steps:

    python validations_combine_algs.py collect -alg DE1
        Run after validations_produce_results.py -alg DE1. Copies the figures
        each case's results.tex uses, the case status file and
        saved_parameters.tex into algs/DE1/, so the next algorithm's run can
        overwrite the case directories.

    python validations_combine_algs.py combine
        Once algs/<alg>/ exists for every algorithm. A figure that is identical
        for all of them (a setup diagram, a mesh) goes to combined/shared/ and
        is shown once; any other goes to combined/<alg>/ and is shown as one
        panel per algorithm. Writes the summary table and typesets
        validations_report_combined.pdf.

The release workflow runs collect in each algorithm's job and combine in a
job of its own, so combine needs only Python and LaTeX, not anuga.
"""
import argparse
import filecmp
import os
import re
import shutil
import subprocess
import sys

ALGS = ['DE0', 'DE1', 'DE_ader2']
HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.dirname(HERE)          # validation_tests/

INCLUDE = re.compile(r'\\includegraphics\s*(?:\[[^\]]*\])?\s*\{([^}]*)\}')
FIGIF = re.compile(r'\\figifexists\s*\{([^}]*)\}')
CASE = re.compile(r'\\inputresults\s*\{\.\./([^}]*)\}')


def uncommented(text):
    """The text with LaTeX comments removed (keeping escaped \\%)."""
    return '\n'.join(re.sub(r'(?<!\\)%.*', '', line) for line in text.splitlines())


def cases():
    """Case paths relative to validation_tests/, in report order."""
    with open(os.path.join(HERE, 'validations_report_body.tex')) as f:
        text = uncommented(f.read())
    # The appendix shows \inputresults{../Directory/Name} as an example
    text = re.sub(r'\\begin\{verbatim\}.*?\\end\{verbatim\}', '', text, flags=re.S)
    return CASE.findall(text)


def figures(case):
    path = os.path.join(TESTS, case, 'results.tex')
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        text = uncommented(f.read())
    names = INCLUDE.findall(text) + FIGIF.findall(text)
    return sorted(set(n.strip() for n in names))


def collect(alg):
    dest = os.path.join(HERE, 'algs', alg)
    if os.path.isdir(dest):
        shutil.rmtree(dest)
    os.makedirs(dest)
    n = 0
    for case in cases():
        for name in figures(case):
            src = os.path.join(TESTS, case, name)
            if os.path.isfile(src):
                os.makedirs(os.path.dirname(os.path.join(dest, case, name)), exist_ok=True)
                shutil.copy2(src, os.path.join(dest, case, name))
                n += 1
    for name in ('case_status_alg_%s.csv' % alg, 'saved_parameters.tex'):
        if os.path.isfile(os.path.join(HERE, name)):
            shutil.copy2(os.path.join(HERE, name), dest)
    print('Collected %d figures for %s into %s' % (n, alg, dest))


def read_status(alg):
    path = os.path.join(HERE, 'algs', alg, 'case_status_alg_%s.csv' % alg)
    status = {}
    if os.path.isfile(path):
        with open(path) as f:
            next(f)
            for line in f:
                case, outcome, secs = line.strip().split(',')
                status[case] = (outcome, float(secs))
    return status


def tex_escape(s):
    return s.replace('\\', r'\textbackslash{}').replace('_', r'\_').replace('&', r'\&')


def summary_table(algs):
    status = {alg: read_status(alg) for alg in algs}
    lines = [
        r'\begin{longtable}{l' + 'r' * len(algs) + '}',
        r'\hline',
        'Case & ' + ' & '.join(tex_escape(a) for a in algs) + r' \\',
        r'\hline\endhead',
    ]
    chapter = None
    for case in cases():
        top, name = case.split('/', 1)
        if top != chapter:
            chapter = top
            lines.append(r'\multicolumn{%d}{l}{\textit{%s}} \\' % (len(algs) + 1, tex_escape(top)))
        cells = []
        for alg in algs:
            outcome, secs = status[alg].get(case, ('missing', None))
            if outcome == 'ok':
                cells.append('%.0f s' % secs)
            else:
                cells.append(r'\textbf{%s}' % outcome)
        lines.append(r'\quad %s & %s \\' % (tex_escape(name), ' & '.join(cells)))
    lines += [r'\hline', r'\end{longtable}']
    return '\n'.join(lines) + '\n'


def algs_macros(algs):
    """The macros validations_report_combined.tex lays the panels out with."""
    width = '%.3f' % (0.99/len(algs))
    names = [tex_escape(a) for a in algs]
    listed = names[0] if len(names) == 1 else ', '.join(names[:-1]) + ' and ' + names[-1]
    panels = r'\hfill'.join(r'\algpanel{%s}{%s}{%s}{#1}' % (width, a, n)
                             for a, n in zip(algs, names))
    exists = '#2'
    for a in reversed(algs):
        exists = r'\IfFileExists{combined/%s/\casepath/#1}{#2}{%s}' % (a, exists if a != algs[-1] else '#3')
    return '\n'.join([
        r'\newcommand{\algnames}{%s}' % listed,
        r'\newcommand{\allalgpanels}[1]{\par\noindent%s\par}' % panels,
        r'\newcommand{\anyalgfigexists}[3]{\IfFileExists{combined/shared/\casepath/#1}{#2}{%s}}' % exists,
    ]) + '\n'


def combine(algs):
    out = os.path.join(HERE, 'combined')
    if os.path.isdir(out):
        shutil.rmtree(out)
    shared = varying = 0
    for case in cases():
        for name in figures(case):
            srcs = {alg: os.path.join(HERE, 'algs', alg, case, name) for alg in algs}
            have = [a for a in algs if os.path.isfile(srcs[a])]
            if not have:
                continue
            same = len(have) == len(algs) and all(
                filecmp.cmp(srcs[have[0]], srcs[a], shallow=False) for a in have[1:])
            targets = ['shared'] if same else have
            for t in targets:
                dst = os.path.join(out, t, case, name)
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                shutil.copy2(srcs[have[0]] if same else srcs[t], dst)
            shared += same
            varying += not same
    print('%d figures shared by all algorithms, %d shown per algorithm' % (shared, varying))

    with open(os.path.join(HERE, 'combined_summary.tex'), 'w') as f:
        f.write(summary_table(algs))

    with open(os.path.join(HERE, 'combined_algs.tex'), 'w') as f:
        f.write(algs_macros(algs))

    # Version and time from any algorithm's run; the algorithm list replaces \alg
    for alg in algs:
        p = os.path.join(HERE, 'algs', alg, 'saved_parameters.tex')
        if os.path.isfile(p):
            shutil.copy2(p, os.path.join(HERE, 'saved_parameters.tex'))
            break

    pdf = os.path.join(HERE, 'validations_report_combined.pdf')
    if os.path.isfile(pdf):
        os.remove(pdf)
    latex = ['pdflatex', '-interaction=batchmode', 'validations_report_combined.tex']
    for cmd in (latex, ['bibtex', 'validations_report_combined'], latex, latex):
        subprocess.call(cmd, cwd=HERE)
    if not os.path.isfile(pdf):
        sys.exit('validations_report_combined.pdf was not produced '
                 '(see validations_report_combined.log)')
    print('Wrote ' + pdf)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    sub = p.add_subparsers(dest='step', required=True)
    c = sub.add_parser('collect')
    c.add_argument('-alg', required=True)
    m = sub.add_parser('combine')
    m.add_argument('-algs', default=','.join(ALGS),
                   help='comma-separated algorithms, in column order')
    a = p.parse_args()
    if a.step == 'collect':
        collect(a.alg)
    else:
        combine(a.algs.split(','))


if __name__ == '__main__':
    main()
