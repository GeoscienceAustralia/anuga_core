"""
Terminal viewer for ANUGA TOML scenario files.

ANUGA scenario configs are mostly repetition — a `towradgi.toml` is ~1000 lines
of which 235 are one-per-file `[[initial_conditions.friction]]` / `[[culverts]]`
/ `[[rainfall]]` blocks. This module syntax-highlights a TOML file and, by
default, *collapses* long runs of identical array-of-table blocks down to the
first block plus a "N more" marker, so the structure of a scenario is legible
at a glance.

Public API:
    collapse_toml(text, threshold=6) -> str
    highlight(text, color=True) -> str
    render(text, collapse=True, color=True, threshold=6) -> str
    page(text)

Powers ``anuga_toml_run <config> --dry-run --format text``.
"""

import re

# An array-of-tables header:  [[a.b.c]]   (optionally trailing comment)
_ARRAY_HDR = re.compile(r'^\s*\[\[(?P<name>[^\]]+)\]\]\s*(#.*)?$')
# A plain table header:  [a.b]
_TABLE_HDR = re.compile(r'^\s*\[(?P<name>[^\[\]][^\]]*)\]\s*(#.*)?$')


def _is_comment_or_blank(line):
    s = line.strip()
    return s == '' or s.startswith('#')


def _split_segments(lines):
    """Split lines into (header_name_or_None, [lines]) segments.

    A segment begins at a table / array-of-tables header and runs to the next
    header. Comment/blank lines immediately preceding a header are re-attached
    to that header's segment (so a section comment travels with its block).
    """
    # indices of header lines
    starts = [i for i, ln in enumerate(lines)
              if _ARRAY_HDR.match(ln) or _TABLE_HDR.match(ln)]
    if not starts:
        return [(None, lines)]

    # Pull a header's leading comment/blank run back to that header.
    adj = []
    for s in starts:
        j = s
        while j - 1 >= 0 and _is_comment_or_blank(lines[j - 1]):
            j -= 1
        adj.append(j)
    # keep boundaries strictly increasing and after the preamble
    bounds = []
    prev = 0
    for s, a in zip(starts, adj):
        b = max(a, prev)
        bounds.append((b, s))
        prev = s + 1

    segments = []
    if bounds and bounds[0][0] > 0:
        segments.append((None, lines[:bounds[0][0]]))     # preamble
    for k, (b, s) in enumerate(bounds):
        end = bounds[k + 1][0] if k + 1 < len(bounds) else len(lines)
        m = _ARRAY_HDR.match(lines[s])
        name = f'[[{m.group("name")}]]' if m else None
        segments.append((name, lines[b:end]))
    return segments


def collapse_toml(text, threshold=6):
    """Collapse runs of >``threshold`` identical ``[[name]]`` blocks.

    The first block of a run is shown in full; the remainder are replaced by a
    single comment marker. Plain tables and non-repeated blocks are untouched.
    """
    lines = text.splitlines()
    segments = _split_segments(lines)

    out = []
    i = 0
    while i < len(segments):
        name, body = segments[i]
        if name is None:
            out.extend(body)
            i += 1
            continue
        # gather a run of the same array-of-tables header
        j = i
        while j < len(segments) and segments[j][0] == name:
            j += 1
        run = segments[i:j]
        if len(run) > threshold:
            out.extend(run[0][1])                       # first block in full
            n = len(run) - 1
            out.append(f'# ⋯ {n} more {name} block'
                       f'{"s" if n != 1 else ""} collapsed '
                       f'(use --full to expand) ⋯')
            out.append('')
        else:
            for _, b in run:
                out.extend(b)
        i = j
    return '\n'.join(out).rstrip('\n') + '\n'


def highlight(text, color=True):
    """Syntax-highlight TOML for the terminal via pygments, if available."""
    if not color:
        return text
    try:
        from pygments import highlight as _hl
        from pygments.lexers import TOMLLexer
        from pygments.formatters import Terminal256Formatter
    except Exception:
        return text
    try:
        return _hl(text, TOMLLexer(), Terminal256Formatter(style='monokai'))
    except Exception:
        return text


def render(text, collapse=True, color=True, threshold=6):
    """Return the view-ready string: optionally collapsed, optionally coloured."""
    if collapse:
        text = collapse_toml(text, threshold=threshold)
    return highlight(text, color=color)


def page(text):
    """Send *text* to a pager when stdout is an interactive terminal.

    Honours ``$PAGER`` (default ``less -R`` so ANSI colour survives); writes
    plainly when not a tty or if the pager cannot be launched.
    """
    import sys
    import os
    if not sys.stdout.isatty():
        sys.stdout.write(text)
        return
    pager = os.environ.get('PAGER', 'less -R')
    try:
        import subprocess
        proc = subprocess.Popen(pager, shell=True, stdin=subprocess.PIPE)
        proc.communicate(text.encode('utf-8', 'replace'))
    except (OSError, BrokenPipeError, KeyboardInterrupt):
        try:
            sys.stdout.write(text)
        except BrokenPipeError:
            pass
