"""Shrink the images stored in notebook outputs without re-running them.

The frames of matplotlib's jshtml animations are re-encoded as JPEG (quality
85, no chroma subsampling): they are mostly smooth hillshade, where a 256-colour
PNG bands visibly. Static image/png figures are recompressed losslessly.
Each image is replaced only if the result is smaller.

Run it after re-running and saving a notebook in examples/notebooks: the
docs render those stored outputs, and the animations otherwise store every
frame as a full-colour PNG (the five notebooks went from 21 to 12 MB).

Usage: python tools/slim_notebook_outputs.py [--dry-run] NOTEBOOK...
"""
import base64, io, json, re, sys
from PIL import Image

# matplotlib's jshtml breaks each frame's base64 with backslash-newline
FRAME = re.compile(r'data:image/png;base64,([A-Za-z0-9+/=\\\n]+)')


def _decode(b64):
    return base64.b64decode(b64.replace('\\\n', '').replace('\n', ''))


def frame(m):
    raw = _decode(m.group(1))
    im = Image.open(io.BytesIO(raw))
    if 'A' in im.getbands():                      # flatten onto white
        bg = Image.new('RGB', im.size, 'white'); bg.paste(im, mask=im.getchannel('A')); im = bg
    out = io.BytesIO(); im.convert('RGB').save(out, 'JPEG', quality=85, optimize=True, subsampling=0)
    if out.tell() >= len(raw):
        return m.group(0)
    return 'data:image/jpeg;base64,' + base64.b64encode(out.getvalue()).decode()


def static(b64):
    raw = _decode(b64)
    out = io.BytesIO(); Image.open(io.BytesIO(raw)).save(out, 'PNG', optimize=True)
    return base64.b64encode(out.getvalue()).decode() if out.tell() < len(raw) else b64.replace('\n', '')


def main():
    dry = '--dry-run' in sys.argv
    total = [0, 0]
    for path in [a for a in sys.argv[1:] if a != '--dry-run']:
        raw = open(path, newline='').read()
        nb = json.loads(raw)
        for c in nb['cells']:
            for o in c.get('outputs', []):
                d = o.get('data', {})
                if 'image/png' in d:
                    v = d['image/png']; s = v if isinstance(v, str) else ''.join(v)
                    d['image/png'] = static(s) + ('\n' if s.endswith('\n') else '')
                if 'text/html' in d:
                    v = d['text/html']; html = v if isinstance(v, str) else ''.join(v)
                    d['text/html'] = FRAME.sub(frame, html).splitlines(keepends=True)
        text = json.dumps(nb, indent=1, ensure_ascii=False) + '\n'
        total[0] += len(raw); total[1] += len(text)
        print('%-46s %5.2f MB -> %5.2f MB' % (path.split('/')[-1], len(raw) / 1e6, len(text) / 1e6))
        if not dry:
            open(path, 'w', newline='').write(text)
    print('%-46s %5.2f MB -> %5.2f MB' % ('total', total[0] / 1e6, total[1] / 1e6))


if __name__ == '__main__':
    main()
