"""Give notebook animations a still image for the PDF build.

The example notebooks store their animations as ``text/html`` output (the
matplotlib JavaScript player), which LaTeX cannot show: the PDF had only the
``<matplotlib.animation.FuncAnimation ...>`` text. The player holds every
frame as a base64 image, so this takes a few of them, evenly spaced from first
to last, and adds them to the output as one ``image/png`` contact sheet.
nbsphinx keeps the player for HTML and uses the image for LaTeX.

Nothing is written to the notebooks: the image is added as the source is read.
"""

import base64
import io
import json
import re

FRAME = re.compile(r'data:image/(?:png|jpeg);base64,([A-Za-z0-9+/=\\n\s]+?)"')

COLUMNS = 2
ROWS = 3
GAP = 12          # pixels between frames
FRAME_WIDTH = 800  # each frame is scaled to this width


def contact_sheet(frames):
    """One PNG (bytes) of up to ROWS x COLUMNS frames, first to last."""

    from PIL import Image

    number = len(frames)
    wanted = min(ROWS * COLUMNS, number)
    if wanted == 1:
        picks = [0]
    else:
        picks = sorted({round(i * (number - 1) / (wanted - 1))
                        for i in range(wanted)})

    images = []
    for k in picks:
        raw = base64.b64decode(frames[k].replace('\\n', ''))
        image = Image.open(io.BytesIO(raw)).convert('RGB')
        height = round(image.height * FRAME_WIDTH / image.width)
        images.append(image.resize((FRAME_WIDTH, height), Image.LANCZOS))

    columns = min(COLUMNS, len(images))
    rows = -(-len(images) // columns)
    height = max(image.height for image in images)
    sheet = Image.new('RGB', (columns * FRAME_WIDTH + (columns - 1) * GAP,
                              rows * height + (rows - 1) * GAP), 'white')
    for i, image in enumerate(images):
        x = (i % columns) * (FRAME_WIDTH + GAP)
        y = (i // columns) * (height + GAP)
        sheet.paste(image, (x, y))

    out = io.BytesIO()
    sheet.save(out, format='PNG', optimize=True)
    return out.getvalue()


def add_stills(app, docname, source):
    """source-read handler: source[0] is the notebook's JSON for an .ipynb."""

    text = source[0]
    if 'FuncAnimation' not in text or not text.lstrip().startswith('{'):
        return
    try:
        notebook = json.loads(text)
    except ValueError:
        return

    changed = False
    for cell in notebook.get('cells', []):
        for output in cell.get('outputs', []):
            data = output.get('data', {})
            if 'text/html' not in data or 'image/png' in data:
                continue
            frames = FRAME.findall(''.join(data['text/html']))
            if not frames:
                continue
            sheet = contact_sheet(frames)
            data['image/png'] = base64.b64encode(sheet).decode('ascii')
            changed = True

    if changed:
        source[0] = json.dumps(notebook)


def setup(app):
    app.connect('source-read', add_stills)
    return {'parallel_read_safe': True, 'parallel_write_safe': True}
