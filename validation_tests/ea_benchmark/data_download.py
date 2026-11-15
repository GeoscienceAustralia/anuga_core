"""
Download the Environment Agency 2D benchmark input data (SC120002,
Neelz and Pender 2013) and unpack each test into data/TestN/.

The data are published under the Open Government Licence v3.0:
https://gov.uk/flood-and-coastal-erosion-risk-management-research-reports/2d-benchmarking-evaluating-the-latest-generation-of-the-hydraulic-models-for-fcrm

Usage: python data_download.py [--force]
"""

import hashlib
import io
import os
import shutil
import sys
import urllib.request
import zipfile

URL = ('https://assets.publishing.service.gov.uk/media/'
       '6033a9bbd3bf7f72182e99ca/Benchmarking_Model_Data.zip')
SHA256 = '388c2789570e5975fe909033e34bedbe5d7bde5e6d7f9b7d7c1eaf04e58bc732'

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, 'data')
ZIP_FILE = os.path.join(HERE, 'Benchmarking_Model_Data.zip')


def download(force=False):
    if os.path.exists(ZIP_FILE) and not force:
        print('Using existing', ZIP_FILE)
    else:
        print('Downloading', URL)
        tmp = ZIP_FILE + '.part'
        with urllib.request.urlopen(URL) as response, open(tmp, 'wb') as f:
            shutil.copyfileobj(response, f)
        os.replace(tmp, ZIP_FILE)

    with open(ZIP_FILE, 'rb') as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    if digest != SHA256:
        raise RuntimeError('Checksum mismatch for %s: got %s, expected %s'
                           % (ZIP_FILE, digest, SHA256))


def _test_dir(zip_name):
    """'Test1_dataset_2010.zip' -> 'Test1'; 'XS info Tests_5&6.zip' -> 'XS_Tests5_6'."""
    if zip_name.startswith('Test'):
        return zip_name.split('_')[0]
    return 'XS_Tests5_6'


def _extract_flat(zf, dest):
    """Extract every file into dest, dropping the archive's own folders.

    The per-test archives are inconsistent: most wrap their files in a
    'TestN dataset 2010/' folder, Test 7's has none.
    """
    os.makedirs(dest, exist_ok=True)
    for info in zf.infolist():
        if info.is_dir():
            continue
        name = os.path.basename(info.filename)
        with zf.open(info) as src, open(os.path.join(dest, name), 'wb') as f:
            shutil.copyfileobj(src, f)


def unpack():
    """The outer ZIP holds one ZIP per test; unpack each into data/TestN/."""
    if os.path.isdir(DATA_DIR):
        shutil.rmtree(DATA_DIR)
    with zipfile.ZipFile(ZIP_FILE) as outer:
        for name in outer.namelist():
            if not name.endswith('.zip'):
                continue
            dest = os.path.join(DATA_DIR, _test_dir(name))
            print('Unpacking', name, '->', os.path.relpath(dest, HERE))
            with zipfile.ZipFile(io.BytesIO(outer.read(name))) as inner:
                _extract_flat(inner, dest)


if __name__ == '__main__':
    download(force='--force' in sys.argv)
    unpack()
    print('Data in', DATA_DIR)
