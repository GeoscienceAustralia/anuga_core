# anugaSed comparison harness

The runnable sediment examples moved to **`examples/sediment/`**. What is left
here is not an example: `compare_with_anugased.py` runs **anugaSed** and this
module on the same case in one process and diffs them, so it needs anugaSed
installed (the py3-modernisation branch, `pip install -e .`) and will not run
without it.

It is kept out of `examples/` for that reason -- it documents and sizes the
known divergences (D1a, D4, D5, D6) rather than showing anyone how to use the
module.
