The module here (order_boundary.py) contains a function that takes a CSV
file of boundary points and orders the points into an output CSV file.

The first point in the file is assumed to be the start point.  All other
points in the file are the points that are closest to each succeeding point.

Note that in pathological cases the ordering may be wrong.  Always check!
CAVEAT EMPTOR.
