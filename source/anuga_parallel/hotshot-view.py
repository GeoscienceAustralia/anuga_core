#!/usr/bin/env python
import sys
import hotshot.stats

if len(sys.argv) == 2:
    stats = hotshot.stats.load(sys.argv[1])
    stats.strip_dirs().sort_stats('cumulative').print_stats(20)
else:
    print "Usage: %s profile_filename" % sys.argv[0]
