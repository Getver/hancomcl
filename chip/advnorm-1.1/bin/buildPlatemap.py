#!/usr/bin/python

import sys

f=open(sys.argv[1],'r')
for line in f:
    if line.startswith("cel_files"):
        break;
t=line.rstrip('\n').split('\t')
column=t.index("affymetrix-plate-barcode")
print "cel_files\tplate"
for line in f:
    t=line.rstrip('\n').split('\t')
    print '\t'.join([t[0],t[column]])
f.close()

