#!/usr/bin/python

import sys

probesets=set()
probesets.add("ProbeSetName")
probesets.add("probeset_id")
probesets.add("id")
probesets.add("GENERIC")
probesets.add("GENERIC:1")
probesets.add("AX-1")

f=open(sys.argv[1],'r')
for line in f:
    probesets.add(line.rstrip())
    probesets.add(line.rstrip()+':1')
    probesets.add(line.rstrip()+'-A')
    probesets.add(line.rstrip()+'-B')
f.close()

f=open(sys.argv[2],'r')
for line in f:
    t=line.rstrip().split('\t',2)
    if t[0] in probesets or t[0].startswith('#') or t[1] in probesets:
        print line.rstrip()
f.close()
