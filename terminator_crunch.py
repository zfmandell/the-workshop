#!/usr/bin/env python

from __future__ import division
from glob import glob
import sys
from collections import Counter

cipher = {'Corynebacterium_SV.bag' : 61.7,'Kurthia_NT.bag' : 36.9, 'Enterococcus_NT.bag' : 37.4, 'Chlamydia_NT.bag' : 41.3, 'Renibacterium_HT.bag' : 56.3, 'Intrasporangium_HT.bag' : 70.75, 'Kitasatospora_NT.bag' : 72, 'Vibrio_SV.bag' : 45.5503, 'Dickeya_SV.bag' : 56.2,'Edwardsiella_SV.bag' : 57.4, 'Cellulosimicrobium_HT.bag' : 74.7}


def median(lst):
    n = len(lst)
    s = sorted(lst)
    return (sum(s[n//2-1:n//2+1])/2.0, s[n//2])[n % 2] if n else None

for fyle in glob('*.bag'):
    results = []
    with open(fyle,'r') as inp:
        for line in inp:
            if len(line.strip().split()) == 14:
                coun = Counter(line.strip().split()[8])
                results.append(((coun['G']+coun['C'])/len(line.strip().split()[8]))*100)

    print fyle
    print median(results)-cipher[fyle]
    print median(results)
