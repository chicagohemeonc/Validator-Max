#!/usr/bin/env python

from PeptideFragmentSingleton import PeptideFragment

f = PeptideFragment( 1 ) # mass type = mono

#                   charge   modifications
f.analyze('ACKRM', 3,       "16.0@M 57.1@3")


print f.sequence(), ": ", f.composition();
print "---- pI: ---- "
print "  ", f.pI()

print "---- b-ions: ---- "
l = f.b_ions();
for i, m in enumerate(l):
    print "  ", m

print "---- y-ions: ---- "
l = f.y_ions();
for i, m in enumerate(l):
    print "  ", m

print "---- mass: ---- "
l = f.peptide_mass();
for i, m in enumerate(l):
    print "  ", m

print

f.analyze('KC')
print f.sequence(), ": ", f.composition();

print "---- pI: ---- "
print "  ", f.pI()

print "---- b-ions: ---- "
l = f.b_ions();
for i, m in enumerate(l):
    print "  ", m

print "---- y-ions: ---- "
l = f.y_ions();
for i, m in enumerate(l):
    print "  ", m

print "---- mass: ---- "
l = f.peptide_mass();
for i, m in enumerate(l):
    print "  ", m
