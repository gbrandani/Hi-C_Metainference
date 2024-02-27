#!/usr/bin/env python

import sys
import math
import numpy as np
import MDAnalysis

# define atoms class
class Atom:
    pass
class Bond:
    pass
class Angle:
    pass
class Dihedral:
    pass

def splitp(l):
    f,i=float,int
    return [l[0:6],i(l[6:11]),l[12:16],l[16],l[17:20],l[21],i(l[22:26]),l[26],f(l[30:38]),f(l[38:46]),f(l[46:54]),f(l[54:60]),f(l[60:66])]
def writeatom(a):
    return '%6s%5d %4s %3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n'%('ATOM  ',a.id,a.type,a.res,a.rid,a.pos[0],a.pos[1],a.pos[2],1.0,1.0)
    '''
Splits a line from a pdb file in the correct fields, and convert Angstrom to nm.
 1 -  6        Record name     "ATOM  "                                            
 7 - 11        Integer         Atom serial number.                   
13 - 16        Atom            Atom name.                            
17             Character       Alternate location indicator.         
18 - 20        Residue name    Residue name.                         
22             Character       Chain identifier.                     
23 - 26        Integer         Residue sequence number.              
27             AChar           Code for insertion of residues.       
31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.                       
39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.                            
47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.                            
55 - 60        Real(6.2)       Occupancy.                            
61 - 66        Real(6.2)       Temperature factor (Default = 0.0).                   
73 - 76        LString(4)      Segment identifier, left-justified.   
77 - 78        LString(2)      Element symbol, right-justified.      
79 - 80        LString(2)      Charge on the atom.
    '''

# center of mass
def com(atoms):
  dim=len(atoms[0].pos)
  cm=np.array([0. for i in range(dim)])
  M=0.
  for a in atoms:
    cm+=1.0*a.pos
    M+=1.0
  return cm/M

# input
print('# syntax: make_conf.py -ncg 200 -in input.pdb -o conf.pdb -d data.input')
ncg=30 # number of beads
input_fname=''
conf_fname ='conf.pdb'
data_fname='data.input'
for i in range(len(sys.argv)):
  if(sys.argv[i]=='-ncg'): ncg=int(sys.argv[i+1])
  elif(sys.argv[i]=='-in'): input_fname=sys.argv[i+1]
  elif(sys.argv[i]=='-o'):  conf_fname =sys.argv[i+1]
  elif(sys.argv[i]=='-d'):  data_fname =sys.argv[i+1]

# if argument is given, take positions from PDB
conf=[]
if(input_fname!=''):
  u = MDAnalysis.Universe(input_fname,input_fname)
  for a in u.atoms:
    conf.append(a.position)
  ncg=len(conf)

# create input structure
sigma=1.0
print('# creating structure ..')
atoms=[]
for i in range(ncg):
  atom=Atom()
  atom.id=i+1
  atom.type='H'
  atom.res='ALA'
  atom.rid=i+1
  if(i==0):
    atom.pos=np.array([0.,0.,0.])
  else:
    if(i<ncg//4):
      dx=np.array([0.,0.,sigma])
    elif(i<2*ncg//4):
      dx=np.array([0.,sigma,0.])
    elif(i<3*ncg//4+1):
      dx=np.array([0.,0.,-sigma])
    else:
      dx=np.array([0.,-sigma,0.])
    atom.pos=atoms[-1].pos+dx
  if(len(conf)>0):
    atom.pos=conf[i]
  atoms.append( atom )
print('# %d atoms found'%len(atoms))

# set center of mass to zero and print pdb
outpdb=open(conf_fname,'w')
center=com(atoms)
for a in atoms: a.pos=a.pos-center
for a in atoms: outpdb.write(writeatom(a))
outpdb.close()

# define bonds
bonds=[]
angles=[]
dihedrals=[]
# make bonds
nbonds=0
nangles=0
for i in range(ncg):
  if(i>0):
    b=Bond()
    nbonds+=1
    b.id=nbonds
    b.type=1
    b.a1=i
    b.a2=i+1
    bonds.append(b)
  if(i>1):
    b=Angle()
    nangles+=1
    b.id=nangles
    b.type=1
    b.a1=i-1
    b.a2=i
    b.a3=i+1
    angles.append(b)
# check
if(len(bonds)!=(ncg-1)):  print('# nbonds!=ncg-1'); sys.exit()
else: print('#',len(bonds),'bonds created')
if(len(angles)!=(ncg-2)): print('# nangles!=ncg-2'); sys.exit()
else: print('#',len(angles),'bonds created')

# print data file
for a in atoms: a.mol=1; a.ix=0; a.iy=0; a.iz=0; a.type=1
print('# printing data file..')
outdata=open(data_fname,'w')
outdata.write('LAMMPS data file\n\n')
outdata.write('%d atoms\n'%len(atoms))
if(len(bonds)>0): outdata.write('%d bonds\n'%len(bonds))
if(len(angles)>0): outdata.write('%d angles\n'%len(angles))
outdata.write('\n')
outdata.write('%d atom types\n'%len(set(a.type for a in atoms)))
if(len(bonds)>0): outdata.write('%d bond types\n'%len(set(b.type for b in bonds)))
if(len(angles)>0): outdata.write('%d angle types\n'%len(set(a.type for a in angles)))
outdata.write('\n')
halfbox=ncg/8.+20.*sigma # enough to fit entire protein fully stretched, ~250/2
outdata.write('%f %f xlo xhi\n'%(-halfbox, halfbox))
outdata.write('%f %f ylo yhi\n'%(-halfbox, halfbox))
outdata.write('%f %f zlo zhi\n'%(-halfbox, halfbox))
outdata.write('\nMasses\n\n')
# the mass of carbon is 12 grams/mole, let's set all masses to 1, as for hydrogen
for atype in sorted(set(a.type for a in atoms)): outdata.write('%d 1.0\n'%atype)

# print atoms
outdata.write('\nAtoms\n\n')
for a in atoms: outdata.write('%d %d %d %f %f %f %d %d %d\n'%(a.id,a.mol,a.type,a.pos[0],a.pos[1],a.pos[2],a.ix,a.iy,a.iz))
# print bonds
if(len(bonds)>0):
  outdata.write('\nBonds\n\n')
  for b in bonds: outdata.write('%d %d %d %d\n'%(b.id,b.type,b.a1,b.a2))
# print angles
if(len(angles)>0):
  outdata.write('\nAngles\n\n')
  for a in angles: outdata.write('%d %d %d %d %d\n'%(a.id,a.type,a.a1,a.a2,a.a3))
outdata.close()

