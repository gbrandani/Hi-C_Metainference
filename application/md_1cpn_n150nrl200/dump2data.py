#!/usr/bin/env python

import sys

# input
print('# syntax: dump2data.py -s in.lammps -f in.dump -o in.traj.lammps -t 0')
topol_fname=''
traj_fname =''
tframe=0
out_fname='in.traj.lammps'
for i in range(len(sys.argv)):
  if(sys.argv[i]=='-f'):   traj_fname =sys.argv[i+1]
  elif(sys.argv[i]=='-s'): topol_fname=sys.argv[i+1]
  elif(sys.argv[i]=='-o'): out_fname  =sys.argv[i+1]
  elif(sys.argv[i]=='-t'): tframe =int(sys.argv[i+1])
if(traj_fname==''):  print('# error: trajectory file not given'); sys.exit()
if(topol_fname==''): print('# error: topology   file not given'); sys.exit()

# read trajectory
positions={}
orientations={}
item=''
for l in open(traj_fname,'r').readlines():
  t=l.split()
  if(t[0]=='ITEM:'):
      item=t[1]
      if(item=='BOX'): i_box=0
  elif(item=='TIMESTEP'):
    timestep=int(l)
    if(timestep>tframe): break
  elif(item=='ATOMS' and timestep==tframe):
    t=l.split()
    ai=int(t[0])
    atype=int(t[1])
    x,y,z=float(t[2]),float(t[3]),float(t[4])
    o1,o2,o3,o4=float(t[5]),float(t[6]),float(t[7]),float(t[8])
    positions[ai]   =[x,y,z]
    orientations[ai]=[o1,o2,o3,o4]

# replace positions and orientations in topology
fout=open(out_fname,'w')
for l in open(topol_fname,'r').readlines():
  t=l.split()
  if(len(t)==0):
    fout.write(l)
  elif(len(t)==1):
    item=t[0]
    fout.write(l)
  elif(item=='Atoms'):
    ai=int(t[0])
    fout.write('%s'%(t[0]))
    for i in range(1,len(t)):
      if(i==2):   fout.write(' %.6f'%(positions[ai][0]))
      elif(i==3): fout.write(' %.6f'%(positions[ai][1]))
      elif(i==4): fout.write(' %.6f'%(positions[ai][2]))
      else:       fout.write(' %s'%(t[i]))
    fout.write('\n')
  elif(item=='Ellipsoids'):
    ai=int(t[0])
    fout.write('%s'%(t[0]))
    for i in range(1,len(t)):
      if(i==4):   fout.write(' %.6f'%(orientations[ai][0]))
      elif(i==5): fout.write(' %.6f'%(orientations[ai][1]))
      elif(i==6): fout.write(' %.6f'%(orientations[ai][2]))
      elif(i==7): fout.write(' %.6f'%(orientations[ai][3]))
      else:       fout.write(' %s'%(t[i]))
    fout.write('\n')
  else:
    fout.write(l)

fout.close()
print('\n# done. conformation coefficients printed to '+out_fname)

