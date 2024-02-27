#!/usr/bin/env python
import sys
import numpy as np
from MDAnalysis.analysis.distances import self_distance_array

def read_dump(topol_fname,traj_fname,tb=-1,te=-1,return_indices=False):
	# read topology
	positions=[]
	indices={}
	item=''
	N=0
	i_box=0
	box=np.array([0.,0.,0.])
	for l in open(topol_fname,'r').readlines():
		t=l.split()
		if(t[0]=='ITEM:'):
			item=t[1]
			if(item=='BOX'): i_box=0
		elif(item=='BOX'):
			box[i_box]=float(t[1])-float(t[0])
			i_box+=1
		elif(item=='ATOMS'):
			t=l.split()
			ai=int(t[0])
			atype=int(t[1])
			if(atype==1):
				x,y,z=float(t[2]),float(t[3]),float(t[4])
				positions.append([x,y,z])
				indices[ai]=N
				N+=1
	positions=np.array(positions)
	traj=[]
	# read trajectory
	timestep=0
	i=0 # nucleosome counter
	i_box=0
	box=np.array([0.,0.,0.])
	for l in open(traj_fname,'r').readlines():
		t=l.split()
		if(t[0]=='ITEM:'):
			item=t[1]
			if(item=='BOX'): i_box=0
		elif(item=='TIMESTEP' and len(t)==1):
			timestep=int(l)
			if(timestep<tb): continue
			if(timestep>te and te>0): break
			i=0
			positions=np.zeros((N,3))
		elif(item=='BOX'):
			box[i_box]=float(t[1])-float(t[0])
			i_box+=1
		elif(item=='ATOMS' and len(t)==9):
			ai=int(t[0])
			atype=int(t[1])
			if(atype==1):
				inuc=indices[ai]
				x,y,z=float(t[2]),float(t[3]),float(t[4])
				if(inuc<0 or inuc>N):
					print('error: invalid inuc')
					print(t)
					return None
				positions[inuc,0]=x
				positions[inuc,1]=y
				positions[inuc,2]=z
				i+=1 #
				if(i==N): # reached the end of the frame
					traj.append(np.array(positions))
	if(return_indices):
		return np.array(traj),list(indices.keys())
	else:
		return np.array(traj)

# kernel for making the Hi-C map
# 50 bp x 3.4 A = 170 A is the length of the linker DNA in our system
def kernel_rational(dist,sigma=113.0,dc=65.,nn=6):
	gamma=1./sigma
	val=gamma*(dist-dc)
	epsilon=0.000001
	mm=2*nn
	contact=np.where(np.fabs(val-1.)<epsilon,0.5,(1.-val**nn)/(1.-val**mm))
	return contact

# takes as input a list of trajectories
# outputs the average hic map
def traj2hic(trajs,tb=0,te=-1,skip=1,every=1,sigma=113.0,dc=None,nn=6):
	if(dc==None): dc=0.5*sigma
	L=trajs[0].shape[1]//every
	hic=np.zeros(L*(L-1)//2)
	obs=np.zeros(L*(L-1)//2)
	for traj in trajs: # loop over the trajectories
		frame=-1
		for positions in traj:
			frame+=1
			if(frame%skip!=0): continue
			if(frame<tb): continue
			if(frame>te and te>=0): break
			dist=self_distance_array(positions[::every])
			cmap=kernel_rational(dist,sigma,dc,nn)
			hic+=cmap
			obs+=1.
	hic/=obs
	return hic

# takes a trajectory, returns a cmap trajectory
def traj2hictraj(traj,tb=0,te=-1,skip=1,every=1,sigma=113.0,dc=None,nn=6):
	if(dc==None): dc=0.5*sigma
	L=traj.shape[1]//every
	hic=[]
	frame=-1
	for positions in traj:
		frame+=1
		if(frame%skip!=0): continue
		if(frame<tb): continue
		if(frame>te and te>=0): break
		dist=self_distance_array(positions[::every])
		cmap=kernel_rational(dist,sigma,dc,nn)
		hic.append(cmap)
	return np.array(hic)

# convert array of offdiagonal elements into matrix
def hic2mat(h,otherhalf=False,diag=False):
	L=int((1+np.sqrt(1+8*len(h)))/2)
	if(L*(L-1)//2!=len(h)):
		print('error in hic2mat: L*(L-1)//2!=len(h)')
		return None
	mat=np.zeros((L,L))
	k=-1
	for i in range(L):
		if(diag): mat[i,i]=diag
		for j in range(i+1,L):
			k+=1
			mat[i,j]=h[k]
			if(otherhalf): mat[j,i]=h[k]
	return mat

