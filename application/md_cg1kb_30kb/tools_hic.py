#!/usr/bin/env python
import numpy as np
from MDAnalysis.analysis.distances import self_distance_array

def kernel_rational(dist,dc=0.5,gamma=1.0,nn=6):
	val=gamma*(dist-dc)
	epsilon=0.000001
	mm=2*nn
	contact=np.where(np.fabs(val-1.)<epsilon,0.5,(1.-val**nn)/(1.-val**mm))
	return contact

# takes as input a list of trajectories
# outputs the average hic map
def traj2hic(trajs,tb=0,te=-1,skip=1,every=1,dc=0.5,gamma=1.0,nn=6,avg=True):
	L=len(trajs[0].atoms)//every
	hic=np.zeros(L*(L-1)//2)
	obs=np.zeros(L*(L-1)//2)
	trajs_hic=[]
	for u in trajs: # loop over the trajectories
		trajs_hic.append([])
		for ts in u.trajectory:
			if(ts.frame%skip!=0): continue
			if(ts.frame<tb): continue
			if(ts.frame>te and te>=0): break
			dist=self_distance_array(u.atoms.positions[::every])
			cmap=kernel_rational(dist,dc,gamma,nn)
			hic+=cmap
			obs+=1.
			if(not avg): trajs_hic[-1].append(np.array(cmap))
		if(not avg): trajs_hic[-1]=np.array(trajs_hic[-1])
	hic/=obs
	if(avg):
		return hic
	else:
		return trajs_hic

# convert array of offdiagonal elements into matrix
def hic2mat(h,otherhalf=False,diag=False):
	L=int((1+np.sqrt(1+8*len(h)))/2)
	if(L*(L-1)//2!=len(h)):
		print('error in hic2mat: L*(L-1)//2!=len(h)')
		return None
	mat=np.zeros((L,L))
	k=-1
	for i in range(L):
		if(diag): mat[i,i]=1.
		for j in range(i+1,L):
			k+=1
			mat[i,j]=h[k]
			if(otherhalf): mat[j,i]=h[k]
	return mat

