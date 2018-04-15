############################################################################
# Author: 
#
#	This script implements the 3D Yee algorithm to simulate the 
#		Maxwell equations, using abrosbing boundary conditions.
#	This class has a method to calculate the time evolution of
#		a field: time_step. It has an option to change both permitivities
# 		in a spot or a region, and it can implement hard and soft sources.
# 	For future work it is needed to implement support for lossy mediums 
# 		and to use proper boundary conditions
#
############################################################################


#necessary imports
import numpy as np 
import scipy as sp 
import matplotlib.pyplot as pl 
from scipy import special
import time
import os


class yee:
	def __init__(self,L_x,L_y,L_z,N_x,N_y,N_z,n_border=80,c=1.,eps_0=1.,mew_0=1,prints=True):

		# makes sure that there is a even number of points in every grid direction

		self.dx=2.*L_x/N_x
		self.dy=2.*L_y/N_y
		self.dz=2.*L_z/N_z

		#calculates the time interval based on the stability condition

		self.dt= 1/(c*np.sqrt(1/self.dx**2+1/self.dy**2+1/self.dz**2))
		
		self.E=[]
		self.H=[]

		#creates the points of a square grid centered on zero
		X,Y,Z=np.arange(-L_x-n_border*self.dx ,L_x+self.dx+n_border*self.dx,self.dx),np.arange(-L_y-n_border*self.dy,L_y+self.dy+n_border*self.dy,self.dy),np.arange(-L_z-n_border*self.dz,L_z+self.dz+n_border*self.dz,self.dz)

		#indexing ij para ter outputs da forma (X,Y,Z), o default iria trocar a ordem do output
		
		X,Y,Z=np.meshgrid(X,Y,Z,indexing='ij')
		self.grid=np.array([X,Y,Z])
	

		#creates the grid for the electric and magnetic field
		self.E=np.complex128(0.*self.grid)
		self.H=np.complex128(0.*self.grid)

	
	
		#default value for eps and mew is 1. 
		self.eps=np.complex128(self.grid[0]*0.+eps_0)
		self.mew=np.complex128(self.E[0]*0+mew_0)
	

		self.old_source_E=0
		self.old_source_H=0

		
		val=0.5*10**(-10)*n_border**6-3.5*10**(-9)*n_border**5  -5.2001*10**(-4)*n_border**4+2.3199999*10**(-2)*n_border**3 + 0.4599999*n_border**2 - 4.2*n_border+6813.9	
		base=1.028

		#creates the boundary ABC conditions ->changes the eps in the boundary region
		for j in range(n_border):

			placeholder=((((base**(val*j-10))-1.1)*1j)**2+0.123)/100
			
			self.eps[n_border-j,:,:]-=placeholder
			self.eps[self.eps.shape[0]-n_border+j,:,:]-=placeholder

			self.eps[:,n_border-j,:]-=placeholder
			self.eps[:,self.eps.shape[1]-n_border+j,:]-=placeholder

			self.eps[:,:,n_border-j]-=placeholder
			self.eps[:,:,self.eps.shape[2]-n_border+j]-=placeholder

		

		#stores the number of points in the border and the index of the last element in the simulation region
		self.box=[[n_border,(len(X)-1)/2-n_border,N_x],[n_border,(len(Y)-1)/2-n_border,N_y],[n_border,(len(Z)-1)/2-n_border,N_z]]
		self.dict={'x':0,'y':1,'z':2}
		if prints==True:
			if N_x%2!=0:
				N_x+=1
				print('The number of points should be even. N_x increaded by one')

			if N_y%2!=0:
				N_y+=1 
				print('The number of points should be even. N_y increaded by one')

			if N_z%2!=0:
				N_z+=1
				print('The number of points should be even. N_z increaded by one')
			

			print('Points in use: %s'%[N_x,N_y,N_z])
			print('dt=%s'%self.dt)
			

		del X,Y,Z,val,base



	def reset(self):
		'''
		Resets the object to the initial conditions
		'''
		self.E*=0.
		self.H*=0.

		self.eps=np.complex128(self.grid[0]*0.)+1
		self.mew=np.complex128(self.grid[0]*0.)+1
		
		self.old_source_E=0
		self.old_source_H=0

		n_border=self.box[0][0]
		val= -5.20*10**(-4)*n_border**4 + 2.32*10**(-2)*n_border**3 + 0.46*n_border**2 - 4.2*n_border+6813.9
		base=1.038
		
		#creates the boundary ABC conditions ->changes the eps in the boundary region
		for j in range(n_border):

			placeholder=((((base**(val*j-10.))*1j)**2))
			self.eps[n_border-j,:,:]-=placeholder
			self.eps[self.eps.shape[0]-n_border+j,:,:]-=placeholder

			self.eps[:,n_border-j,:]-=placeholder
			self.eps[:,self.eps.shape[1]-n_border+j,:]-=placeholder

			self.eps[:,:,n_border-j]-=placeholder
			self.eps[:,:,self.eps.shape[2]-n_border+j]-=placeholder

		del val,base,n_border

	def set_params(self,param,region,dist=[0,0,0],sides=[0,0,0]): 
		'''
		This function allows the user to change the values of epsilon and mew. 
		The input param specifies the permitivity to change: eps for epsilon and mew for mew
	
		dist: distance from the origin to the center of the back side of the box
		sides: lenght of the sides.

		region: the values for the parameter to change. It can be a number and then it changes all the values in the defined region to that or
		it can be an array, allowing the user to individually define the permitivity. However, it's shape should match the one from the 
		desired region
		'''

		#stores the physical boundaries of the simulation area
		p0=int(self.box[0][2])
		p1=int(self.box[1][2])
		p2=int(self.box[2][2])
		n_border=self.box[0][0]
		
			
		if (abs(dist[0])<self.dx and dist[0]!=0) or (abs(dist[1])<self.dy and dist[1]!=0) or(abs(dist[2])!=0 and dist[2]<self.dz):
			raise ValueError('Vector components smaller than the spacial step')

		# finds the center of the back side of the box in which the field will be changed
		x=(self.eps.shape[0]-1)/2.+dist[0]/self.dx
		y=(self.eps.shape[1]-1)/2.+dist[1]/self.dy
		z=(self.eps.shape[2]-1)/2.+dist[2]/self.dz
		
		
		#finds the y and z corner of the box
		y-=0.5*sides[1]/self.dy
		z-=0.5*sides[2]/self.dz
		
		x=int(x)
		y=int(y)
		z=int(z)

		#finds the opposite corner of the box
		x1=int(x+sides[0]/self.dx)
		y1=int(y+sides[1]/self.dy)
		z1=int(z+sides[2]/self.dz)

		#makes sure that the region does not overlap with the boundary region
		if x>self.box[0][2]+n_border or y>self.box[1][2]+n_border or z>self.box[2][2]+n_border:
			raise ValueError('Box corner outside the simulation zone')
		if x <n_border or y<n_border or z<n_border:
			raise ValueError('Box corner outside the simulation zone')
		if x1>self.box[0][2]+n_border or y1>self.box[1][2]+n_border or z1>self.box[2][2]+n_border:
			raise ValueError('Box outside the simulation zone')
		if x1 <n_border or y1<n_border or z1<n_border:
			raise ValueError('Box corner outside the simulation zone')


		if param=='eps':
			self.eps[x:x1+1,y:y1+1,z:z1+1]=region

		elif param=='mew':
			self.mew[x:x1+1,y:y1+1,z:z1+1]=region		
		else:
				raise TypeError('Such parameter does not exist')

		del n_border,p0,p1,p2,x,y,z,x1,y1,z1

	
	def force_field(self,field,axis,region):

		'''
		Sets the eletrical or magnetical field to the input array (region). 
		field specifies the field to be changed: E(eletric) or H (magnetic)

		val should be an array with the same shape as the grid of the axis in question. 
		This function only allows to change the fields on the simulation zone: can't change on the boundary area.
		It is not possible to replace/force specific areas of the simulation region
		
		'''
		p0=int(self.box[0][2])
		p1=int(self.box[1][2])
		p2=int(self.box[2][2])

		n_border=self.box[0][0]
		i=0
		for reg in region:
			
			if reg.shape!=self.E[self.dict[axis],0,n_border:n_border+int(p1)+1,n_border:n_border+1+int(p2)].shape:
				raise TypeError('Make sure that the input regions have the same shape as the simulation area')
		
		
		if field=='E':
			for j in range(n_border,n_border+p0+1):
				self.E[self.dict[axis],j,n_border:n_border+int(p1)+1,n_border:n_border+1+int(p2)]=region[i]
				i+=1

			
		elif field=='H':
			for j in range(n_border,n_border+p0+1):
				self.H[self.dict[axis],j,n_border:n_border+int(p1)+1,n_border:n_border+1+int(p2)]=region[i]
				i+=1

		else:
			raise TypeError('Input field does not exist')


		del p0,p1,p2,n_border,i


	def dot_source(self,field,axis,value,dist=[0,0,0],warnings=True,verifications=True):

		'''
		Changes the eletrical or magnetic field value at the chosen point to the chosen value.
		dist -> distance to the origin of the -axis- component of -field-.

		This is a hard source, since the wave will "see it as a wall", being reflected by it. It should be turned off after it's in the system, since
		it will cause unwanted reflections.

		If warnings==False then it does not alert the user of approximation issues
		If verifications==False then it doesn't check if the dot is inside the simulation zone
		'''
		if (dist[0]!=0 and dist[0]<self.dx) or (dist[1]<self.dy and dist[1]!=0) or (dist[2]<self.dz and dist[2]!=0):
			if dist!=[0,0,0]:
				raise ValueError('Vector components smaller than the spacial step')

		n_border=self.box[0][0]
		if dist[0]>self.dx*(self.box[0][2]+n_border) or dist[1]>self.dy*(self.box[1][2]+n_border) or dist[2]>self.dz*(self.box[2][2]+n_border):
			raise ValueError('Vector component bigger than the simulation zone')

		
		x=int(self.box[0][0]+self.box[0][2]/2.+dist[0]/self.dx)
		y=int(self.box[1][0]+self.box[1][2]/2.+dist[1]/self.dy)
		z=int(self.box[2][0]+self.box[2][2]/2.+dist[2]/self.dz)

		if verifications==True:

			if x<n_border or y<n_border or z<n_border:
				raise ValueError('Dot outside the simulation range')

			if x>n_border+2*self.box[0][2] or y>n_border+2*self.box[1][2] or z>n_border+2*self.box[2][2]:
				raise ValueError('Dot outside the simulation range')

		if field=='E':
			self.E[self.dict[axis],x,y,z]=value
		elif field=='H':
			self.H[self.dict[axis],x,y,z]=value


		else:
			raise TypeError('No such field')
		if dist[0]%self.dx!=0 or dist[1]%self.dy!=0 or dist[2]%self.dz!=0:
			if warnings==True:
				print('WARNING:since the distance is not a multiple of the spacial steps there will be rounding of numbers.Placement of the dot source may not be on the desired spot')
		del x,y,z


	def soft_source(self,field,axis,region,dist=[0,0,0],sides=[0,0,0],verifications=True,warnings=True):

		'''
		Changes the eletrical or magnetic field value at the chosen point to the chosen value.
		dist -> distance to the origin of the -axis- component of -field-.

		This is a hard source, since the wave will "see it as a wall", being reflected by it. It should be turned off after it's in the system, since
		it will cause unwanted reflections.

		If verications==False then it does not check if the region is inside the simulation zone
		'''

		p0=int(self.box[0][2])
		p1=int(self.box[1][2])
		p2=int(self.box[2][2])
		n_border=self.box[0][0]
		

		# finds the center of the back side of the box in which the field will be changed
		x=(self.E[self.dict[axis]].shape[0]-1)/2.+dist[0]/self.dx
		y=(self.E[self.dict[axis]].shape[1]-1)/2.+dist[1]/self.dy
		z=(self.E[self.dict[axis]].shape[2]-1)/2.+dist[2]/self.dz
		
		#finds the y and z corner of the box
		y-=0.5*sides[1]/self.dy
		z-=0.5*sides[2]/self.dz
		
		x=int(x)
		y=int(y)
		z=int(z)

		#finds the opposite corner of the box
		x1=int(x+sides[0]/self.dx)
		y1=int(y+sides[1]/self.dy)
		z1=int(z+sides[2]/self.dz)
		
		if verifications:
			if (abs(dist[0])<self.dx and dist[0]!=0) or (abs(dist[1])<self.dy and dist[1]!=0) or(abs(dist[2])!=0 and dist[2]<self.dz):
				raise ValueError('Vector components smaller than the spacial step')
			if x>self.box[0][2]+n_border or y>self.box[1][2]+n_border or z>self.box[2][2]+n_border:
				raise ValueError('Box corner outside the simulation zone')
			if x <n_border or y<n_border or z<n_border:
				raise ValueError('Box corner outside the simulation zone')

			if x1>self.box[0][2]+n_border or y1>self.box[1][2]+n_border or z1>self.box[2][2]+n_border:
				raise ValueError('Box outside the simulation zone')
			if x1 <n_border or y1<n_border or z1<n_border:
				raise ValueError('Box corner outside the simulation zone')

		if warnings:
			if dist[0]%self.dx!=0 or dist[1]%self.dy!=0 or dist[2]%self.dz!=0:
				print('WARNING:since the distance is not a multiple of the spacial steps there will be rounding of numbers.Placement of the source may not be on the desired spot')
			if sides[0]%self.dx!=0 or sides[1]%self.dy!=0 or sides[2]%self.dz!=0:
				print('WARNING:since one of the sides is not a multiple of the spacial steps there will be rounding of numbers.Placement of the source may not be on the desired spot')


		if field=='E':
			self.E[self.dict[axis],x:x1+1,y:y1+1,z:z1+1]+=(region-self.old_source_E)
			self.old_source_E=region


		elif field=='H':
			self.H[self.dict[axis],x:x1+1,y:y1+1,z:z1+1]+=(region-self.old_source_H)
			self.old_source_H=region
		else:
			raise TypeError('Input field does not exist')

		del p0,p1,p2,n_border,x,y,z,x1,y1,z1

	#calculates the step in E and H  
	def time_step(self,field,dt=[]):
		'''
		Calculates the time step for the input field. Before doing that, it's verified that the values of permitivities were changed from the default ones.
		'''
		if dt==[]:
			dt=self.dt

		if field=='E':	
			coef=np.complex128(dt)/self.eps
			self.E[0]+= coef*(((self.H[2]-np.roll(self.H[2],1,axis=1))/self.dy)-((self.H[1]-np.roll(self.H[1],1,axis=2))/self.dz))
			self.E[1]+= coef*(((self.H[0]-np.roll(self.H[0],1,axis=2))/self.dz)-((self.H[2]-np.roll(self.H[2],1,axis=0))/self.dx))
			self.E[2]+= coef*(((self.H[1]-np.roll(self.H[1],1,axis=0))/self.dx)-((self.H[0]-np.roll(self.H[0],1,axis=1))/self.dy))

		elif field=='H':
			coef=np.complex128(dt)/self.mew
			self.H[0]+= coef*(((np.roll(self.E[1],-1,axis=2)-self.E[1])/self.dz)-((np.roll(self.E[2],-1,axis=1)-self.E[2])/self.dy))
			self.H[1]+= coef*(((np.roll(self.E[2],-1,axis=0)-self.E[2])/self.dx)-((np.roll(self.E[0],-1,axis=2)-self.E[0])/self.dz))
			self.H[2]+= coef*(((np.roll(self.E[0],-1,axis=1)-self.E[0])/self.dy)-((np.roll(self.E[1],-1,axis=0)-self.E[1])/self.dx))
		else:
			raise TypeError('No such field')

	def result(self,entry):
		'''
		return the vale of the attributes. If entry== name of a method then it returns the method
		'''

		return getattr(self,entry)




	

	



