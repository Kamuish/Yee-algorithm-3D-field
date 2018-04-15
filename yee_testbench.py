#############################################################
# Author: 
#This script works as a testing environment for the class yee.
#We will study the effects of changes of different phisical
#conditions (different permitivities and placement/types of 
#sources)
#
#############################################################

from yee import yee
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.patches as patches
import os
import time



##############################################################################
#path to save the folders
path=''


#wave parameter:
lambd=0.3
"should be between 10 and 20 samples for wavelength"
sampl=lambd/20

#Box construction:
Lx=0.5
Ly=2
Lz=1

#number of points in the boundary region
n_border=80

#number of iterations, each corresponding to a time step of dt
simul_time=201


#Dot source
dot_sources=False

#Soft sinusoidal source
soft_source_1=False

#soft constant source
soft_source_2=False

#soft source on a simulation area with a box of different permitivity
different_parameters=True

#Time analysis of 1 iteration for different parameters
time_analysis=False


# Lenght sizes of the Y component and different border sizes for the time analysis
L_y_time_simulation=[0.5,1,1.5]
border_time_simulation=[0,80,100,120,150]


#############################################################################################

Nx=int(2*Lx/sampl)
Ny=int(2*Ly/sampl)
Nz=int(2*Lz/sampl)



for name in ['E','H']:
	if os.path.isdir(path+name):
		for name_2 in ['x','y','z']:
				if not os.path.isdir(path+name+'/'+name_2):
					    os.makedirs(path+name+'/'+name_2)

	else:
		for name_2 in ['x','y','z']:
			os.makedirs(path+name+'/'+name_2)


def save_images(path,field='',t=0,boundaries=True,param=''):

	''''
	Function to plot and save the data corresping the eletrical or magnetical field.magnetical
	If boundaries==True then the plot uses the entire region, including the boundary conditions. If it's False
	then it uses only the computational region

	'''

	if param=='eps':
		Y,Z=np.linspace(-Ly,Ly,Ny+1),np.linspace(-Lz,Lz,Nz+1)

		pl.figure()
		pl.title('Eletrical permitivity')
		pl.contourf(Z,Y,np.absolute(v.eps[int(n_border+Nx/2)][int(n_border):int(n_border+Ny+1),int(n_border):int(n_border+Nz+1)]),1)
		pl.colorbar()
		pl.savefig(path+'/eps.png')
		pl.close()

	

	elif param=='mew':
		Y,Z=np.linspace(-Ly,Ly,Ny+1),np.linspace(-Lz,Lz,Nz+1)
		pl.figure()
		pl.title('magnetical permitivity')
		pl.contourf(Z,Y,np.absolute(v.mew[int(n_border+Nx/2)][int(n_border):int(n_border+Ny+1),int(n_border):int(n_border+Nz+1)]),1)
		pl.colorbar()
		pl.savefig(path+'mew.png')
		pl.close()
	elif param!='':
		raise ValueError('No such field')

	if boundaries==True and param=='':
		Y,Z=np.linspace(-Ly-n_border*dy,Ly+n_border*dy,Ny+1+2*n_border),np.linspace(-Lz-n_border*dz,Lz+n_border*dz,Nz+1+2*n_border)		

		fig = pl.figure()
		ax2=fig.add_subplot(111)
		pl.title('X component')
		im=np.real(v.result(field)[0][int(n_border+Nx/2)])
		pl.imshow(im)
		cbar=pl.colorbar()
		cbar.set_clim(-0.0005, 0.0005)	
		rect=patches.Rectangle((n_border+Nz,n_border+Ny),-Nz,-Ny,fill=False)
		ax2.add_patch(rect)		

		fig.savefig(path+field+'/x/%s.png'%t)
		pl.close()

		fig = pl.figure()
		ax2=fig.add_subplot(111)
		pl.title('Y component')
		im=np.real(v.result(field)[1][int(n_border+Nx/2)])
		pl.imshow(im)
		cbar=pl.colorbar()
		cbar.set_clim(-0.0005, 0.0005)	
		rect=patches.Rectangle((n_border+Nz,n_border+Ny),-Nz,-Ny,fill=False)
		ax2.add_patch(rect)	
		fig.savefig(path+field+'/y/%s.png'%t)
		pl.close()	

		fig = pl.figure()
		ax2=fig.add_subplot(111)
		pl.title('Z component')
		im=np.real(v.result(field)[2][int(n_border+Nx/2)])
		pl.imshow(im)
		cbar=pl.colorbar()
		cbar.set_clim(-0.0005, 0.0005)	
		rect=patches.Rectangle((n_border+Nz,n_border+Ny),-Nz,-Ny,fill=False)
		ax2.add_patch(rect)		
		fig.savefig(path+field+'/z/%s.png'%t)
		pl.close()

		del Y,Z,rect

	elif boundaries==False and param=='':
		Y,Z=np.linspace(-Ly,Ly,Ny+1),np.linspace(-Lz,Lz,Nz+1)


		contorno=np.absolute(v.eps[int(n_border+Nx/2)][int(n_border):int(n_border+Ny+1),int(n_border):int(n_border+Nz+1)])
		
		pl.figure()
		pl.title('X component')
		pl.contour(Z,Y,contorno,1,linewidths=2,colors=('green'))
		pl.contourf(Z,Y,np.real(v.result(field)[0][int(n_border+Nx/2)][int(n_border):int(n_border+Ny+1),int(n_border):int(n_border+Nz+1)]),200)
		cbar=pl.colorbar()
		cbar.set_clim(-0.0005, 0.0005)		
		pl.savefig(path+field+'/x/no_b%s.png'%t)

		pl.clf()
		pl.title('Y component')
		pl.contour(Z,Y,contorno,1,linewidths=2,colors=('green'))
		pl.contourf(Z,Y,np.real(v.result(field)[1][int(n_border+Nx/2)][int(n_border):int(n_border+Ny+1),int(n_border):int(n_border+Nz+1)]),200)
		cbar=pl.colorbar()
		cbar.set_clim(-0.0005, 0.0005)			
		pl.savefig(path+field+'/y/no_b%s.png'%t)



		pl.clf()
		pl.title('Z component')
		pl.contour(Z,Y,contorno,1,linewidths=2,colors=('green'))
		pl.contourf(Z,Y,np.real(v.result(field)[2][int(n_border+Nx/2)][int(n_border):int(n_border+Ny+1),int(n_border):int(n_border+Nz+1)]),200)
		cbar=pl.colorbar()
		cbar.set_clim(-0.0005, 0.0005)		
		pl.savefig(path+field+'/z/no_b%s.png'%t)
		pl.clf()
		pl.close()
		del contorno

	elif param=='' and boundaries not in (True,False):
		raise TypeError('Boundaries should be True or False')


v=yee(Lx,Ly,Lz,Nx,Ny,Nz,n_border=n_border)

dx=v.dx
dy=v.dy
dz=v.dz

dt=v.dt


if dot_sources:

	dist=[0,0,0]
	sides=[0.05,0.05,0.05]
	save_images(path,param='eps')
	save_images(path,param='mew')


	#v.set_params('eps',region=2,dist=dist,sides=sides)
	#X,Y,Z=np.arange(-Lx,Lx+dx,dx),np.arange(-Ly,Ly+dy,dy),np.arange(-Lz,Lz+dz,dz)

	v.dot_source('E','z',1.,dist=[0,0,0])
	v.time_step('H')

	file=open(path+'info.txt','w')
	file.write('box:%s'%[Lx,Ly,Lz])
	file.write('steps:%s'%[v.dx,v.dy,v.dz])
	file.write('lambda:%s'%lambd)
	file.close()

	for t in range(simul_time):
		v.time_step('E')
		v.time_step('H')
		
		if t%5==0 or t==0:
			
			save_images(path,'E',t)
			save_images(path,'H',t)
		print(t)
		
		
if soft_source_1:

	dist=[0,0,0]
	sides=[0.05,0.10,0.25]
	
	X,Y,Z=np.arange(-Lx,Lx+dx,dx),np.arange(-Ly,Ly+dy,dy),np.arange(-Lz,Lz+dz,dz)

	#v.set_params('eps',region=2,dist=dist,sides=sides)


	save_images(path,param='eps')
	save_images(path,param='mew')

	
	val=np.complex128(0.5*np.sin((((2.*np.pi)/lambd))*0))	
	v.soft_source('E','z',val,dist=[0,0,0],verifications=True)
	v.time_step('H',dt/2)
	time=0

	file=open(path+'info.txt','w')
	file.write('box:%s'%[Lx,Ly,Lz])
	file.write('lambda:%s'%lambd)
	file.close()

	for t in range(1,simul_time):
		
		val=np.complex128(0.5*np.sin((((2.*np.pi)/lambd))*time))	
		v.soft_source('E','z',val,dist=[0,0,0],verifications=False,warnings=False)

		v.time_step('E')
		v.time_step('H')
		
		if t%5==0 or t==1:
	
			save_images(path,'E',t,boundaries=True)
			save_images(path,'H',t,boundaries=True)

		print(t)
		time+=v.dt

	file=open(path+'info.txt','w')
	file.write('Simulation time:%s seconds,dt=%s'%(time,v.dt))
	file.close()
		
	
if different_parameters:

	dist=[0,0,0]
	sides=[0.05,1,0.3]
	
	X,Y,Z=np.arange(-Lx,Lx+dx,dx),np.arange(-Ly,Ly+dy,dy),np.arange(-Lz,Lz+dz,dz)

	v.set_params('eps',region=2,dist=dist,sides=sides)


	save_images(path,param='eps')
	save_images(path,param='mew')

	
	val=np.complex128(0.5*np.sin((((2.*np.pi)/lambd))*0))	
	v.soft_source('E','z',val,dist=[0,0,0],verifications=True)
	v.time_step('H',dt/2)
	time=0

	file=open(path+'info.txt','w')
	file.write('box:%s'%[Lx,Ly,Lz])
	file.write('lambda:%s'%lambd)
	file.close()

	for t in range(1,simul_time):
		
		val=np.complex128(0.5*np.sin((((2.*np.pi)/lambd))*time))	
		v.soft_source('E','z',val,dist=[0,0,0],verifications=False,warnings=False)

		v.time_step('E')
		v.time_step('H')
		
		if t%5==0 or t==1:
	
			save_images(path,'E',t,boundaries=True)
			save_images(path,'H',t,boundaries=True)

		print(t)
		time+=v.dt

	file=open(path+'info.txt','w')
	file.write('Simulation time:%s seconds,dt=%s'%(time,v.dt))
	file.close()
		
	
if soft_source_2:

	dist=[0,0,0]
	sides=[0.05,0.10,0.25]
	
	X,Y,Z=np.arange(-Lx,Lx+dx,dx),np.arange(-Ly,Ly+dy,dy),np.arange(-Lz,Lz+dz,dz)

	#v.set_params('eps',region=2,dist=dist,sides=sides)


	save_images(path,param='eps')
	save_images(path,param='mew')

	val=1
	v.soft_source('E','z',val,dist=[0,0,0],verifications=False)
	v.time_step('H',dt/2)
	time=0

	for t in range(1,simul_time):
		
		val=1	
		v.soft_source('E','z',val,dist=[0,0,0],verifications=False)

		v.time_step('E')
		v.time_step('H')
		
		if t%5==0 or t==1:
	
			save_images(path,'E',t,boundaries=True)
			save_images(path,'H',t,boundaries=True)

		print(t)
		time+=v.dt

if time_analysis:

	times=[]
	times_2=[]

	for j in range(len(L_y)):
		times.append([])
		
	
	for Ly in L_y_time_simulation:
		
		for n_border in border_time_simulation:
			print(Ly,n_border)
			Lx=0.5
			
			Lz=1

			Nx=int(2*Lx/sampl)
			Ny=int(2*Ly/sampl)
			Nz=int(2*Lz/sampl)

			v=yee(Lx,Ly,Lz,Nx,Ny,Nz,n_border=n_border,prints=False)


			dist=[0,0,0]
			sides=[0.05,0.05,0.05]
			
			v.dot_source('E','z',1.,dist=[0,0,0])
			v.time_step('H')


			t1=time.time()
			for t in range(1):
				v.time_step('E')
				v.time_step('H')	
			tf=time.time()

			
			times[L_y.index(Ly)].append(tf)
		times_2.append(Ny)
		

	fig=pl.figure()
	print(times_2)
	for j in range(len(times)): 
		pl.plot(border,times[j],label='Ly= %s points'% str(times_2[j]))

	pl.xlabel('Border points')
	pl.ylabel('Time(s)')
	pl.legend()

	pl.show()


	


