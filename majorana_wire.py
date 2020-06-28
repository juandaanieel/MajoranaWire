from numpy import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from random import *

#Define parameters of the chain
N = 200	#length
mu = -2	#chemical potential
a = 0.3	#spin orbit
t = 1	#hopping
d = 0.3	#sc gap
z = 0.4	#Zeeman field

#Define paratemers of the quantum dot
mu0 = -0	#dot chemical potential	
mu_n = -0	#normal part chemical potential
dis = 0		#sc length
centerQD = 0	#center of the quantum dot
lengthQD = 0	#lengt/2 
dQD = 0		#distance to the normal part

#SC disable
SCdis = range(dis)


#Plots

def createSpectra(eigenvalues):
	"""Data for the plot of the spectra"""
	plot=[[],[]]
	for h in range(len(eigenvalues)):
		plot[0].append(h)
		plot[1].append(eigenvalues[h])
	return plot

def plotSpectra(array1):
	"""Plot spectra"""
	#plt.title('Dot: '+str(centerQD - leng//2)+' to '+str(centerQD + leng//2 + distance)+'\nQD Well in:'+str(mu0)+' out:'+str(mu_n)+'\nVz: '+str(z))
	#plt.subplot(2,1,1)
	plt.grid(True,which='both')
	#plt.ylabel('Eigenvalue')
	plt.plot(array1[0],array1[1],'ro')
	#plt.subplots_adjust(left=0.15)
	#axes = plt.gca()
	#axes.set_xlim([490,510])
	#axes.set_ylim([-0.1,0.1])
	#plt.subplot(2,1,2)
	#plt.grid(True,which='both')
	plt.ylabel('LDOS')
	plt.xlabel('site')
	#plt.plot(array2[0],array2[1])
	plt.show()

def plot3D(array):
	"""Plot 3D information"""
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1, projection='3d')
	array[0], array[1] = meshgrid(array[0],array[1])
	print shape(array[0]),shape(array[1]),shape(array[2])
	surf = ax.plot_wireframe(array[0], array[1], array[2],cmap=cm.coolwarm,linewidth=0.5, antialiased=False)
	#fig.colorbar(surf, shrink=0.5, aspect=5)
	ax.set_xlabel('Site')
	ax.set_ylabel('Energy')
	ax.set_zlabel('LDOS')
	plt.show()

def surfLDOS(eigenvalues,eigenvectors,N,q):
	"""Obtain the local density of states for a chain of N sites
	in an energy range of \pm q"""
	surf_ldos = [[],[],[]]
	surf_ldos[0] = arange(N)
	for i in range(q,q+N/2):
		print eigenvalues[i]
		surf_ldos[1].append(eigenvalues[i])
		surf_ldos[2].append(inner_product(i,eigenvectors))
	surf_ldos[2]=array(surf_ldos[2])
	return surf_ldos	


def Hamiltonian(N):
	"""Create an empty Hamiltonian"""
	return zeros((4*N,4*N))

def Set_Hamiltonian(ham,N):
	"""
	Fill the Hamiltonian with the following interactions:
		a) on-site chemical potential
		b) Rashba (sigma y) spin orbit interaction
		c) Zeeman field (sigma x)
		d) s-wave superconductivity (tau x)
	"""
	for j in range(N):
		#Diagonal terms
		for i in range(4):
			#rand = uniform(0.010,0.020)
			if i<2:
				ham[i+4*j][i+4*j]=-mu #+ rand
			else:
				ham[i+4*j][i+4*j]=mu #+ rand
		#Hopping terms
		for k in range(4):
			if j<N-1:
				if k<2:			
					ham[k+4*j][k+4*j+4]=-t
					ham[k+4*j+4][k+4*j]=-t
				else:
		                        ham[k+4*j][k+4*j+4]=t
        		                ham[k+4*j+4][k+4*j]=t
		#Spin Orbit
		for f in range(4):
			if j<N-1:
				if f==0:
					ham[f+4*j][f+4*j+4+1]=-a	
					ham[f+4*j+4][f+4*j+1]=a
	
				elif f==1:
        	                        ham[f+4*j][f+4*j+4-1]=a
        	                        ham[f+4*j+4][f+4*j-1]=-a
	
        	                elif f==2:
        	                        ham[f+4*j][f+4*j+4+1]=a
        	                        ham[f+4*j+4][f+4*j+1]=-a
	
        	                elif f==3:
        	                        ham[f+4*j][f+4*j+4-1]=-a
        	                        ham[f+4*j+4][f+4*j-1]=a
	
		#Superconductin Term
		for s in range(4):
			if j not in SCdis:
				if s<2:
					ham[s+4*j][s+4*j+2]=d
        		        else:
        		                ham[s+4*j][s+4*j-2]=d
	
		#Zeeman
		for c in range(4):
			if c==0 or c==2:
				ham[c+4*j][c+4*j+1]=z
			elif c==1 or c==3:
				ham[c+4*j][c+4*j-1]=z
	return ham
				
	

def Print_Hamiltonian(hamiltonian):
	"""Print Hamiltonian Matrix"""
	for row in hamiltonian:
		for val in row:		
			print '{:4}'.format(val),
		print

def inner_product(pos,eigenvectors):
	"""Inner Product in Nambu Space"""
	array=[]
	for i in range(N):
		site1 = 0
		site2 = 0
		for j in range(4):
			site1 += eigenvectors[:,pos][j+4*i]**2
			site2 += eigenvectors[:,pos-1][j+4*i]**2
		#print site1,site2,i
		#array.append(i)
		array.append(site1+site2)
	return array

def eigenv(ham):
	"""Find Eigenvalues and Eigenvectors"""
	eigen = linalg.eig(ham)
	eigenvalues=eigen[0]
	eigenvectors=eigen[1]

	idx = eigenvalues.argsort()[::-1]
	eigenvalues = eigenvalues[idx]
	eigenvectors = eigenvectors[:,idx]
	
	return eigenvalues, eigenvectors


def createDot(hamiltonian,mu_dot,mu_barrier,center,length,distance):
	"""
	Create a quantum dot centered at a point and with a given lengt, and separated of the topological region with a barrier of chemical potential mu_barrier.
	"""
	for j in range(dis):
		if j<center+length+distance and j>center-length:
			for i in range(4):
				if i<2:
					hamiltonian[i+4*j][i+4*j]=-mu_dot
					hamiltonian[i+4*j][i+4*j+2]=0
				else:
					hamiltonian[i+4*j][i+4*j]=mu_dot
					hamiltonian[i+4*j][i+4*j-2]=0
		else:
			for i in range(4):
				if i<2:
					hamiltonian[i+4*j][i+4*j]=-mu_barrier
					hamiltonian[i+4*j][i+4*j+2]=0
				else:
					hamiltonian[i+4*j][i+4*j]=mu_barrier
					hamiltonian[i+4*j][i+4*j-2]=0
	return hamiltonian
