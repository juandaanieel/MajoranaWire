from numpy import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from random import *

######################################################
#QUANTUM#DOT#SPECTROSCOPY#TOPOLOGICAL#RASHBA#NANOWIRE#
######################################################
######Developed by: Juan Daniel Torres################
######################################################

"""1) Initialize the parameters of the Hamiltonian class
        N: number of sites
        mu: chemical potential
        delta: SC gap
        alpha: SO magnitude
        Vz: magnetic field
        t: hopping"""

class Hamiltonian(object):
    """Hamiltonian for the nanowire"""
    def __init__(self,N,mu,delta,alpha,Vz,t):
        
        self.H = []

        self.N = N 
        self.mu = mu
        self.delta = delta
        self.alpha = alpha 
        self.Vz = Vz
        self.t = t
        
        self.spectra = []
        self.wavefuctions = []
        self.site_dos = [[],[]]
        self.density_of_states =[]

    def initialize_H(self):
        """Create H as 4Nx4N zero matrix"""
        self.H = zeros((4*self.N,4*self.N))

    """2) Build the interactions in the Nanowire """
    def set_self_energy(self,N):
        for i in range(4):
            if i<2:
                self.H[N*4+i][N*4+i] = -self.mu
            else:
                self.H[N*4+i][N*4+i] = self.mu

    def set_hopping(self,N):
        for i in range(4):
            if N<self.N-1:
                if i<2:
                    self.H[4*N+i][4*N+i+4] = -self.t
                    self.H[4*N+i+4][4*N+i] = -self.t
                else:
                    self.H[4*N+i][4*N+i+4] = self.t
                    self.H[4*N+i+4][4*N+i] = self.t
    
    def set_zeeman(self,N):
        for i in range(4):
            if i == 0 or i == 2:
                self.H[4*N+i][4*N+i+1] = self.Vz
            else:
                self.H[4*N+i][4*N+i-1] = self.Vz

    def set_SC(self,N,control):
        for i in range(4):
            if i<2:
                self.H[4*N+i][4*N+i+2] = self.delta*control 
            else:
                self.H[4*N+i][4*N+i-2] = self.delta*control
    
    def set_SOI(self,N):
        for i in range(4):
            if N<self.N-1:
                if i == 0:
                    self.H[4*N+i][4*N+i+4+1] = -self.alpha
                    self.H[4*N+i+4][4*N+i+1] = self.alpha
                elif i == 1:
                    self.H[4*N+i][4*N+i+4-1] = self.alpha
                    self.H[4*N+i+4][4*N+i-1] = -self.alpha
                elif i == 2:
                    self.H[4*N+i][4*N+i+4+1] = self.alpha
                    self.H[4*N+i+4][4*N+i+1] = -self.alpha
                elif i == 3:
                    self.H[4*N+i][4*N+i+4-1] = -self.alpha
                    self.H[4*N+i+4][4*N+i-1] = self.alpha

    """3) Set the full Hamiltonian"""
    def set_H(self):
        for i in range(self.N):
            self.set_self_energy(i)
            self.set_hopping(i)
            self.set_SC(i,1)
            self.set_zeeman(i)
            self.set_SOI(i)
    """4) A spectroscopy measurement will be performed on the nanowire
          by removing partially the superconductivity, creating a normal part NP          and building locally a quantum dot."""
    def create_normal_part(self,length,diff):
        """Disable Superconductivity until site=length
           Increase energy self energy/barrier"""
        self.mu = self.mu - diff
        for i in range(length):
            self.set_self_energy(i)
            self.set_SC(i,0)
        self.mu = self.mu + diff

    def create_dot(self,length,start,depth):
        """External potential applied from start to start+length"""
        self.mu = self.mu - depth
        for i in range(start,start+length):
            self.set_self_energy(i)
        self.mu = self.mu +depth

    """5) Set of functions to analyze the spectra of the Hamiltonian"""
    def get_dos(self):
        """DOS at each site at a given energy"""
        for i in range(self.N):
            self.density_of_states.append(self.inner_product(4*i))

    def print_H(self):
        for row in self.H:
            for val in row:
                print '{:4}'.format(val),
            print

    def diagonalize(self):
        """Get eigenvalues and wavefuctions of H"""
        eigen = linalg.eig(self.H)
        self.spectra = eigen[0]
        self.wavefuctions = eigen[1]
        
        idx = self.spectra.argsort()[::-1]
        self.spectra = self.spectra[idx]
        self.wavefuctions = self.wavefuctions[:,idx]

    def inner_product(self,pos):
    	array=[]
    	for i in range(self.N):
    		site1 = 0
    		for j in range(4):
    			site1 += self.wavefuctions[:,pos][j+4*i]**2
    		array.append(site1)
	self.site_dos[1] = array
	self.site_dos[0] = range(self.N)

    def frange(self,start,stop,step):
        array = []
        i = start
        while i<stop:
            array.append(i)
            i += step
        return array
        
    def plot_spectra(self):
        plt.title('spectra')
        plt.grid(True,which='both')
        plt.plot(range(4*self.N),self.spectra,'ro')
        plt.xlabel('site*')
        plt.ylabel('eigenvalue')
        plt.show()

    def plot_dos(self):
        plt.title('LDOS')
        plt.grid(True)
        plt.plot(range(self.N),self.density_of_states,'ro')
        plt.xlabel('site')
        plt.ylabel('amplitude')
        plt.show()

    def plotWF(self,array1,array2,energy):
        plt.title('Density of States at Energy: '+str(energy))
        plt.grid(True,which='both')
        plt.xlabel('Site')
        plt.ylabel('Amplitude')
        plt.plot(array1,array2)
        plt.show()

    def plotDepts(self,array1,array2):
        plt.title('DOS vs QD Chemical Potential')
        plt.grid(True,which='both')
        plt.xlabel('Chemical Potential')
        plt.ylabel('Energy')
        plt.plot(array1,array2,'b+')
        plt.show()

