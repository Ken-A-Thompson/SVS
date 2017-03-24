import numpy as np
import matplotlib.pyplot as plt

K = 1000 #max number of parents (positive integer)
n = 2 #number of traits (positive integer)
B = 2 #number of offspring per generation per parent (positive integer)
u = 0.001 #mutation probability per genome (0<u<1)
alpha = 0.1 #mu

sim_id = 'K%d_n%d_B%d_u%r_alpha%r' %(K,n,B,u,alpha)

#get data from simulation
time,pop = np.loadtxt('pop_%s.dat' %(sim_id), skiprows=1, unpack=True)
time,mut = np.loadtxt('mut_%s.dat' %(sim_id), skiprows=1, unpack=True)

with open('pop_%s.dat' %(sim_id)) as f: #trait values
	pop = np.array(f.read().split())
