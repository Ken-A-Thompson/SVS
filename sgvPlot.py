import numpy as np
import matplotlib.pyplot as plt
import pickle

K = 1000 #max number of parents (positive integer)
n = 2 #number of traits (positive integer)
B = 2 #number of offspring per generation per parent (positive integer)
u = 0.001 #mutation probability per genome (0<u<1)
alpha = 0.1 #mu

sim_id = 'K%d_n%d_B%d_u%r_alpha%r' %(K,n,B,u,alpha)


# load pop data
f = open('pop_%s.pkl' %(sim_id), 'rb')
pop = []
while 1:
    try:
        pop.append(pickle.load(f))
    except EOFError:
        break

# load mut data
g = open('mut_%s.pkl' %(sim_id), 'rb')
mut = []
while 1:
    try:
        mut.append(pickle.load(g))
    except EOFError:
        break

# load gen data
f = open('gen_%s.pkl' %(sim_id), 'rb')
gen = []
while 1:
    try:
        gen.append(pickle.load(f))
    except EOFError:
        break

# make phenotype data
phenos = []
i = 0
while i < len(gen): 
	phenos.append(np.dot(pop[i],mut[i]))
	i = i + 1


# plot phenotypes
i = 0
while i < len(phenos): 
	plt.scatter(phenos[i][:,0],phenos[i][:,1])
	plt.show()
	i = i + 1
plt.show()



