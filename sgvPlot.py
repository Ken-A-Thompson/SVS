import numpy as np
import matplotlib.pyplot as plt
import pickle
import matplotlib.cm as cm

K = 1000 #max number of parents (positive integer)
n = 2 #number of traits (positive integer)
B = 2 #number of offspring per generation per parent (positive integer)
u = 0.001 #mutation probability per genome (0<u<1)
alpha = 0.1 #mu

sim_id = 'K%d_n%d_B%d_u%r_alpha%r' %(K,n,B,u,alpha)

data_dir = 'data'

# load pop data
f = open('%s/pop_%s.pkl' %(data_dir,sim_id), 'rb')
pop = []
while 1:
    try:
        pop.append(pickle.load(f))
    except EOFError:
        break

# load mut data
g = open(('%s/mut_%s.pkl' %(data_dir,sim_id), 'rb')
mut = []
while 1:
    try:
        mut.append(pickle.load(g))
    except EOFError:
        break

# load gen data
f = open('%s/gen_%s.pkl' %(data_dir,sim_id), 'rb')
gen = []
while 1:
    try:
        gen.append(pickle.load(f))
    except EOFError:
        break

# make phenotype data
phenos = []
for i in range(len(gen)):
	phenos.append(np.dot(pop[i],mut[i]))


# plot phenotypes
colors = cm.rainbow(np.linspace(0, 1, len(phenos)))
for i, c in zip(range(len(phenos)), colors):
    plt.scatter(phenos[i][:,0],phenos[i][:,1], color=c)

plt.show()



