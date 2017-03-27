import numpy as np
import matplotlib.pyplot as plt
import pickle
import matplotlib.cm as cm

K = 1000 #max number of parents (positive integer)
n = 2 #number of traits (positive integer)
B = 2 #number of offspring per generation per parent (positive integer)
u = 0.001 #mutation probability per genome (0<u<1)
alpha = 0.01 #mutation SD

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
g = open('%s/mut_%s.pkl' %(data_dir,sim_id), 'rb')
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

# mean phenotypes
mean_phenos = []
for i in range(len(phenos)):
    # plt.plot(np.mean(phenos[i][:,0]),np.mean(phenos[i][:,1]), color='black')
    # plt.plot(np.mean(phenos[i][:,0]),np.mean(phenos[i][:,1]), color='black')
    mean_phenos.append(np.mean(phenos[i], axis=0))

mean_pheno = np.array(mean_phenos)

# save pheno data as CSV for R (can take a little time!)
# import csv
# with open("%s/phenos_%s.csv" %(data_dir,sim_id), "w") as f:
#     writer = csv.writer(f)
#     writer.writerows(phenos)

# plot all phenotypes
colors = cm.rainbow(np.linspace(0, 1, len(phenos)))
for i, c in zip(range(len(phenos)), colors):
    plt.scatter(phenos[i][:,0],phenos[i][:,1], color=c)

# plot mean phenotypes
plt.plot(mean_pheno[:,0],mean_pheno[:,1], color='black')

# save plot
plt.savefig('%s/plot_%s.png' %(data_dir,sim_id))



