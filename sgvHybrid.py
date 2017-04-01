#Author: Matthew Osmond <mmosmond@zoology.ubc.ca>
#Description: Adaptation from standing genetic variance (SGV) in Fisher's geometric model, implications for hybrids
#Make hybrids from populations that have adapted

import numpy as np
import matplotlib.pyplot as plt
import pickle
import matplotlib.cm as cm

######################################################################
##LOAD DATA FROM ADAPTED POPNS##
######################################################################

K = 1000 #max number of parents (positive integer)
n = 2 #number of traits (positive integer)
B = 2 #number of offspring per generation per parent (positive integer)
u = 0.001 #mutation probability per genome (0<u<1)
alpha = 0.01 #mutation SD
maxgen = 1000

pop = []
mut = []
opt1s = [[0.5] * n, [-0.5] * n]
for i in range(len(opt1s)):
    # optimum for simulation i
    opt1 = opt1s[i]
    # filename and directory of data
    sim_id = 'K%d_n%d_B%d_u%r_alpha%r_gens%r_opt%s_adapt' %(K,n,B,u,alpha,maxgen,'-'.join(str(e) for e in opt1))
    data_dir = 'data'
    # load pop data
    f = open('%s/pop_%s.pkl' %(data_dir,sim_id), 'rb')
    while 1:
        try:
            pop.append(pickle.load(f))
        except EOFError:
            break
    # load mut data
    g = open('%s/mut_%s.pkl' %(data_dir,sim_id), 'rb')
    while 1:
        try:
            mut.append(pickle.load(g))
        except EOFError:
            break

######################################################################
##PARENTAL PHENOTYPES##
######################################################################

#Make parental phenotype data
phenos = [  ]
for i in range(len(opt1s)):
	phenos.append(np.dot(pop[i],mut[i]))

# mean parental phenotypes
mean_phenos = []
for i in range(len(opt1s)):
    mean_phenos.append(np.mean(phenos[i], axis=0))

mean_pheno = np.array(mean_phenos) #reformat to numpy array

# save pheno data as CSV for R (can take a little time!)
import csv
sim_id = 'K%d_n%d_B%d_u%r_alpha%r_gens%r' %(K,n,B,u,alpha,maxgen)
with open("%s/parent_phenos_%s.csv" %(data_dir,sim_id), "w") as f:
    writer = csv.writer(f)
    writer.writerows(phenos)

# ######################################################################
# ##MAKE HYBRIDS##
# ######################################################################

nHybrids = 100 #number of hybrids to make

i=0
j=1
offphenos = []
for k in range(nHybrids):
    # random parents
    randpari = pop[i][np.random.choice(len(pop[i]))] 
    randparj = pop[j][np.random.choice(len(pop[j]))]
    # random parent phenotypes
    phenpari = np.dot(randpari,mut[i]) 
    phenparj = np.dot(randparj,mut[j])
    # mutations held by random parents
    mutpari = mut[i] * randpari[:,None]
    mutparj = mut[j] * randparj[:,None]
    setA = set(tuple(x) for x in mutpari)
    setB = set(tuple(x) for x in mutparj)
    # mutations shared by two parents (all in offspring)
    sharedmuts = np.array([x for x in setA & setB])
    # mutations not shared by two parents
    unsharedmuts = np.array([x for x in setA ^ setB])
    # which unshared mutations in offspring (free recombination between all loci, therefore gets each with 0.5 probability)
    randmuts = np.random.randint(2, size=(len(unsharedmuts)))
    offmuts = unsharedmuts * randmuts[:,None]
    if len(offmuts) < 1:
        offmuts = np.array([[0]*n]) #give something in case empty
    # offspring phenotype
    offphenos.append(sum(np.append(sharedmuts,offmuts,axis=0)))

offphenos01 = np.array(offphenos) #reformat correctly
mean_offpheno01 = np.mean(offphenos,axis=0) #mean

sim_id = 'K%d_n%d_B%d_u%r_alpha%r_gens%r' %(K,n,B,u,alpha,maxgen)
with open("%s/hybrid_phenos_%s.csv" %(data_dir,sim_id), "w") as f:
    writer = csv.writer(f)
    writer.writerows(offphenos01)

# i=0
# j=2
# offphenos = []
# for k in range(nHybrids):
#     # random parents
#     randpari = pop[i][np.random.choice(len(pop[i]))] 
#     randparj = pop[j][np.random.choice(len(pop[j]))]
#     # random parent phenotypes
#     phenpari = np.dot(randpari,mut[i]) 
#     phenparj = np.dot(randparj,mut[j])
#     # mutations held by random parents
#     mutpari = mut[i] * randpari[:,None]
#     mutparj = mut[j] * randparj[:,None]
#     A = set(tuple(x) for x in mutpari)
#     B = set(tuple(x) for x in mutparj)
#     # mutations shared by two parents (all in offspring)
#     sharedmuts = np.array([x for x in A & B])
#     # mutations not shared by two parents
#     unsharedmuts = np.array([x for x in A ^ B])
#     # which unshared mutations in offspring (free recombination between all loci, therefore gets each with 0.5 probability)
#     randmuts = np.random.randint(2, size=(len(unsharedmuts)))
#     offmuts = unsharedmuts * randmuts[:,None]
#     # offspring phenotype
#     offphenos.append(sum(np.append(sharedmuts,offmuts,axis=0)))

# offphenos02 = np.array(offphenos) #reformat correctly
# mean_offpheno02 = np.mean(offphenos,axis=0) #mean

# i=1
# j=2
# offphenos = []
# for k in range(nHybrids):
#     # random parents
#     randpari = pop[i][np.random.choice(len(pop[i]))] 
#     randparj = pop[j][np.random.choice(len(pop[j]))]
#     # random parent phenotypes
#     phenpari = np.dot(randpari,mut[i]) 
#     phenparj = np.dot(randparj,mut[j])
#     # mutations held by random parents
#     mutpari = mut[i] * randpari[:,None]
#     mutparj = mut[j] * randparj[:,None]
#     A = set(tuple(x) for x in mutpari)
#     B = set(tuple(x) for x in mutparj)
#     # mutations shared by two parents (all in offspring)
#     sharedmuts = np.array([x for x in A & B])
#     # mutations not shared by two parents
#     unsharedmuts = np.array([x for x in A ^ B])
#     # which unshared mutations in offspring (free recombination between all loci, therefore gets each with 0.5 probability)
#     randmuts = np.random.randint(2, size=(len(unsharedmuts)))
#     offmuts = unsharedmuts * randmuts[:,None]
#     # offspring phenotype
#     offphenos.append(sum(np.append(sharedmuts,offmuts,axis=0)))

# offphenos12 = np.array(offphenos) #reformat correctly
# mean_offpheno12 = np.mean(offphenos,axis=0) #mean

# ######################################################################
# ##PLOT##
# ######################################################################

# plot all phenotypes in parental populations
colors = cm.rainbow(np.linspace(0, 1, len(phenos)))
for i, c in zip(range(len(phenos)), colors):
    plt.scatter(phenos[i][:,0],phenos[i][:,1], color=c)

# plot mean phenotypes of parental population
plt.scatter(mean_pheno[:,0],mean_pheno[:,1], color='black')

# plot hybrid phenos
plt.scatter(offphenos01[:,0],offphenos01[:,1], color='gray')
plt.scatter(mean_offpheno01[0],mean_offpheno01[1], color='black')

# plt.scatter(offphenos02[:,0],offphenos02[:,1], color='green')
# plt.scatter(mean_offpheno02[0],mean_offpheno02[1], color='black')

# plt.scatter(offphenos12[:,0],offphenos12[:,1], color='gray')
# plt.scatter(mean_offpheno12[0],mean_offpheno12[1], color='black')

# show plot
# plt.show()

# save plot
plt.savefig('%s/plot_%s_hybrids.png' %(data_dir,sim_id))
