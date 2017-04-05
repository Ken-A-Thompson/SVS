#Author: Matthew Osmond <mmosmond@zoology.ubc.ca>
#Description: Adaptation from standing genetic variance (SGV) in Fisher's geometric model, implications for hybrids
#Make hybrids from populations that have adapted

import numpy as np
import matplotlib.pyplot as plt
import pickle
import matplotlib.cm as cm
import csv

######################################################################
##LOAD DATA FROM BURN IN##
######################################################################

K = 1000 #max number of parents (positive integer)
n = 2 #number of traits (positive integer)
B = 2 #number of offspring per generation per parent (positive integer)
u = 0.001 #mutation probability per genome (0<u<1)
alpha = 0.02 #mutation SD
maxgen = 10000 #number of gens in burn-in (positive integer)

sim_id = 'K%d_n%d_B%d_u%r_alpha%r_gens%r_burn' %(K,n,B,u,alpha,maxgen)

data_dir = 'data'

# load pop data
f = open('%s/pop_%s.pkl' %(data_dir,sim_id), 'rb')
popall = []
while 1:
    try:
        popall.append(pickle.load(f))
    except EOFError:
        break

# load mut data
g = open('%s/mut_%s.pkl' %(data_dir,sim_id), 'rb')
mutall = []
while 1:
    try:
        mutall.append(pickle.load(g))
    except EOFError:
        break

# make csv of last time point of population list (mutations in each individual)
with open("%s/ancestor_pop_%s.csv" %(data_dir,sim_id), "w") as f:
    writer = csv.writer(f)
    writer.writerows(popall[-1])

######################################################################
##LOAD DATA FROM PARENTAL POPNS##
######################################################################

maxgen = 10000 #number of generations during parent adaptation post-burnin (positive integer)
nfounders = len(popall[-1]) #number of individuals from burn-in used to found parntal population
opt1s = [[0.5] * n, [-0.5] * n] #simulations to load

data_dir = 'data'

pops = dict()
muts = dict()
for i in range(len(opt1s)):
    # filename of data
    sim_id = 'K%d_n%d_B%d_u%r_alpha%r_gens%r_founders%d_opt%s_adapt' %(K,n,B,u,alpha,maxgen,nfounders,'-'.join(str(e) for e in opt1s[i]))
    # load pop data
    pop = []
    f = open('%s/pop_%s.pkl' %(data_dir,sim_id), 'rb')
    while 1:
        try:
            pop.append(pickle.load(f))
        except EOFError:
            break
    pops[i] = pop
    # load mut data
    mut = []
    g = open('%s/mut_%s.pkl' %(data_dir,sim_id), 'rb')
    while 1:
        try:
            mut.append(pickle.load(g))
        except EOFError:
            break
    muts[i] = mut

# make csv
for i in range(len(pops)):
    sim_id = 'K%d_n%d_B%d_u%r_alpha%r_gens%r_opt%s' %(K,n,B,u,alpha,maxgen,'-'.join(str(e) for e in opt1s[i]))
    with open("%s/parent_pop_%s.csv" %(data_dir,sim_id), "w") as f:
        writer = csv.writer(f)
        writer.writerows(pops[i])

######################################################################
##PARENTAL PHENOTYPES##
######################################################################

#Make parental phenotype data
phenos = dict()
for i in range(len(pops)):
	pheno = []
	for j in range(len(pops[i])):
		pheno.append(np.dot(pops[i][j],muts[i][j]))
	phenos[i] = pheno

# mean parental phenotypes
mean_phenos = dict()
for i in range(len(phenos)):
    mean_phenos[i] = np.mean(phenos[i], axis=1)

# save pheno data as CSV for R (can take a little time!)
for i in range(len(phenos)):
	sim_id = 'K%d_n%d_B%d_u%r_alpha%r_gens%r_opt%s' %(K,n,B,u,alpha,maxgen,'-'.join(str(e) for e in opt1s[i]))
	with open("%s/parent_phenos_%s.csv" %(data_dir,sim_id), "w") as f:
		writer = csv.writer(f)
		writer.writerows(phenos[i])

# ######################################################################
# ##MAKE HYBRIDS##
# ######################################################################

nHybrids = 1000 #number of hybrids to make

i=0
j=1
offphenos = []
for k in range(nHybrids):
    # random parents
    randpari = pops[i][-1][np.random.choice(len(pops[i][-1]))] 
    randparj = pops[j][-1][np.random.choice(len(pops[j][-1]))]
    # random parent phenotypes
    phenpari = np.dot(randpari,muts[i][-1]) 
    phenparj = np.dot(randparj,muts[j][-1])
    # mutations held by random parents
    mutpari = muts[i][-1] * randpari[:,None]
    mutparj = muts[j][-1] * randparj[:,None]
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
for i in range(len(phenos)):
	colors = cm.rainbow(np.linspace(0, 1, len(phenos[i])))
	for j, c in zip(range(len(phenos[i])), colors):
		plt.scatter(phenos[i][j][:,0],phenos[i][j][:,1], color=c)

# plot mean phenotypes of parental population
for i in range(len(phenos)):
    plt.scatter(mean_phenos[i][:,0],mean_phenos[i][:,1], color='black')

# plot hybrid phenos
plt.scatter(offphenos01[:,0],offphenos01[:,1], color='gray')
plt.scatter(mean_offpheno01[0],mean_offpheno01[1], color='black')

# # plt.scatter(offphenos02[:,0],offphenos02[:,1], color='green')
# # plt.scatter(mean_offpheno02[0],mean_offpheno02[1], color='black')

# # plt.scatter(offphenos12[:,0],offphenos12[:,1], color='gray')
# # plt.scatter(mean_offpheno12[0],mean_offpheno12[1], color='black')

# # show plot
# plt.show()

# # save plot
plt.savefig('%s/plot_%s_hybrids.png' %(data_dir,sim_id))
