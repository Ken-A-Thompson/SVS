#Author: Matthew Osmond <mmosmond@zoology.ubc.ca>
#Description: Adaptation from standing genetic variance (SGV) in Fisher's geometric model, implications for hybrids
#Make hybrids from populations that have adapted

import numpy as np
import matplotlib.pyplot as plt
import pickle
import matplotlib.cm as cm
import csv
import itertools

######################################################################
##SGV (and DNM) or DNM##
######################################################################

# style = 'both' #standing genetic variation and de novo mutation
# style = 'sgv' #standing genetic variation only
style = 'dnm' #de novo mutation only

######################################################################
##ANCESTOR DATA##
######################################################################

K = 10000 #max number of parents (positive integer)
n = 2 #number of traits (positive integer)
B = 2 #number of offspring per generation per parent (positive integer)
u = 0.001 #mutation probability per genome (0<u<1)
alpha = 0.02 #mutation SD
sigma = 0.1
rep = 1

if style == 'both' or 'sgv':
#
    # which ancestral population
    maxgen = 5000 #SGV number of gens in burn-in (positive integer)
    sim_id = 'K%d_n%d_B%d_u%r_sigma%r_alpha%r_gens%r_burn_rep%d' %(K,n,B,u,sigma,alpha,maxgen,rep)
    data_dir = 'data'
#
    # load pop data
    f = open('%s/pop_%s.pkl' %(data_dir,sim_id), 'rb')
    popall = []
    while 1:
        try:
            popall.append(pickle.load(f))
        except EOFError:
            break

    # # load mut data
    # g = open('%s/mut_%s.pkl' %(data_dir,sim_id), 'rb')
    # mutall = []
    # while 1:
    #     try:
    #         mutall.append(pickle.load(g))
    #     except EOFError:
    #         break

if style == 'dnm':
#
    popall = np.array([[1]] * K)
    # mutfound = np.array([[0] * n])

# make csv of population list (mutations in each individual)
with open("%s/ancestor_pop_%s_%s.csv" %(data_dir,sim_id,style), "w") as f:
    writer = csv.writer(f)
    writer.writerows(popall) #write for all timepoints
    # writer.writerows(popall[-1]) #write for just last timepoint

######################################################################
##PARENT DATA##
######################################################################

maxgenAdapt = 2000 #number of generations during parent adaptation post-burnin (positive integer)
KAdapt = 2000 # carrying capacity of adapting populations (positive integer)
nfounders = KAdapt
rep = 1
uadapt = 0.001
alpha = 0.02
#sigma = 0.1

opt1s = [[0.100] * n, [0.101] * n] #which parental optima to use (Parallel Near)
# opt1s = [[0.100] * n, [-0.100] * n] #which parental optima to use (Parallel Near)
# opt1s = [[0.25] * n, [0.251] * n] #which parental optima to use (Parallel Near)
# opt1s = [[-0.249] * n, [0.251] * n] #which parental optima to use (Divergent)

data_dir = 'data'

pops = dict()
muts = dict()
for i in range(len(opt1s)):
#    
    # filename of data
    sim_id = 'K%d_n%d_B%d_u%r_alpha%r_gens%d_founders%d_opt%s_adapt_%s_rep%d' %(KAdapt,n,B,uadapt,alpha,maxgenAdapt,nfounders,'-'.join(str(e) for e in opt1s[i]),style,rep)
#
    # load pop data
    pop = []
    f = open('%s/pop_%s.pkl' %(data_dir,sim_id), 'rb')
    while 1:
        try:
            pop.append(pickle.load(f))
        except EOFError:
            break
    pops[i] = pop
#    
    # load mut data
    mut = []
    g = open('%s/mut_%s.pkl' %(data_dir,sim_id), 'rb')
    while 1:
        try:
            mut.append(pickle.load(g))
        except EOFError:
            break
    muts[i] = mut

# make csv of parent chromosomes
for i in range(len(pops)):
    sim_id = 'K%d_n%d_B%d_u%r_alpha%r_gens%r_founders%d_opt%s' %(KAdapt,n,B,u,alpha,maxgenAdapt,nfounders,'-'.join(str(e) for e in opt1s[i]))
    with open("%s/parent_pop_%s_%s.csv" %(data_dir,sim_id,style), "w") as f:
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

# make csv of parent phenotypes
for i in range(len(phenos)):
	sim_id = 'K%d_n%d_B%d_u%r_alpha%r_gens%r_opt%s' %(KAdapt,n,B,u,alpha,maxgenAdapt,'-'.join(str(e) for e in opt1s[i]))
	with open("%s/parent_phenos_%s_%s.csv" %(data_dir,sim_id,style), "w") as f:
		writer = csv.writer(f)
		writer.writerows(phenos[i])

#######################################################################
###MAKE HYBRIDS##
#######################################################################

nHybrids = 1000 #number of hybrids to make in crosses between each pair of parent populations

offphenos = dict()
mean_offphenos = dict()
pairs = list(itertools.combinations(range(len(pops)),2)) #all parent population combinations
for h in range(len(pairs)):
    i, j = pairs[h]
    offpheno = []
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
        offpheno.append(sum(np.append(sharedmuts,offmuts,axis=0)))
#
    offpheno = np.array(offpheno) #reformat correctly
    mean_offpheno = np.mean(offpheno,axis=0) #mean

    # save as csv
    sim_id = 'K%d_n%d_B%d_u%r_alpha%r_gens%r_opt1_%s_opt2_%s' %(KAdapt,n,B,u,alpha,maxgenAdapt,'-'.join(str(e) for e in opt1s[i]),'-'.join(str(e) for e in opt1s[j]))
    with open("%s/hybrid_phenos_%s_%s.csv" %(data_dir,sim_id,style), "w") as f:
        writer = csv.writer(f)
        writer.writerows(offpheno)

    #save for plotting in python
    # offphenos[h] = offpheno
    # mean_offphenos[h] = mean_offpheno

########################################################################
####PLOT##
########################################################################

# plot all phenotypes in parental populations
# for i in range(len(phenos)):
# 	colors = cm.rainbow(np.linspace(0, 1, len(phenos[i])))
# 	for j, c in zip(range(len(phenos[i])), colors):
# 		plt.scatter(phenos[i][j][:,0],phenos[i][j][:,1], color=c)

# # plot mean phenotypes of parental population
# for i in range(len(mean_phenos)):
#     plt.scatter(mean_phenos[i][:,0],mean_phenos[i][:,1], color='black')

# # plot hybrid phenos
# for i in range(len(pairs)):
#     plt.scatter(offphenos[i][:,0],offphenos[i][:,1], color='gray')
#     plt.scatter(mean_offphenos[i][0],mean_offphenos[i][1], color='black')

# # show plot
# plt.show()

# save plot
# plt.savefig('%s/plot_%s_hybrids.png' %(data_dir,sim_id))
