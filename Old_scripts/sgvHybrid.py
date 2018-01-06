#Author: Matthew Osmond <mmosmond@zoology.ubc.ca> & Ken A. Thompson <ken.thompson@zoology.ubc.ca>
#Description: Adaptation from standing genetic variance (SGV) in Fisher's geometric model, implications for hybrids
#Make hybrids from populations that have adapted

import numpy as np
import matplotlib.pyplot as plt
import pickle
import matplotlib.cm as cm
import csv
import itertools

######################################################################
##PARAMETERS##
######################################################################

nHybrids = 1000 #number of hybrids to make in crosses between each pair of parent populations

#style of parental adaptation
style = 'both' #standing genetic variation and de novo mutation (if u_adapt > 0)
# style = 'sgv' #standing genetic variation only (no mutation)
# style = 'dnm' #de novo mutation only

#style of ancestor creation
# style2 = 'burn' #ancestor created by burn-in
style2 = 'art' #artificially created ancestor

data_dir = 'data' #where parental and common ancestor data is, and where hybrid data will be deposited

save_CSVs = False #save CSVs?

######################################################################
##PARENTAL PARAMETERS##
######################################################################

K_adapt = 1000 #max number of parents (positive integer)
n_adapt = 2 #number of traits (positive integer)
B_adapt = 2 #number of offspring per generation per parent (positive integer)
u_adapt = 0.01 #mutation probability per genome (0<u<1) (set to zero for sgv only)
alpha_adapt = 0.02 #mutational sd (positive real number)
sigma_adapt = 1 #strength of selection (positive real number)
theta_adapt = 0 #drift parameter in Ornstein-Uhlenbeck movement of the phenotypic optimum (positive real number; setting to zero makes Brownian motion)
sigma_opt_adapt = 0 #diffusion parameter in Ornstein-Uhlenbeck movement of the phenotypic optimum (positive real number; setting to zero makes constant optimum at opt0)

#meta-parameters
maxgen_adapt = 1000 #maximum number of generations (positive integer)
rep = 1 #which replicate? (positive integer)

#change in optima that parental populations experienced
# opt1s = [[0.1] * n, [0.1] * n] #which parental optima to use (Parallel)
opt1s = [[0.1] * n, [-0.1] * n] #which parental optima to use (Divergent)

######################################################################
##COMMON ANCESTOR PARAMETERS##
######################################################################

if style == 'sgv' or 'both':

    K = K_adapt #max number of parents (positive integer)
    n = n_adapt #number of traits (positive integer)
    alpha = alpha_adapt #mutational sd (positive real number)
    rep = 1 #which replicate (positive integer)

    if style2 == 'real':

        B = B_adapt #number of offspring per generation per parent (positive integer)
        u = u_adapt #mutation probability per genome (0<u<1)
        sigma = sigma_adapt #strength of selection (positive real number)
        theta = 1 #drift parameter in Ornstein-Uhlenbeck movement of the phenotypic optimum (positive real number; setting to zero makes Brownian motion)
        sigma_opt = 0.1 #diffusion parameter in Ornstein-Uhlenbeck movement of the phenotypic optimum (positive real number; setting to zero makes constant optimum at opt0)
        maxgen = 1000 #number of gens in ancestral burn-in (positive integer)

        #file to load
        sim_id = 'K%d_n%d_B%d_u%r_alpha%r_sigma%r_theta_%r_sigmaopt%r_gens%r_burn_rep%d' %(K,n,B,u,alpha,sigma,theta,sigma_opt,maxgen,rep)

    if style2 == 'art':

        n_muts = 1000 #number of mutations (positive integer)
        p_mut = 0.5 #probability of having mutation at any one locus (0<p<1)
        
        #file to load
        sim_id = 'K%d_n%d_alpha%r_nmuts%r_pmut%r_create_rep%d' %(K, n, alpha, n_muts, p_mut, rep)


######################################################################
##COMMON ANCESTOR DATA##
######################################################################

if style == 'both' or 'sgv':

    # load pop data
    f = open('%s/pop_%s.pkl' %(data_dir,sim_id), 'rb')
    popall = []
    while 1:
        try:
            popall.append(pickle.load(f))
        except EOFError:
            break
    nfounders = min(K_adapt,len(popall[-1]))

    # load mut data
    g = open('%s/mut_%s.pkl' %(data_dir,sim_id), 'rb')
    mutall = []
    while 1:
        try:
            mutall.append(pickle.load(g))
        except EOFError:
            break

if style == 'dnm':

    popall = np.array([[1]] * K)
    # mutfound = np.array([[0] * n])
    # nfounders = min(K_adapt,len(popall))

######################################################################
##PARENTAL DATA##
######################################################################

#get parental data
pops = dict()
muts = dict()
for i in range(len(opt1s)):
    
    # filename of data
    sim_id = 'K%d_n%d_B%d_u%r_alpha%r_sigma%r_theta%r_sigmaopt%r_gens%d_founders%d_opt%s_adapt_%s_rep%d' %(K_adapt, n_adapt, B_adapt, u_adapt, alpha_adapt, sigma_adapt, theta_adapt, sigma_opt_adapt, maxgen_adapt, nfounders,'-'.join(str(e) for e in opt1s[i]), style, rep)
    
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

if save_CSVs == True:
    
    # make csv of parent chromosomes
    for i in range(len(pops)):
        sim_id = 'K%d_n%d_B%d_u%r_alpha%r_sigma%r_theta%r_sigmaopt%r_gens%d_founders%d_opt%s_adapt_%s_rep%d' %(K_adapt, n_adapt, B_adapt, u_adapt, alpha_adapt, sigma_adapt, theta_adapt, sigma_opt_adapt, maxgen_adapt, nfounders,'-'.join(str(e) for e in opt1s[i]), style, rep)
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
    mean_phenos[i] = [np.mean(j) for j in phenos[i]]

if save_CSVs == True:

    # make csv of parent phenotypes
    for i in range(len(phenos)):
        sim_id = 'K%d_n%d_B%d_u%r_alpha%r_sigma%r_theta%r_sigmaopt%r_gens%d_founders%d_opt%s_adapt_%s_rep%d' %(K_adapt, n_adapt, B_adapt, u_adapt, alpha_adapt, sigma_adapt, theta_adapt, sigma_opt_adapt, maxgen_adapt, nfounders,'-'.join(str(e) for e in opt1s[i]), style, rep)
        with open("%s/parent_phenos_%s.csv" %(data_dir,sim_id), "w") as f:
            writer = csv.writer(f)
            writer.writerows(phenos[i])

#######################################################################
###MAKE HYBRIDS##
#######################################################################

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

print(mean_offpheno)

if save_CSVs == True:

    # save as csv
    sim_id = 'K%d_n%d_B%d_u%r_alpha%r_gens%r_opt1_%s_opt2_%s' %(K_adapt,n,B,u,alpha,maxgen_adapt,'-'.join(str(e) for e in opt1s[i]),'-'.join(str(e) for e in opt1s[j]))
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
