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
##HELPER FUNCTIONS##
######################################################################

def open_output_files(K_adapt, n_adapt, B_adapt, u_adapt, alpha_adapt, sigma_adapt, maxgen_adapt, opt1s, style, style2):
    """
    This function opens the output files and returns file
    handles to each.
    """

    sim_id = 'K%d_n%d_B%d_u%r_alpha%r_sigma%r_gens%r_opt1_%s_opt2_%s_%s_%s' %(K_adapt, n_adapt, B_adapt, u_adapt, alpha_adapt, sigma_adapt, maxgen_adapt,'-'.join(str(e) for e in opt1s[i]),'-'.join(str(e) for e in opt1s[j]), style, style2)
    data_dir = 'data'

    outfile_A = open("%s/hybrids_%s.pkl" %(data_dir,sim_id),"wb")

    return [outfile_A]

def write_data_to_output(fileHandles, data):
    """
    This function writes a (time, data) pair to the
    corresponding output file. We write densities
    not abundances.
    """
    
    for i in range(0,len(fileHandles)):
        pickle.dump(data[i],fileHandles[i])

    # for i in range(0,len(fileHandles)):
    #   pickle.dump(data[i],fileHandles[i])

def close_output_files(fileHandles):
    """
    This function closes all output files.
    """

    for i in range(0,len(fileHandles)):
        fileHandles[i].close()

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

save_CSV = True #save CSV of hybrids?

######################################################################
##PARENTAL PARAMETERS##
######################################################################

K_adapt = 1000 #max number of parents (positive integer)
n_adapt = 2 #number of traits (positive integer)
B_adapt = 2 #number of offspring per generation per parent (positive integer)
u_adapt = 0.001 #mutation probability per genome (0<u<1) (set to zero for sgv only)
alpha_adapt = 0.02 #mutational sd (positive real number)
sigma_adapt = 10 #strength of selection (positive real number)
theta_adapt = 0 #drift parameter in Ornstein-Uhlenbeck movement of the phenotypic optimum (positive real number; setting to zero makes Brownian motion)
sigma_opt_adapt = 0 #diffusion parameter in Ornstein-Uhlenbeck movement of the phenotypic optimum (positive real number; setting to zero makes constant optimum at opt0)
nfounders = K_adapt

#meta-parameters
maxgen_adapt = 1000 #maximum number of generations (positive integer)
rep = 1 #which replicate? (positive integer)

#change in optima that parental populations experienced
opt1s = [[0.1] * n_adapt, [0.101] * n_adapt] #which parental optima to use (Parallel)
# opt1s = [[0.1] * n_adapt, [-0.1] * n_adapt] #which parental optima to use (Divergent)

######################################################################
##PARENTAL DATA##
######################################################################

#get parental data
pops = dict()
muts = dict()
for i in range(len(opt1s)):
    
    # filename of data
    sim_id = 'K%d_n%d_B%d_u%r_alpha%r_sigma%r_theta%r_sigmaopt%r_gens%d_founders%d_opt%s_adapt_%s_%s_rep%d' %(K_adapt, n_adapt, B_adapt, u_adapt, alpha_adapt, sigma_adapt, theta_adapt, sigma_opt_adapt, maxgen_adapt, nfounders,'-'.join(str(e) for e in opt1s[i]), style, style2, rep)
    
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

#######################################################################
###MAKE HYBRIDS##
#######################################################################

offphenos = dict()
# mean_offphenos = dict()
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
            offmuts = np.array([[0]*n_adapt]) #give something in case empty
        # offspring phenotype
        offpheno.append(sum(np.append(sharedmuts,offmuts,axis=0)))

    offpheno = np.array(offpheno) #reformat correctly
    # mean_offpheno = np.mean(offpheno,axis=0) #mean

    fileHandles = open_output_files(K_adapt, n_adapt, B_adapt, u_adapt, alpha_adapt, sigma_adapt, maxgen_adapt, opt1s, style, style2) 
    write_data_to_output(fileHandles, [offpheno])
    close_output_files(fileHandles)

if save_CSV == True:

    # save as csv
    sim_id = 'K%d_n%d_B%d_u%r_alpha%r_sigma%r_gens%r_opt1_%s_opt2_%s_%s_%s' %(K_adapt, n_adapt, B_adapt, u_adapt, alpha_adapt, sigma_adapt, maxgen_adapt,'-'.join(str(e) for e in opt1s[i]),'-'.join(str(e) for e in opt1s[j]), style, style2)
    with open("%s/hybrid_%s.csv" %(data_dir,sim_id), "w") as f:
        writer = csv.writer(f)
        writer.writerows(offpheno)
