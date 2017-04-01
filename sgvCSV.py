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
maxgen = 10000

pop = []
mut = []
opt1s = [[0] * n, [0.5] * n, [-0.5] * n]
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
            pop.append(pickle.load(f, encoding = 'latin1'))
        except EOFError:
            break
    # load mut data
    g = open('%s/mut_%s.pkl' %(data_dir,sim_id), 'rb')
    while 1:
        try:
            mut.append(pickle.load(g))
        except EOFError:
            break

# make phenotype data
phenos = []
for i in range(len(opt1s)):
    phenos.append(np.dot(pop[i],mut[i]))


# save pheno data as CSV for R (can take a little time!)
import csv
with open("%s/phenos_%s.csv" %(data_dir,sim_id), "w") as f:
    writer = csv.writer(f)
    writer.writerows(phenos)