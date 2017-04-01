#Author: Matthew Osmond <mmosmond@zoology.ubc.ca>
#Description: Adaptation from standing genetic variance (SGV) in Fisher's geometric model, implications for hybrids
#Adapt using standing genetic variance from burn-in as well as new mutations

import numpy as np
import time
import pickle

######################################################################
##LOAD DATA FROM BURN_IN##
######################################################################

K = 1000 #max number of parents (positive integer)
n = 2 #number of traits (positive integer)
B = 2 #number of offspring per generation per parent (positive integer)
u = 0.001 #mutation probability per genome (0<u<1)
alpha = 0.01 #mutation SD
maxgen = 1000 #number of gens in burn-in

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

######################################################################
##HELPER FUNCTIONS##
######################################################################

def open_output_files(K, n, B, u, alpha, maxgen, opt1):
    """
    This function opens the output files and returns file
    handles to each.
    """

    sim_id = 'K%d_n%d_B%d_u%r_alpha%r_gens%r_opt%s_adapt' %(K,n,B,u,alpha,maxgen,'-'.join(str(e) for e in opt1))
    data_dir = 'data'

    outfile_A = open("%s/pop_%s.pkl" %(data_dir,sim_id),"ab")
    outfile_B = open("%s/mut_%s.pkl" %(data_dir,sim_id),"ab")
    outfile_C = open("%s/gen_%s.pkl" %(data_dir,sim_id),"ab")

    return [outfile_A, outfile_B, outfile_C]

def write_data_to_output(fileHandles, data):
    """
    This function writes a (time, data) pair to the
    corresponding output file. We write densities
    not abundances.
    """
	
    for i in range(0,len(fileHandles)):
        pickle.dump(data[i],fileHandles[i])

	# for i in range(0,len(fileHandles)):
	# 	pickle.dump(data[i],fileHandles[i])

def close_output_files(fileHandles):
    """
    This function closes all output files.
    """

    for i in range(0,len(fileHandles)):
        fileHandles[i].close()

######################################################################
##PARAMETERS##
######################################################################

maxgen = 1000 #maximum number of generations (positive integer)
# opt1 = [0] * n #optimum phenotype 
# opt1 = [0.5] * n #optimum phenotype 
opt1 = [-0.5] * n #optimum phenotype 

outputFreq = 100 #record and print update this many generations

remove_lost = True #remove mutations that are lost?

######################################################################
##SIMULATION##
######################################################################

def main():

	# intialize
	pop = popall[-1] #list of mutations held by each individual (all start with same mutation at first locus)
	mut = mutall[-1] #list of phenotypic effect of initial mutations (mutation at first locus puts all individuals at origin)
		
	# open output files
	fileHandles = open_output_files(K, n, B, u, alpha, maxgen, opt1) 
	
	# optimum phenotype
	opt = opt1

	gen = 1 #generation
	while gen < maxgen + 1:

		# genotype to phenotype
		phenos = np.dot(pop,mut) #sum mutations held by each individual

		# viability selection
		dist = np.linalg.norm(phenos - opt, axis=1) #phenotypic distance from optimum
		w = np.exp(-dist**2) #probability of survival
		rand1 = np.random.uniform(size = len(pop)) #random uniform number in [0,1] for each individual
		surv = pop[rand1 < w] #survivors
		if len(surv) > K:
			surv = surv[np.random.randint(len(surv), size = K)] #randomly choose K individuals if more than K
		
		#end simulation if extinct        
		if len(surv) == 0: 
			print("Extinct")              
			break 
            
        #otherwise continue
        # dump data every outputFreq iteration
        # also print a short progess message (generation and number of parents)
		if gen > 0 and (gen % outputFreq) == 0:
			write_data_to_output(fileHandles, [surv,mut,gen])
			print("gen %d    N %d" %(gen, len(surv)))  
			
		# birth
		# off = np.repeat(surv, B, axis=0) #offspring of survivors (asexual)
		# sex: i.e., make diploid from random haploid parents then segregate to haploid offspring
		# pairs = np.transpose(np.array([np.arange(len(surv)),np.random.randint(len(surv),size=len(surv))])) #random mate pairs (can mate multiple times; all mate)
		pairs = np.resize(np.random.choice(len(surv), size=len(surv), replace=False), (int(len(surv)/2), 2)) #random mate pairs (each mates at most once and not with self)
		rand2 = np.random.randint(2, size=(len(pairs), len(surv[0]))) #from which parent each offspring inherits each allele (free recombination, fair transmission)
		rec = np.resize(np.append(rand2,1-rand2,axis=1),(len(rand2),2,len(rand2[0]))) #reshape
		off1 = np.sum(surv[pairs] * rec, axis=1) #one product of meiosis
		off2 = np.sum(surv[pairs] * (1-rec), axis=1) #other product of meiosis
		off = np.repeat(np.append(off1, off2, axis=0), B, axis=0) #each product of meiosis produced B times

		# mutation
		rand3 = np.random.uniform(size = len(off)) #random uniform number in [0,1] for each offspring
		nmuts = sum(rand3 < u) # mutate if random number is below mutation rate; returns number of new mutations
		whomuts = np.where(rand3 < u) #indices of mutants
		newmuts = np.random.normal(0, alpha, (nmuts,n)) #phenotypic effect of new mutations

		# update
		pop = np.append(off, np.transpose(np.identity(len(off),dtype=int)[whomuts[0]]), axis=1) #add new loci and identify mutants
		mut = np.append(mut,newmuts,axis=0) #append effect of new mutations to mutation list

		# remove lost mutations (all zero columns)
		if remove_lost:
			mut = np.delete(mut, np.where(~pop.any(axis=0))[0], axis=0)
			pop = pop[:, ~np.all(pop==0, axis=0)] 
		
		# go to next generation
		gen += 1

	# cleanup
	close_output_files(fileHandles)

######################################################################
##RUNNING##
######################################################################    
    
#run (with timer)
start = time.time()
main()
end = time.time()
print(end-start)
