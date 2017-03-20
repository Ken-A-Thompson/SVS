#Author: Matthew Osmond <mmosmond@zoology.ubc.ca>
#Description: Adaptation from standing genetic variance (SGV) in Fisher's geometric model, implications for hybrids

import numpy as np
import time

######################################################################
##HELPER FUNCTIONS##
######################################################################

def open_output_files(K, n, B, u, alpha):
    """
    This function opens the output files and returns file
    handles to each.
    """

    sim_id = 'K%d_n%d_B%d_u%r_alpha%r' %(K,n,B,u,alpha)

    outfile_A = open("pop_%s.dat" %(sim_id),"w")
    outfile_B = open("mut_%s.dat" %(sim_id),"w")

    return [outfile_A, outfile_B]

def write_data_to_output(fileHandles, gen, data):
    """
    This function writes a (time, data) pair to the
    corresponding output file. We write densities
    not abundances.
    """
    
    for i in range(0,len(fileHandles)):
        fileHandles[i].write("%d  %r\n" %(gen, data[i]))
    
def close_output_files(fileHandles):
    """
    This function closes all output files.
    """

    for i in range(0,len(fileHandles)):
        fileHandles[i].close()

######################################################################
##PARAMETERS##
######################################################################

K = 1000 #max population size (positive integer)
n = 2 #number of traits (positive integer)
B = 2 #number of offspring per generation per surviving individual
u = 0.01 #mutation probability
alpha = 0.1 #mutational sd

N0 = K #initial population size
maxgen = 1000 #maximum number of generations (positive integer)
opt = [0] * n #optimum phenotype
outputFreq = 100 #record and print update this many generations

######################################################################
##SIMULATION##
######################################################################

def main():

	# intialize
	pop = np.array([[1]] * N0) #list of mutations held by each individual (all start with same mutation at first locus)
	mut = np.array([[0] * n]) #list of phenotypic effect of initial mutations (mutation at first locus puts all individuals at origin)
		
	# open output files
	fileHandles = open_output_files(K, n, B, u, alpha) 
	
	gen = 0 #generation
	while gen < maxgen:

		# genotype to phenotype
		phenos = np.dot(pop,mut)

		# viability selection
		dist = np.linalg.norm(phenos - opt, axis=1) #phenotypic distance from optimum
		w = np.exp(-dist**2) #probability of survival
		rand = np.random.uniform(size = len(pop)) #random uniform number in [0,1] for each individual
		surv = pop[rand < w] #survivors
		if len(surv) > K:
			surv = surv[np.random.randint(len(surv), size = K)] #randomly choose K individuals if more than K
		
		# birth
		off = np.repeat(surv, B, axis=0) #offspring of survivors (asexual)
		
		# mutation
		rand = np.random.uniform(size = len(off)) #random uniform number in [0,1] for each offspring
		nmuts = sum(rand < u) #number of new mutations
		whomuts = np.where(rand < u) #indices of mutants
		newmuts = np.random.normal(0, alpha, (nmuts,n)) #phenotypic effect of new mutations

		# update
		pop = np.append(off, np.transpose(np.identity(len(off),dtype=int)[whomuts[0]]), axis=1) #add new loci and identify mutants
		mut = np.append(mut,newmuts,axis=0) #append effect of new mutations to mutation list

		# remove lost mutations (all zero columns)
		mut = np.delete(mut, np.where(~pop.any(axis=0))[0], axis=0)
		pop = pop[:, ~np.all(pop==0, axis=0)]

		#end simulation if extinct        
		if len(pop) == 0: 
			print("Extinct")              
			break 
            
        #otherwise continue
        # dump data every outputFreq iteration
        # also print a short progess message 
		if (gen % outputFreq) == 0:
			write_data_to_output(fileHandles, gen, [pop,mut])
			print("gen %d    N %d" %(gen, len(pop)/B))   

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

























