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

    outfile_A = open("dat_%s.dat" %(sim_id),"w")

    return [outfile_A]

def write_data_to_output(fileHandles, gen, data):
    """
    This function writes a (time, data) pair to the
    corresponding output file. We write densities
    not abundances.
    """
    
    for i in range(0,len(fileHandles)):
        fileHandles[i].write("%d  %r\n" %(gen, data))
    
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
maxgen = 10000 #maximum number of generations (positive integer)
opt = [0] * n #optimum phenotype
outputFreq = 100 #record and print update this many generations

######################################################################
##SIMULATION##
######################################################################

def main():

	# intialize
	gen = 0 #generation
	pop = np.array([[0]] * N0]) #list of mutations held by each individual (0 is mutation that does nothing - ie a placeholder)
	mut = np.array([[0] * n]) #list of phenotypic effect of mutations
		
	# open output files
	fileHandles = open_output_files(K, n, B, u, alpha) 
	
	while gen < maxgen:

		# genotype to phenotype
		phenos = np.array([sum(mut[i]) for i in pop]) #add mutation effects up within each individual (try to improve by getting rid of loop!)

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
		newmuts = np.random.normal(0, alpha, (nmuts,n)) #phenotypic effect of new mutations

		# update
		pop = 
		mut = np.append(mut,newmuts,axis=0) #append effect of new mutations to mutation list

		#end simulation if extinct        
		if len(pop) == 0: 
			print("Extinct")              
			break 
            
        #otherwise continue
        # dump data every outputFreq iteration
        # also print a short progess message 
		if (gen % outputFreq) == 0:
			write_data_to_output(fileHandles, gen, pop)
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

























