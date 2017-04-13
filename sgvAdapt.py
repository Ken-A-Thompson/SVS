#Author: Matthew Osmond <mmosmond@zoology.ubc.ca>
#Description: Adaptation from standing genetic variance (SGV) in Fisher's geometric model, implications for hybrids
#Adapt using standing genetic variance from burn-in as well as new mutations

import numpy as np
import time
import pickle

######################################################################
##HELPER FUNCTIONS##
######################################################################

def open_output_files(K, n, B, u, alpha, maxgenAdapt, nfounders, opt1, style):
	"""
	This function opens the output files and returns file
	handles to each.
	"""

	sim_id = 'K%d_n%d_B%d_u%r_alpha%r_gens%d_founders%d_opt%s_adapt_%s' %(KAdapt,n,B,u,alpha,maxgenAdapt,nfounders,'-'.join(str(e) for e in opt1),style)
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
	#   pickle.dump(data[i],fileHandles[i])

def close_output_files(fileHandles):
	"""
	This function closes all output files.
	"""

	for i in range(0,len(fileHandles)):
		fileHandles[i].close()


######################################################################
##PARAMETERS FOR ADAPTING POPULATIONS##
######################################################################

maxgenAdapt = 1000 #maximum number of generations (positive integer)
KAdapt = 1000 # maximum population size of adapting populations

outputFreq = 100 #record and print update this many generations

remove_lost = True #If true, remove mutations that are lost (0 for all individuals)
# remove_lost = False
#if remove_lost = True, remove...
# remove = 'any' #... any mutation that is lost
remove = 'derived' #.. any derived (not from ancestor) mutation that is lost 

style = 'sgv' #standing genetic variance and de novo mutation
# style = 'dnm' #de novo mutation only

######################################################################
##SGV (and DNM)##
######################################################################

if style == 'sgv':

	# which ancestor (burn-in) to get data from
	K = 1000 #max number of parents (positive integer)
	n = 2 #number of traits (positive integer)
	B = 2 #number of offspring per generation per parent (positive integer)
	u = 0.001 #mutation probability per genome (0<u<1)
	alpha = 0.02 #mutation SD
	maxgen = 1000 #SGV number of gens in burn-in (positive integer)

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

	# choose founding population from ancestor
	nfounders = min(KAdapt,len(popall[-1])) #size of founding population
	# nfounders = len(popall[-1])
	whofounds = np.random.choice(len(popall[-1]),size=nfounders) #random choice of nfounders from ancestral population 
	popfound = popall[-1][whofounds] #list of mutations held by each founding individual
	if remove_lost and remove == 'any': #if removing ancestral mutations when lost
		mutfound = np.delete(mutall[-1], np.where(~popfound.any(axis=0))[0], axis=0) #remove lost mutation effects (from random choice of founders)
		popfound = popfound[:, ~np.all(popfound==0, axis=0)] #remove lost loci
	else:
		mutfound = mutall[-1]

######################################################################
##DNM only##
######################################################################

if style == 'dnm':

	#choose parameters
	n = 2 #number of traits (positive integer)
	B = 2 #number of offspring per generation per parent (positive integer)
	u = 0.001 #mutation probability per genome (0<u<1)
	alpha = 0.02 #mutation SD

	nfounders = KAdapt #initial population size
	popfound = np.array([[1]] * nfounders) #list of mutations held by each individual (all start with same mutation at first locus)
	mutfound = np.array([[0] * n]) #list of phenotypic effect of initial mutations (mutation at first locus puts all individuals at origin)

######################################################################
##NEW OPTIMUM##
######################################################################

opt1 = [0] * n #optimum phenotype 
# opt1 = [0.51] * n #optimum phenotype 
# opt1 = [-0.51] * n #optimum phenotype 
# opt1 = [0.49] * n #optimum phenotype 
# opt1 = [-0.5] * n #optimum phenotype 

######################################################################
##SIMULATION##
######################################################################

def main():

	pop = popfound
	mut = mutfound
		
	# open output files
	fileHandles = open_output_files(KAdapt, n, B, u, alpha, maxgenAdapt, nfounders, opt1, style) 
	
	# optimum phenotype
	opt = opt1

	gen = 1 #generation
	while gen < maxgenAdapt + 1:

		# genotype to phenotype
		phenos = np.dot(pop,mut) #sum mutations held by each individual

		# viability selection
		dist = np.linalg.norm(phenos - opt, axis=1) #phenotypic distance from optimum
		w = np.exp(-dist**2) #probability of survival
		rand1 = np.random.uniform(size = len(pop)) #random uniform number in [0,1] for each individual
		surv = pop[rand1 < w] #survivors
		if len(surv) > KAdapt:
			surv = surv[np.random.randint(len(surv), size = KAdapt)] #randomly choose KAdapt individuals if more than KAdapt
		
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

		# remove lost mutations (all zero columns in pop)
		if remove_lost:
			if remove == 'any':
				keep = pop.any(axis=0)
				mut = mut[keep]
				pop = pop[:, keep]
			if remove == 'derived':
				keep = pop.any(axis=0) | np.array(range(len(mut))) < len(mutfound) #keep mutation if not lost or ancestral 
				mut = mut[keep]
				pop = pop[:, keep]
		
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
