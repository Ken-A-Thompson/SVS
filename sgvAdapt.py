#Author: Matthew Osmond <mmosmond@zoology.ubc.ca> & Ken A. Thompson <ken.thompson@zoology.ubc.ca>

#Description: Adaptation from standing genetic variance (SGV) in Fisher's geometric model, implications for hybrids
#Adapt using standing genetic variance from burn-in as well as new mutations

import numpy as np
import time
import pickle

######################################################################
##HELPER FUNCTIONS##
######################################################################

def open_output_files(K, n, B, u, alpha, maxgenAdapt, nfounders, opt1, style, rep):
	"""
	This function opens the output files and returns file
	handles to each.
	"""

	sim_id = 'K%d_n%d_B%d_u%r_alpha%r_gens%d_founders%d_opt%s_adapt_%s_rep%d' %(KAdapt,n,B,u,alpha,maxgenAdapt,nfounders,'-'.join(str(e) for e in opt1),style,rep)
	data_dir = 'data'

	outfile_A = open("%s/pop_%s.pkl" %(data_dir,sim_id),"wb")
	outfile_B = open("%s/mut_%s.pkl" %(data_dir,sim_id),"wb")
	outfile_C = open("%s/gen_%s.pkl" %(data_dir,sim_id),"wb")

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

maxgenAdapt = 2000 #maximum number of generations (positive integer)
KAdapt = 1000 # maximum population size of adapting populations
sigmaAdapt = 2 #strength of selection

outputFreq = 500 #record and print update this many generations

remove_lost = True #If true, remove mutations that are lost (0 for all individuals)
# remove_lost = False

remove = 'derived' #.. any derived (not from ancestor) mutation that is lost 

######################################################################
##SIMULATION GENETICS##
######################################################################

# style = 'both' #standing genetic variation and de novo mutation
# style = 'sgv' #standing genetic variation only
style = 'dnm' #de novo mutation only

nReps = 1

######################################################################
##SGV (and DNM)##
######################################################################

if style == 'both':
#
	# which ancestor (burn-in) to get data from
	K = 10000 #max number of parents (positive integer)
	n = 2 #number of traits (positive integer)
	B = 2 #number of offspring per generation per parent (positive integer)
	u = 0.001 #mutation probability per genome (0<u<1)
	sigma = 0.1
	alpha = 0.02 #mutation SD
	maxgen = 5000 #number of gens in ancestral burn-in (positive integer)
	rep = 1
#
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
#
	# load mut data
	g = open('%s/mut_%s.pkl' %(data_dir,sim_id), 'rb')
	mutall = []
	while 1:
		try:
			mutall.append(pickle.load(g))
		except EOFError:
			break
#
	# choose founding population from ancestor
	nfounders = min(KAdapt,len(popall[-1])) #size of founding population
	# nfounders = len(popall[-1])
	whofounds = np.random.choice(len(popall[-1]),size=nfounders) #random choice of nfounders from ancestral population 
	popfound = popall[-1][whofounds] #list of mutations held by each founding individual
	if remove_lost and remove == 'any': #if removing ancestral mutations when lost
		keep = popfound.any(axis=0)
		mutfound = mutall[-1][keep]
		popfound = popfound[:, keep]
	else:
		mutfound = mutall[-1]

######################################################################
##SGV only##
######################################################################

if style == 'sgv':
#
	# which ancestor (burn-in) to get data from
	K = 10000 #max number of parents (positive integer)
	n = 2 #number of traits (positive integer)
	B = 2 #number of offspring per generation per parent (positive integer)
	u = 0 #mutation probability per genome (0 because SGV only)
	uburn = 0.001 #ancestral mutation probability per genome (0<u<1)
	sigma = 0.1 #ancestral strength of selection
	alpha = 0.02 #mutation SD
	alphaburn = 0.02 #ancestral mutation SD
	maxgen = 5000 #SGV number of gens in burn-in (positive integer)
	rep = 1


#
	sim_id = 'K%d_n%d_B%d_u%r_alpha%r_gens%r_burn_rep%d' %(K,n,B,uburn,alphaburn,maxgen,rep)
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
#
	# load mut data
	g = open('%s/mut_%s.pkl' %(data_dir,sim_id), 'rb')
	mutall = []
	while 1:
		try:
			mutall.append(pickle.load(g))
		except EOFError:
			break
#
	# choose founding population from ancestor
	nfounders = min(KAdapt,len(popall[-1])) #size of founding population
	# nfounders = len(popall[-1])
	whofounds = np.random.choice(len(popall[-1]),size=nfounders) #random choice of nfounders from ancestral population 
	popfound = popall[-1][whofounds] #list of mutations held by each founding individual
	if remove_lost and remove == 'any': #if removing ancestral mutations when lost
		keep = popfound.any(axis=0)
		mutfound = mutall[-1][keep]
		popfound = popfound[:, keep]
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

# opt1 = [0.1] * n #optimum phenotype 
# opt1 = [0.249] * n #optimum phenotype 
# opt1 = [0.25] * n #optimum phenotype
# opt1 = [0.251] * n #optimum phenotype
opt1 = [0.100] * n #optimum phenotype 
# opt1 = [0.101] * n #optimum phenotype 
# opt1 = [-0.100] * n #optimum phenotype


######################################################################
##SIMULATION##
######################################################################

def main():

	rep = 1
	while rep < nReps + 1:

		pop = popfound
		mut = mutfound
			
		# open output files
		fileHandles = open_output_files(KAdapt, n, B, u, alpha, maxgenAdapt, nfounders, opt1, style, rep) 
		
		# optimum phenotype
		opt = opt1

		gen = 1 #generation
		while gen < maxgenAdapt + 1:

			# genotype to phenotype
			phenos = np.dot(pop,mut) #sum mutations held by each individual

			# viability selection
			dist = np.linalg.norm(phenos - opt, axis=1) #phenotypic distance from optimum
			w = np.exp(-sigmaAdapt*dist**2) #probability of survival
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
				print("rep %d    gen %d    N %d" %(rep, gen, len(surv)))  
				
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
					segregating = pop.any(axis=0)
					ancestral = np.array(range(len(mut))) < len(mutfound)
					keep = np.add(segregating,ancestral)
					mut = mut[keep]
					pop = pop[:, keep]
			
			# go to next generation
			gen += 1

		# cleanup
		close_output_files(fileHandles)

		# next replicate run
		rep += 1

######################################################################
##RUNNING##
######################################################################    
	
#run (with timer)
start = time.time()
main()
end = time.time()
print(end-start)
