#Author: Matthew Osmond <mmosmond@zoology.ubc.ca> & Ken A. Thompson <ken.thompson@zoology.ubc.ca>
#Description: Adaptation from standing genetic variance (SGV) in Fisher's geometric model, implications for hybrids
#Adapt using standing genetic variance from burn-in as well as new mutations

import numpy as np
import time
import pickle

######################################################################
##HELPER FUNCTIONS##
######################################################################

def open_output_files(K_adapt, n_adapt, B_adapt, u_adapt, alpha_adapt, sigma_adapt, theta_adapt, sigma_opt_adapt, maxgen_adapt, nfounders, opt1, style, rep):
	"""
	This function opens the output files and returns file
	handles to each.
	"""

	sim_id = 'K%d_n%d_B%d_u%r_alpha%r_sigma%r_theta%r_sigmaopt%r_gens%d_founders%d_opt%s_adapt_%s_rep%d' %(K_adapt, n_adapt, B_adapt, u_adapt, alpha_adapt, sigma_adapt, theta_adapt, sigma_opt_adapt, maxgen_adapt, nfounders,'-'.join(str(e) for e in opt1), style, rep)
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
##WHICH ANCESTOR TO START FROM##
######################################################################

K = 1000 #max number of parents (positive integer)
n = 2 #number of traits (positive integer)
B = 2 #number of offspring per generation per parent (positive integer)
u = 0.01 #mutation probability per genome (0<u<1)
alpha = 0.02 #mutational sd (positive real number)
sigma = 1 #strength of selection (positive real number)
theta = 1 #drift parameter in Ornstein-Uhlenbeck movement of the phenotypic optimum (positive real number; setting to zero makes Brownian motion)
sigma_opt = 0.1 #diffusion parameter in Ornstein-Uhlenbeck movement of the phenotypic optimum (positive real number; setting to zero makes constant optimum at opt0)

#meta-parameters
maxgen = 1000 #number of gens in ancestral burn-in (positive integer)
rep = 1 #which replicate

#file to load
sim_id = 'K%d_n%d_B%d_u%r_alpha%r_sigma%r_theta_%r_sigmaopt%r_gens%r_burn_rep%d' %(K,n,B,u,alpha,sigma,theta,sigma_opt,maxgen,rep)
data_dir = 'data'

######################################################################
##COMMON PARAMETERS FOR ADAPTING POPULATIONS##
######################################################################

K_adapt = K #max number of parents (positive integer)
n_adapt = n #number of traits (positive integer)
B_adapt = B #number of offspring per generation per parent (positive integer)
u_adapt = u #mutation probability per genome (0<u<1) (set to zero for sgv only)
alpha_adapt = alpha #mutational sd (positive real number)
sigma_adapt = sigma #strength of selection (positive real number)
theta_adapt = 0 #drift parameter in Ornstein-Uhlenbeck movement of the phenotypic optimum (positive real number; setting to zero makes Brownian motion)
sigma_opt_adapt = 0 #diffusion parameter in Ornstein-Uhlenbeck movement of the phenotypic optimum (positive real number; setting to zero makes constant optimum at opt0)
# opt1 = np.array([0.1] * n_adapt) #optimum phenotype 
opt1 = np.array([-0.1] * n_adapt)

#meta-parameters
maxgen_adapt = 1000 #maximum number of generations (positive integer)
outputFreq = 100 #record and print update this many generations
nReps = 1 #number of replicates

remove_lost = True #If true, remove mutations that are lost (0 for all individuals)
# remove_lost = False
remove = 'derived' #.. any derived (not from ancestor) mutation that is lost 

######################################################################
##SIMULATION GENETICS##
######################################################################

# style = 'both' #standing genetic variation and de novo mutation
style = 'dnm' #de novo mutation only

######################################################################
##SGV (and DNM if u_adapt > 0)##
######################################################################

if style == 'both':

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
	nfounders = min(K_adapt,len(popall[-1])) #size of founding population
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

	nfounders = K_adapt #initial population size
	popfound = np.array([[1]] * nfounders) #list of mutations held by each individual (all start with same mutation at first locus and no others)
	mutfound = np.array([[0] * n_adapt]) #list of phenotypic effect of initial mutations (mutation at first locus puts all individuals at origin)

######################################################################
##SIMULATION##
######################################################################

def main():

	rep = 1
	while rep < nReps + 1:

		#intialize
		pop = popfound
		mut = mutfound
		opt = opt1 #optimum

		# open output files
		fileHandles = open_output_files(K_adapt, n_adapt, B_adapt, u_adapt, alpha_adapt, sigma_adapt, theta_adapt, sigma_opt_adapt, maxgen_adapt, nfounders, opt1, style, rep) 

		gen = 0 #generation
		while gen < maxgen_adapt + 1:

			# genotype to phenotype
			phenos = np.dot(pop,mut) #sum mutations held by each individual

			#optimum
			opt = opt + theta_adapt * (opt1 - opt) + sigma_opt_adapt * np.random.normal(size=n_adapt) # allow optimum to change by Ornstein-Uhlenbeck process, always pulled back toward 0

			# viability selection
			dist = np.linalg.norm(phenos - opt, axis=1) #phenotypic distance from optimum
			w = np.exp(-sigma_adapt*dist**2) #probability of survival
			rand1 = np.random.uniform(size = len(pop)) #random uniform number in [0,1] for each individual
			surv = pop[rand1 < w] #survivors
			if len(surv) > K_adapt:
				surv = surv[np.random.randint(len(surv), size = K_adapt)] #randomly choose K_adapt individuals if more than K_adapt
					
			#end simulation if extinct        
			if len(surv) == 0: 
				print("Extinct")              
				break 
				
			#otherwise continue
			# dump data every outputFreq iteration
			# also print a short progess message (generation and number of parents)
			if gen % outputFreq == 0:
				write_data_to_output(fileHandles, [surv,mut,gen])
				print("rep %d    gen %d    N %d    n_muts %d    mean_muts %.2f    mean_pheno [%.2f, %.2f]" %(rep, gen, len(surv), len(mut), np.mean(np.sum(pop,axis=1)), np.mean(phenos,axis=0)[0], np.mean(phenos,axis=0)[1]))
			
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
			nmuts = sum(rand3 < u_adapt) # mutate if random number is below mutation rate; returns number of new mutations
			whomuts = np.where(rand3 < u) #indices of mutants
			newmuts = np.random.normal(0, alpha, (nmuts,n_adapt)) #phenotypic effect of new mutations

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
