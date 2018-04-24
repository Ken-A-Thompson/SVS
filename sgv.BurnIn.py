#Authors: Matthew Osmond <mmosmond@zoology.ubc.ca> & Ken A. Thompson <ken.thompson@zoology.ubc.ca>
#Description: The role of standing genetic variance in speciation

import numpy as np
import time
import csv
import math

######################################################################
##FUNCTIONS##
######################################################################

def open_output_file(n, K, alpha, B, u, sigma, data_dir):
	"""
	This function opens the output files and returns file
	handles to each.
	"""
	sim_id = 'n%d_K%d_alpha%.1f_B%d_u%.4f_sigma%.1f' %(n, K, alpha, B, u, sigma)
	outfile_A = open("%s/PlotBurn_%s.csv" %(data_dir, sim_id), "w")
	return outfile_A

def write_data_to_output(fileHandles, data):
	"""
	This function writes to the corresponding output file.
	"""
	writer = csv.writer(fileHandles)
	writer.writerow(data)

def close_output_files(fileHandles):
	"""
	This function closes all output files.
	"""
	fileHandles.close()

def save_arrays(n, K, alpha, B, u, sigma, rep, data_dir, mut, pop):
	"""
	Save numpy arrays of mutations and their frequencies.
	"""
	sim_id = 'n%d_K%d_alpha%.1f_B%d_u%.4f_sigma%.1f_rep%d' %(n, K, alpha, B, u, sigma, rep)
	
	filename = "%s/Muts_%s.npy" %(data_dir, sim_id)
	np.save(filename, mut[1:]) #save mutations
	
	filename = "%s/Freqs_%s.npy" %(data_dir, sim_id)
	np.save(filename, np.sum(pop[:,1:],axis=0)/len(pop)) #save frequencies

def survival(sigma, dist):
	"""
	This function gives the probability of survival
	"""
	return np.exp(-sigma * dist**2) #probability of survival

def viability(phenos, theta, pop, K, sigma):
	"""
	This function determines which individuals survive viability selection
	"""
	dist = np.linalg.norm(phenos - theta, axis=1) #phenotypic distance from optimum
	w = survival(sigma, dist) #probability of survival
	rand = np.random.uniform(size = len(pop)) #random uniform number in [0,1] for each individual
	surv = pop[rand < w] #survivors
	if len(surv) > K:
		surv = surv[np.random.randint(len(surv), size = K)] #randomly choose K_adapt individuals if more than K_adapt
	return surv

def recomb(surv, B):
	"""
	This function creates offspring through pairing of parents (haploid) and recombination (i.e, meiosis)
	"""
	pairs = np.resize(np.random.choice(len(surv), size=len(surv), replace=False), (int(len(surv)/2), 2)) #random mate pairs (each mates at most once and not with self)
	rand2 = np.random.randint(2, size=(len(pairs), len(surv[0]))) #from which parent each offspring inherits each allele (free recombination, fair transmission)
	rec = np.resize(np.append(rand2, 1-rand2, axis=1),(len(rand2), 2, len(rand2[0]))) #reshape
	off_1 = np.sum(surv[pairs] * rec, axis=1) #one product of meiosis
	off_2 = np.sum(surv[pairs] * (1-rec), axis=1) #other product of meiosis
	off = np.repeat(np.append(off_1, off_2, axis=0), B, axis=0) #each product of meiosis produced B times
	return off

def mutate(off, u, alpha, n, mut):
	"""
	This function creates mutations and updates population
	"""
	rand3 = np.random.uniform(size = len(off)) #random uniform number in [0,1] for each offspring
	nmuts = sum(rand3 < u) # mutate if random number is below mutation rate; returns number of new mutations
	whomuts = np.where(rand3 < u) #indices of mutants
	newmuts = np.random.normal(0, alpha, (nmuts, n)) #phenotypic effect of new mutations
	pop = np.append(off, np.transpose(np.identity(len(off), dtype=int)[whomuts[0]]), axis=1) #add new loci and identify mutants
	mut = np.append(mut, newmuts, axis=0) #append effect of new mutations to mutation list
	return [pop, mut]

def remove_muts(remove_lost, remove_fixed, pop, mut):
	"""
	This function creates mutations and updates population
	"""
	if remove_lost:
		keep = pop.any(axis=0)
		mut = mut[keep]
		pop = pop[:, keep]
	if remove_fixed:
		keep = np.concatenate((np.array([True]), np.any(pop[:,1:]-1, axis=0))) #note this is a little more complicated because we want to keep the first, base, mutation despite it being fixed
		mut = mut[keep]
		pop = pop[:, keep]
	return [pop, mut]

######################################################################
##PARAMETERS##
######################################################################

n = 2 #phenotypic dimensions (positive integer >=1)
K = 1000 #number of individuals (positive integer >=1)
alpha = 0.1 #mutational sd (positive real number)
B = 2 #number of offspring per generation per parent (positive integer)
u = 0.001 #mutation probability per generation per genome (0<=u<=1)
sigma = 0.05 #strength of selection (positive real number)

theta = np.array([0]*n) #optimum phenotype (n real numbers)

maxgen = 1000 #total number of generations population adapts for (positive integer)
gen_rec = 100 #print every this many generations (positve integer <=maxgen)

remove_lost = True #If true, remove mutations that are lost
remove_fixed = True #If true, remove mutations that are fixed

reps = 10 #number of replicates (positive integer)

data_dir = 'data/burnins' #where to save data

######################################################################
##MAIN SIMULATION##
######################################################################

def main():

	fileHandles = open_output_file(n, K, alpha, B, u, sigma, data_dir)

	rep = 1
	while rep < reps + 1:

		pop = np.array([[1]] * K) #start all individuals with one "mutation"
		mut = np.array([[0] * n]) #but let the "mutation" do nothing (ie put all phenotypes at origin)

		#intitialize generation counter
		gen = 0

		#run until maxgen
		while gen < maxgen + 1:

			# genotype to phenotype
			phenos = np.dot(pop, mut) #sum mutations held by each individual (ie additivity in phenotype)

			# viability selection
			surv = viability(phenos, theta, pop, K, sigma)

			# end simulation if population extinct (or unable to produce offspring)        
			if len(surv) < 2: 
				print("Extinct")              
				break 

			# meiosis
			off = recomb(surv, B)

			# mutation and population update
			[pop, mut] = mutate(off, u, alpha, n, mut)

			# remove lost mutations (all zero columns in pop)
			[pop, mut] = remove_muts(remove_lost, remove_fixed, pop, mut)

			if gen % gen_rec == 0:
				
				#print update
				print('rep=%d   gen=%d   seg=%d' %(rep, gen, len(mut)))

				#save for plotting approach to MS balance (number of mutations and avg frequency)
				write_data_to_output(fileHandles, [rep, gen, len(mut), np.mean(np.sum(pop[:,1:],axis=0)/len(pop))])

			# go to next generation
			gen += 1

		#save mutation and frequency data
		save_arrays(n, K, alpha, B, u, sigma, rep, data_dir, mut, pop)

		# go to next rep
		rep += 1

	close_output_files(fileHandles)

######################################################################
##RUNNING ADAPTATION FUNCTION##
######################################################################    
	
start = time.time()
main()
end = time.time()
print('this took %.2f seconds to complete' %(end-start))
