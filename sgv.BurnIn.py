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
	sim_id = 'n%d_K%d_alpha%.1f_B%d_u%.4f_sigma%.3f' %(n, K, alpha, B, u, sigma)
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
	return np.exp(-0.5 * sigma * dist**2) #probability of survival

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

def remove_lost_muts(remove_lost, pop, mut):
	"""
	This function removes lost mutations
	"""
	if remove_lost:
		keep = pop.any(axis=0)
		mut = mut[keep]
		pop = pop[:, keep]
	return [pop, mut]

def remove_fixed_muts(remove_fixed, pop, mut):
	"""
	This function removes fixed mutations
	"""
	if remove_fixed:
		keep = np.concatenate((np.array([True]), np.any(pop[:,1:]-1, axis=0))) #note this is a little more complicated because we want to keep the first, base, mutation despite it being fixed
		mut = mut[keep]
		pop = pop[:, keep]
	return [pop, mut]

def histogram_files(mut, theta, pop, alpha, n, K, B, u, sigma, rep, data_dir):
	"""
	Save csv of mutation sizes (in SGV and de novo) for plotting histograms
	"""
	if make_histogram_files:

		#segregating mutations
		keep = np.any(pop-1, axis=0)

		sgv_dist = np.linalg.norm(mut[keep] - theta, axis=1) #phenotypic distance from optimum for each individual mutation in sgv
		sgv_freq = np.sum(pop[:, keep], axis=0) #number of copies of each mutation
		newmuts = np.random.normal(0, alpha, (len(sgv_dist), n)) #phenotypic effect of new mutations (make same number as in sgv)
		dist_denovo = np.linalg.norm(newmuts - theta, axis=1) #phenotypic distance from optimum for each individual de novo mutation
		
		sim_id = 'n%d_K%d_alpha%.1f_B%d_u%.4f_sigma%.1f_rep%d' %(n, K, alpha, B, u, sigma, rep) #sim info
		filename_sgv = "%s/sgv_muts_%s.csv" %(data_dir, sim_id) #filename for sgv mutations
		filename_denovo = "%s/denovo_muts_%s.csv" %(data_dir, sim_id) #filename for de novo mutations

		#write sgv csv
		with open(filename_sgv, 'w') as csvfile:
			writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
			writer.writerow(sgv_dist)
			writer.writerow(sgv_freq)

	    #write de novo csv
		with open(filename_denovo, 'w') as csvfile:
			writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
			writer.writerow(dist_denovo)

######################################################################
##PARAMETERS##
######################################################################

n = 2 #phenotypic dimensions (positive integer >=1)
K = 10000 #number of individuals (positive integer >=1)
alpha = 0.1 #mutational sd (positive real number)
B = 2 #number of offspring per generation per parent (positive integer)
u = 0.001 #mutation probability per generation per genome (0<=u<=1)
sigma = 0.01 #strength of selection (positive real number)

theta = np.array([0]*n) #optimum phenotype (n real numbers)

maxgen = 50 #total number of generations population adapts for (positive integer)
gen_rec = 10 #print every this many generations (positve integer <=maxgen)

remove_lost = True #If true, remove mutations that are lost
remove_fixed = False #If true, remove mutations that are fixed

make_histogram_files = True #if true ouput mutation sizes for plotting 

reps = 2 #number of replicates (positive integer)

data_dir = 'data/test' #where to save data

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

			# remove lost mutations if doing so
			[pop, mut] = remove_lost_muts(remove_lost, pop, mut)

			if gen % gen_rec == 0:
				
				#print update
				# print(np.mean(pop[:,1:],axis=0))
				print('rep=%d   gen=%d   seg=%d   mfreq=%.3f   mdist=%.3f' %(rep, gen, len(mut), np.mean(np.mean(pop[:,1:],axis=0)), np.mean(np.linalg.norm(mut - theta, axis=1))))

				#save for plotting approach to MS balance (number of segregating mutations and avg frequency)
				seg = np.any(pop-1, axis=0) #segregating mutations only
				mut_seg = mut[seg]
				pop_seg = pop[:, seg]
				write_data_to_output(fileHandles, [rep, gen, len(mut_seg), np.mean(np.mean(pop_seg,axis=0)), np.mean(np.linalg.norm(mut_seg - theta, axis=1))])

			# go to next generation
			gen += 1

		# remove fixed mutations if doing so
		[pop, mut] = remove_fixed_muts(remove_fixed, pop, mut)

		#save mutation and frequency data
		save_arrays(n, K, alpha, B, u, sigma, rep, data_dir, mut, pop)

		#save mutation sizes for plotting histogram
		histogram_files(mut, theta, pop, alpha, n, K, B, u, sigma, rep, data_dir)

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
