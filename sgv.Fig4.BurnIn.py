#Authors: Matthew Osmond <mmosmond@zoology.ubc.ca> & Ken A. Thompson <ken.thompson@zoology.ubc.ca>
#Description: The role of standing genetic variance in speciation

import numpy as np
import time
# import matplotlib.pyplot as plt
import csv
import math

######################################################################
##FUNCTIONS##
######################################################################

def open_output_files(n, K, alpha, B, data_dir):
	"""
	This function opens the output files and returns file
	handles to each.
	"""
	sim_id = 'n%d_K%d_alpha%.1f_B%d' %(n, K, alpha, B)
	outfile_A = open("%s/heatmap_%s.csv" %(data_dir, sim_id), "w")
	return outfile_A

def write_data_to_output(fileHandles, data):
	"""
	This function writes a (time, data) pair to the
	corresponding output file. We write densities
	not abundances.
	"""
	writer = csv.writer(fileHandles)
	writer.writerow(data)

def close_output_files(fileHandles):
	"""
	This function closes all output files.
	"""
	fileHandles.close()

def found(n_muts, nmuts_max, ancestor_muts, ancestor_freqs, K, n):
	"""
	This function creates a founding population from an ancestral one
	"""

	#make ancestor
	if n_muts > 0:
		probs = ancestor_freqs/ancestor_freqs.sum() #probability of choosing each mutation (make pdf)
		mut_choice = np.random.choice(nmuts_max, size=n_muts, replace=False, p=probs) #choose n_muts different mutations
		mutfound = ancestor_muts[mut_choice] #mutational effects
		p_mut = ancestor_freqs[mut_choice] #expected frequency of these mutations
		popfound = np.random.binomial(1, p_mut, (K, n_muts)) #p_mut chance of having each of n_muts mutations, for all K individuals
	else: #de novo only, even if p_mut>0
		popfound = np.array([[1]] * K)
		mutfound = np.array([[0] * n])
	return [popfound, mutfound]

def survival(dist):
	"""
	This function gives the probability of survival
	"""
	return np.exp(-0.5 * dist**2) #probability of survival

######################################################################
##UNIVERSAL PARAMETERS##
######################################################################

n = 2 #phenotypic dimensions (positive integer >=1)
data_dir = 'data'

######################################################################
##PARAMETERS OF ANCESTOR##
######################################################################

n_reps = 10 #number of reps of ancestor that exist
K_ancestor = 10000 #number of individuals (positive integer >=1)
alpha = 0.1 #mutational sd (positive real number)
B = 2 #number of offspring per generation per parent (positive integer)
u = 0.001 #mutation probability per generation per genome (0<u<1)
sigma = 0.1 #selection strength
burn_dir = 'data/burnins'
rrep = np.random.choice(n_reps, 1, replace=False) #randomly assign each rep an ancestor

######################################################################
##PARAMETERS TO MAKE HYBRID##
######################################################################

K = 1000 #number of individuals (positive integer >=1)
n_mut_list = [5 * x for x in range(2)] #number of mutations (positive integer >=1)

######################################################################
##PARAMETERS FOR PARENTAL POPULATIONS##
######################################################################

opt_dist = 0.6 #distance to optima (eg, 1 puts you on the unit circle)
theta1 = np.array([opt_dist, 0]) #set one optima to be fixed

n_angles = 10 #number of angles between optima to simulate (including 0 and 180)
angles = [math.pi * x / (n_angles - 1) for x in range(n_angles)] #angles to use (in radians)
theta2_list = np.array([[opt_dist*math.cos(x), opt_dist*math.sin(x)] for x in angles]) #optima to use

######################################################################
##MAIN FUNCTION##
######################################################################

def main():

	# open output files
	fileHandles = open_output_files(n, K, alpha, B, data_dir) 

	#loop over all n_muts values
	i = 0
	while i < len(n_mut_list):

		#set number of mutations in hybrid
		n_muts = n_mut_list[i] 

		#load ancestor
		burn_id = 'n%d_K%d_alpha%.1f_B%d_u%.4f_sigma%.1f_rep%d' %(n, K, alpha, B, u, sigma, rrep+1)

		filename = "%s/Muts_%s.npy" %(burn_dir, burn_id)
		ancestor_muts = np.load(filename) #load mutations

		filename = "%s/Freqs_%s.npy" %(burn_dir, burn_id)
		ancestor_freqs = np.load(filename) #load frequencies
		nmuts_max = len(ancestor_freqs) #number of mutations in ancestor

		#make "hybrids"
		[pop, mut] = found(n_muts, nmuts_max, ancestor_muts, ancestor_freqs, K, n)

		#calculate hybrid load
		phenos = np.dot(pop, mut) #phenotypes
		meanpheno = np.mean(phenos, axis=0) #mean pheno
		dist = np.linalg.norm(phenos - meanpheno, axis=1) #phenotypic distance from mean hybrid
		hyload = np.log(1*B) - np.mean(np.log(survival(dist)*B)) #hybrid load as defined by Chevin et al 2014

		#loop over optima
		j = 0
		while j < len(theta2_list):
		
			#set optima
			theta2 = theta2_list[j]

			#move mean hybrid pheno to mean of parents
			hymean = np.mean([theta1,theta2], axis=0) #the mean we want
			offphenos = phenos - meanpheno + hymean	#shift all phenos to achieve wanted mean

			dist1 = np.linalg.norm(offphenos - theta1, axis=1) #phenotypic distance from parental 1 optimum
			dist2 = np.linalg.norm(offphenos - theta2, axis=1) #phenotypic distance from parental 2 optimum
			dist = np.minimum(dist1, dist2) #distance to closest optimum
			w = survival(dist) #viability
			fitness = np.log(w*B) #continuous time growth rate
			meanfit = np.mean(fitness) #mean fitness
			maxfit = np.amax(fitness) #max fitness over all hybrids
			percentilefit = np.percentile(fitness, 90) #max fitness over all hybrids

			#print an update
			print('nmuts=%d, hybrid load=%.3f, angle=%r, mean fitness=%.2f, max fitness=%.2f' %(n_muts, hyload, round(angles[j]*180/math.pi,2), meanfit, maxfit)) 

			#save data
			write_data_to_output(fileHandles, [n_muts, hyload, round(angles[j]*180/math.pi,2), meanfit, maxfit, percentilefit])

			#go to next optima
			j += 1

		#go to next n_muts value
		i += 1

	# cleanup
	close_output_files(fileHandles)

######################################################################
##RUNNING ADAPTATION FUNCTION##
######################################################################    
	
start = time.time()
main()
end = time.time()
print('this took %.2f seconds to complete' %(end-start))
