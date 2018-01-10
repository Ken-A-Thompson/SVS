#Authors: Matthew Osmond <mmosmond@zoology.ubc.ca> & Ken A. Thompson <ken.thompson@zoology.ubc.ca>
#Description: The role of standing genetic variance in speciation

import numpy as np
import time
import matplotlib.pyplot as plt
import csv
import math

######################################################################
##FUNCTIONS##
######################################################################

def open_output_files(n, K, p_mut, alpha, B, data_dir):
	"""
	This function opens the output files and returns file
	handles to each.
	"""
	sim_id = 'n%d_K%d_pmut%.1f_alpha%.1f_B%d' %(n, K, p_mut, alpha, B)
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
##PARAMETERS TO MAKE HYBRID##
######################################################################

K = 1000 #number of individuals (positive integer >=1)
n_mut_list = [10 * x for x in range(11)] #number of mutations (positive integer >=1)
p_mut = 0.1 #probability of having mutation at any one locus (0<=p<=1) #set this to zero for de novo only
alpha = 0.1 #mutational sd (positive real number)

B = 2 #to determine growth rates (Malthusian fitness)

######################################################################
##PARAMETERS FOR PARENTAL POPULATIONS##
######################################################################

opt_dist = 1 #distance to optima (eg, 1 puts you on the unit circle)
theta1 = np.array([opt_dist, 0]) #set one optima to be fixed

n_angles = 5 #number of angles between optima to simulate (including 0 and 180)
angles = [math.pi * x / (n_angles - 1) for x in range(n_angles)] #angles to use (in radians)
theta2_list = np.array([[opt_dist*math.cos(x), opt_dist*math.sin(x)] for x in angles]) #optima to use

######################################################################
##MAIN FUNCTION##
######################################################################

def main():

	# open output files
	fileHandles = open_output_files(n, K, p_mut, alpha, B, data_dir) 

	#loop over all n_muts values
	i = 0
	while i < len(n_mut_list):

		#set number of mutations in hybrid
		n_muts = n_mut_list[i] 

		#make hybrid
		if n_muts > 0: #at least 1 mutation
			pop = np.random.binomial(1, p_mut, (K, n_muts)) #p_mut chance of having each of n_muts mutations, for all K individuals
			mut = np.random.normal(0, alpha, (n_muts, n)) #create n_muts mutations, each with a random normal phenotypic effect in each n dimension with mean 0 and sd alpha
		else: #no mutations
			pop = np.array([[1]] * K)
			mut = np.array([[0] * n])

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

			#print an update
			print('nmuts=%d, hybrid load=%.3f, angle=%r, mean fitness=%.2f, max fitness=%.2f' %(n_muts, hyload, round(angles[j]*180/math.pi,2), meanfit, maxfit)) 

			#save data
			write_data_to_output(fileHandles, [n_muts, hyload, round(angles[j]*180/math.pi,2), meanfit, maxfit])

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
