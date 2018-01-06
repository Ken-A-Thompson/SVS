#Author: Matthew Osmond <mmosmond@zoology.ubc.ca> & Ken A. Thompson <ken.thompson@zoology.ubc.ca>
#Description: Adaptation from standing genetic variance (SGV) in Fisher's geometric model, implications for hybrids
#Artificially create standing genetic variance in ancestor

import numpy as np
import time
import pickle
import csv

######################################################################
##HELPER FUNCTIONS##
######################################################################

# def open_output_files(K, n, alpha, n_muts, rep):
def open_output_files(K, n, alpha, n_muts, p_mut, rep):
    """
    This function opens the output files and returns file
    handles to each.
    """

    sim_id = 'K%d_n%d_alpha%r_nmuts%r_pmut%r_create_rep%d' %(K, n, alpha, n_muts, p_mut, rep)
    # sim_id = 'K%d_n%d_alpha%r_nmuts%r_create_rep%d' %(K, n, alpha, n_muts, rep)

    data_dir = 'data'

    outfile_A = open("%s/pop_%s.pkl" %(data_dir,sim_id),"wb")
    outfile_B = open("%s/mut_%s.pkl" %(data_dir,sim_id),"wb")

    return [outfile_A, outfile_B]

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

K = 1000 #number of individuals (positive integer)
n = 2 #phenotypic dimensions (positive integer)
alpha = 0.02 #mutational sd (positive real number)
n_muts = 30 #number of mutations (positive integer)
p_mut = 0.1 #probability of having mutation at any one locus (0<p<1)

#meta-parameters
nReps = 1 #number of replicates to run (postitive integer)
run = True #run simulation?
make_CSV = True #make a CSV of output?

######################################################################
##SIMULATION##
######################################################################

def main():

	rep = 1
	while rep < nReps + 1:

		# open output files
		fileHandles = open_output_files(K, n, alpha, n_muts, p_mut, rep) 
		# fileHandles = open_output_files(K, n, alpha, n_muts, rep) 

		
		#create ancestor chromosomes and mutation effects
		pop = np.random.binomial(1, p_mut, (K, n_muts)) #p_mut chance of having each of n_muts mutations, for all K individuals
		# pop = np.random.binomial(1, np.random.uniform(0,1), (K, n_muts)) #p_mut chance of having each of n_muts mutations, for all K individuals
		mut = np.random.normal(0, alpha, (n_muts, n)) #create n_muts mutations, each with a random normal phenotypic effect in each n dimension with mean 0 and sd alpha

		#write and close
		write_data_to_output(fileHandles, [pop,mut])
		close_output_files(fileHandles)

		# run next replicate
		rep += 1

######################################################################
##RUNNING##
######################################################################    

if run == True:    
	#run (with timer)
	start = time.time()
	main()
	end = time.time()
	print(end-start)

######################################################################
##MAKE CSV OF MUTATIONS
######################################################################   

if make_CSV == True: 

	rep=1
	sim_id = 'K%d_n%d_alpha%r_nmuts%r_pmut%r_create_rep%d' %(K, n, alpha, n_muts, p_mut, rep)
	# sim_id = 'K%d_n%d_alpha%r_nmuts%r_create_rep%d' %(K, n, alpha, n_muts, rep)

	data_dir = 'data'

	# load pop data
	f = open('%s/pop_%s.pkl' %(data_dir,sim_id), 'rb')
	popall = []
	while 1:
	    try:
	        popall.append(pickle.load(f))
	    except EOFError:
	        break

	# make csv of population list (mutations in each individual)
	with open("%s/pop_%s.csv" %(data_dir,sim_id), "w") as f:
	    writer = csv.writer(f)
	    # writer.writerows(popall) #write for all timepoints
	    writer.writerows(popall[-1]) #write just last timepoint