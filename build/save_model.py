import dihmm_ext as dm
import matplotlib as mpl
import numpy as np
import dihmm_ext

n_bin_states = 30
n_domain_states = 30
domain_size = 20
tolerance = 1e-6
max_iter = 500

if __name__ == "__main__":

	#training_data_file_path = "/Users/stephanostsoucas/Documents/CDiHmm/CDiHmm/Test.txt"
	#annotation_data_file_path = "/Users/stephanostsoucas/Documents/CDiHmm/CDiHmm/Test.txt"
	annotation_data_file_path = "/Users/stephanostsoucas/Documents/CDiHmm/CDiHmm/GM12878_chr17_binary.txt"
	training_data_file_path = "/Users/stephanostsoucas/Documents/CDiHmm/CDiHmm/GM12878_chr10_binary.txt"
	x = dm.run_dihmm(n_bin_states, n_domain_states, domain_size, max_iter, tolerance, training_data_file_path)

	print(dm.save_model(x, "/Users/stephanostsoucas/c/py_dihmm/build/"))

