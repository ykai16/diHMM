import dihmm_ext as dm
import matplotlib as mpl
import numpy as np

n_bin_states = 5
n_domain_states = 5
domain_size = 4
tolerance = 1e-6
max_iter = 500

if __name__ == "__main__":

	#training_data_file_path = "/Users/stephanostsoucas/Documents/CDiHmm/CDiHmm/Test.txt"
	#annotation_data_file_path = "/Users/stephanostsoucas/Documents/CDiHmm/CDiHmm/Test.txt"
	annotation_data_file_path = "/Users/stephanostsoucas/Documents/CDiHmm/CDiHmm/GM12878_chr17_binary.txt"
	training_data_file_path = "/Users/stephanostsoucas/Documents/CDiHmm/CDiHmm/GM12878_chr10_binary.txt"

	x = dm.load_model("/Users/stephanostsoucas/c/py_dihmm/build/", domain_size)

	a_s = dm.annotate(x, [annotation_data_file_path, training_data_file_path])

	print(a_s)

	for a in a_s:
		print(a)