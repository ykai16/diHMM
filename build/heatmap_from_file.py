import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import dihmm_ext
import time

def plot_state_coverage(n_figures, state_array, n_states, title):
	plt.figure(n_figures)
	plt.title(title)
	#plt.yticks(np.arange(0, n_domain_states, 1.0))
	a = np.expand_dims(state_array, axis=0)
	plt.imshow(a, cmap='Blues', interpolation='nearest', norm=mpl.colors.LogNorm())

def plot_bin_state_transitions():
	pass

def plot_emission_probabilities():
	pass


if __name__ == "__main__":

	n_bin_states = 30
	n_domain_states = 30
	domain_size = 20
	tolerance = 1e-6
	max_iter = 500
	#training_data_file_path = "/Users/stephanostsoucas/Documents/CDiHmm/CDiHmm/Test.txt"
	#annotation_data_file_path = "/Users/stephanostsoucas/Documents/CDiHmm/CDiHmm/Test.txt"
	annotation_data_file_path = "/Users/stephanostsoucas/Documents/CDiHmm/CDiHmm/GM12878_chr17_binary.txt"
	training_data_file_path = "/Users/stephanostsoucas/Documents/CDiHmm/CDiHmm/GM12878_chr10_binary.txt"

	start = time.time()
	x = dihmm_ext.run_dihmm(n_bin_states, n_domain_states, domain_size, max_iter, tolerance, [training_data_file_path])
	end = time.time()

	print('DiHMM training took {} seconds'.format(end - start))
	#x_m = dihmm_ext.run_dihmm(n_bin_states, n_domain_states, domain_size, max_iter, tolerance, training_data_file_path)

	#dihmm_ext.save_model(x_m, "/Users/stephanostsoucas/c/py_dihmm/build/")

	#x = dihmm_ext.load_model("/Users/stephanostsoucas/c/py_dihmm/build/", n_bin_states, n_domain_states, domain_size)

	t_d = x.domain_transition_probabilities

	n_figures = 1

	plt.figure(n_figures)
	n_figures += 1

	#plt.subplot(110)
	#plt.xticks(np.arange(0, n_domain_states, 1.0))
	#plt.yticks(np.arange(0, n_domain_states, 1.0))
	#plt.axis([0, n_domain_states - 1, 0, n_domain_states - 1])
	plt.imshow(t_d, cmap='Blues', interpolation='nearest')

	t_b = x.bin_transition_probabilities

	plt.figure(n_figures)
	n_figures += 1
	plt.imshow(x.emission_probabilities, cmap='Blues', interpolation='nearest')
	#plt.axis([0, n_bin_states - 1, 0, n_bin_states - 1])

	#TODO: draw emissions probabilities

	n_plots = t_b.shape[0]
	
	subplot_key = n_plots * 100 + 10
	
	for i in range(0, n_plots):
		plt.figure(n_figures)
		n_figures += 1
		plt.title('Domain %d' % i)
		plt.yticks(np.arange(0, n_bin_states, 1.0))
		plt.xticks(np.arange(0, n_bin_states, 1.0))
		plt.imshow(t_b[i], cmap='Blues', interpolation='nearest')

	a = dihmm_ext.annotate(x, annotation_data_file_path)
	print(a.domain_state_distributions)
	print(a.bin_state_distributions)
	plot_state_coverage(n_figures, a.domain_state_distributions, n_domain_states, "Domain state coverage")
	n_figures += 1
	plot_state_coverage(n_figures, a.bin_state_distributions, n_bin_states, "Bin state coverage")
	plt.show()

	