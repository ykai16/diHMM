import sys
import argparse
import dihmm_ext
import bed_writer

#Constant arguments
domain_size = 20
tolerance = 1e-6
max_iter = 500

description = """A command-line tool to use DiHMM.
This can be used to train, save, and load models, as well as annotate chromosomes."""

def get_dihmm_arg_parser():
	parser = argparse.ArgumentParser(description=description)

	commands = parser.add_mutually_exclusive_group(required=True)
	commands.add_argument('-t', '--train', nargs=2, type=int)
	commands.add_argument('-a', '--annotate', action='store_true')
	parser.add_argument('-i', '--input_files', nargs='+', type=str, required=True)
	parser.add_argument('-o', '--output_dir', nargs=1, type=str, required=True, default='.')
	parser.add_argument('-chr', '--chromosome', nargs=1, type=int)
	parser.add_argument('-cell', '--cell_type', nargs=1, type=int)

	return parser

if __name__ == '__main__':
	
	parser = get_dihmm_arg_parser()

	args = parser.parse_args()

	if (args.train):
		
		n_states = args.train

		n_bin_states = n_states[0]
		n_domain_states = n_states[1]

		input_files = args.input_files

		output_dir = args.output_dir

		model = dihmm_ext.run_dihmm(n_bin_states, n_domain_states, domain_size, max_iter, tolerance, input_files)

		dihmm_ext.save_model(model, output_dir)

	elif (args.annotate):
		
		input_files = args.input_files

		output_dir = args.output_dir

		if (len(input_files) == 1):
			raise ValueError('Only one file supplied. At least two files must be supplied: the directory from which to load the trained model, and the chromosome data files to annotate.')

		model = dihmm_ext.load_model(input_files[0], domain_size)

		a = dihmm_ext.annotate(model, input_files[1:])

		bw = bed_writer.BedWriter(a, model)

		if (args.chromosome and args.cell_type):

			#TODO: get the list of chromosomes. better yet just annotate all chromosomes or read chromosome name from file
			bw.write_bed_files(output_dir, args.cell_type, args.chromosome)
		else:
			raise ValueError('Cell type or chromosome not supplied.')

	else:
		print("NOTHING")