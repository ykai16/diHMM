import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import math

class BedWriter(object):

	def __init__(self, annotations, model):

		self.n_bin_states = model.n_bin_states
		self.n_domain_states = model.n_domain_states

		jet = cm = plt.get_cmap('jet') 
		cNorm  = colors.Normalize(vmin=0, vmax=self.n_bin_states)
		scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
		self.color_map = scalarMap
		self.model = model
		self.bin_annotations = self._create_bin_state_ranges(annotations.annotations)
		self.domain_annotations = self._create_domain_state_ranges(annotations.annotations)

	def _create_bin_state_ranges(self, annotations):
		return self._create_annotation_ranges(annotations, 0)

	def _create_domain_state_ranges(self, annotations):
		return self._create_annotation_ranges(annotations, 1)

	# Compress the chromosome annotations into a list of 3-tuples: (state, start index position, end position (exclusive))
	# bin_or_domain is 0 for bins and 1 for domains
	def _create_annotation_ranges(self, annotations, bin_or_domain):

		t_idx = bin_or_domain

		ranges = []

		current_state = annotations[0][t_idx]

		#start of current range, inclusive
		current_range_start = 0

		#end of current range, exclusive
		current_range_end = 0

		for (idx, state) in enumerate(annotations):
			
			if (state[t_idx] == current_state):

				current_range_end += 1
			else:

				ranges.append((current_state, current_range_start, current_range_end))

				current_range_start = current_range_end
				current_range_end += 1
				current_state = state[t_idx]

		ranges.append((current_state, current_range_start, current_range_end))

		return ranges

	def write_bed_files(self, base_dir, cell_line, chr_name):
		self._write_bin_state_bed_file(base_dir, cell_line, chr_name)
		self._write_domain_state_bed_file(base_dir, cell_line, chr_name)

	def _write_bed_header(self, f, cell_line, file_title):
		f.write('track\tname="{}"\tdescription="{}"\tvisibility=1 itemRgb="on"\n'.format(file_title, file_title))

	def _get_color_string(self, state):

		color = self.color_map.to_rgba(state)

		color_lst = []

		for i in range(3):
			color_lst.append(int(math.ceil(color[i] * 255)))

		return ','.join(str(component) for component in color_lst)

	def _write_bin_state_bed_file(self, base_dir, cell_line, chr_name):

		file_title = chr_name + '_bin_states'

		with open(base_dir + '/' + file_title + '.bed', 'w') as f:

			self._write_bed_header(f, cell_line, file_title)

			for r in self.bin_annotations:
				self._write_bed_row(f, r, chr_name)


	def _write_domain_state_bed_file(self, base_dir, cell_line, chr_name):

		file_title = chr_name + '_domain_states'

		with open(base_dir + '/' + file_title + '.bed', 'w') as f:

			self._write_bed_header(f, cell_line, file_title)

			for r in self.domain_annotations:
				self._write_bed_row(f, r, chr_name)		

	def _write_bed_row(self, f, r, chr_name):
		
		state, start, end = r

		# chromosome name, start, end, bin state, 0, ., start, end, color
		f.write('{0}\t{1}\t{2}\t{3}\t0\t.\t{1}\t{2}\t{4}\n'.format(chr_name, (start * 200), (end * 200), 'N' + str(state), self._get_color_string(state)))


