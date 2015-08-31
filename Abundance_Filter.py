#!/usr/bin/env python

# Remove all sequences that occur less then X times 

# Arguments: Abundance_Filter.py -s [sequence file] -m [min copy number] -c

# author: Youri Lammers
# contact: youri.lammers@naturalis.nl / youri.lammers@gmail.com

# import the modules used by the script
import os, argparse, itertools, sys, collections

# Retrieve the commandline arguments
parser = argparse.ArgumentParser(description = 'Split a sequence file based on a list of primers.')

parser.add_argument('-f', '--sequence_file', metavar='Sequence file', dest='sequence', type=str,
			help='The sequence file in either fastq or fasta format.')
parser.add_argument('-m', '--min', metavar='Min copy number', dest='min', type=int,
			help='The minimum number of reads per duplicate group (default = 10)', default=1)
parser.add_argument('-c', '--copy_num', dest='copy_num', action='store_true',
			help=('Preserve the number of sequences per duplicate group (quality per sequence is normalized), '
			'if not set only one sequence per group is kept'))
args = parser.parse_args()

def extract_sequences():

	# open the sequence file submitted by the user, get 
	# the file format and rewind the file
	sequence_file = open(args.sequence)
	file_format = sequence_file.readline()[0]
	sequence_file.seek(0)
	
	# create a iterative index of all the headers
        lines = (x[1] for x in itertools.groupby(sequence_file, key=lambda line: line[0] == file_format))

	# walk through the header and obtain the sequences (and quality score if applicable)
        for headers in lines:
		header = headers.next().strip()
		if file_format == '>': sequence = [''.join(line.strip() for line in lines.next())]
		else:
			temporary_list, sequence, quality = [line.strip() for line in lines.next()], [], []
		
			# get the multi line sequences and break at the sequence - quality
			# seperator symbol (+)
			while len(temporary_list) > 0:
				line = temporary_list.pop(0)
				if line[0] == '+':
					break
				sequence.append(line)
			quality = temporary_list

			# if the length of the sequences differs from the length of the
			# quality scores (because the quality line starts with a '@') get
			# the next quality line and append it to the previous one
			while len(quality) < len(sequence):
				if len(quality) == 0 and len(sequence) == 1:
					quality.append(headers.next().strip())
				else:
					quality += [line.strip() for line in lines.next()]
			
			# join the sequence lines and quality lines together
			sequence = [''.join(sequence), [ord(var) for var in ''.join(quality)]]
				
		# yield the header + sequence
		yield [header, sequence]


def write_read(count, read):

	# write the read to the output_file in either fasta or fastq
	# format depending on the read
	if len(read) > 1:
		print '@{0}\n{1}\n+\n{2}'.format(count, read[0], read[1])
	else:
		print '>{0}\n{1}'.format(count, '\n'.join([read[0][i:i+60] for i in range(0, len(read[0]), 60)]))


def main ():

	# create empty nested dictionary for the reads + quality
	# scores if present
	sequence_dic = collections.defaultdict(list)

	# get a read from the input file
	for read in extract_sequences():
		# if the read is in fastq format:
		if len(read[1]) > 1:
			sequence_dic[read[1][0]].append(read[1][1])

		# if the read is in fasta format:
		else:
			sequence_dic[read[1][0]].append('')
	
	count = 1
	# walk through the list of sequences
	for sequence in sequence_dic:
		# get the abundance for each sequence group, if the
		# abundance >= the filter cutoff normalize the quality 
		# if needed and write the sequences
		abundance = len(sequence_dic[sequence])
		if abundance >= args.min:

			# if the sequences are in fasta format
			if sequence_dic[sequence][0] == '':
				# if the copy number is preserved				
				if args.copy_num == True:
					# set the sequence names
					for temp in ['{0}_{1}'.format(count, i) for i in range(1,(abundance+1))]:
						# write the sequences
						write_read(temp, [sequence])
				# if the number is not preserved, just write one sequence
				else:
					write_read(count, [sequence])

			# if the sequences are in fastq format
			else:
				# create an empty string for the consensus quality score
				consensus_qual = ''
				# walk through the sequence base by base
				for pos in range(0, len(sequence)):
					# for each base, get the average value for all sequences in the group
					consensus_qual += chr(sum([q[pos] for q in sequence_dic[sequence]])/abundance)

				# if the copy number is preserved
				if args.copy_num == True:
					# set the sequence names
					for temp in ['{0}_{1}'.format(count,i) for i in range(1,(abundance+1))]:
						# write the sequences
						write_read(temp, [sequence, consensus_qual])
				# if the number is not preserved, just write one sequence
				else:
					write_read(count, [sequence, consensus_qual])
		count += 1
	

if __name__ == '__main__':
	main()
