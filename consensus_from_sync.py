import sys
from random import shuffle
from argparse import ArgumentParser

# command line arguments
parser = ArgumentParser(prog="consensus_from_sync", description="Generate a consensus FASTA sequence from pool-seq data in SYNC format")
parser.add_argument("-f", "--file", required=True,
					help="Input SYNC file to be processed", metavar="sync_file")

args = parser.parse_args()


###################
#### FUNCTIONS ####
###################

# # split all count columns and record counts as integers
# order in SYNC files is A, T, C, G, N, del
# here, I ignore Ns
def get_sync_counts(sync_site_counts):
	count_list = []
	for pop_count in sync_site_counts:
		pop_count = pop_count.split(":")
		pop_count = list(map(int, pop_count))
		n_count = pop_count[4]
		non_n_count = pop_count[0:4] + pop_count[5:]
		count_list.append(non_n_count)

	return(count_list)


##############
#### MAIN ####
##############

# read number of pops from SYNC file
with open(args.file) as f:
	line = f.readline()
	# strip trailing newline
	line = line.strip("\n")
	line = line.split("\t")
	# all colon delimited count columns
	site_counts = line[3:]
	n_pops = len(site_counts)

pop_names = [i for i in range(1, n_pops + 1)]

pop_seq_dict = {}
for pop in pop_names:
	pop_seq_dict[pop] = []

with open(args.file) as f:
	for line in f:
		# strip trailing newline
		line = line.strip("\n")
		line = line.split("\t")
		# name of scaffold / chromosome
		scaff = line[0]
		# genomic position
		pos = line[1]
		# base call in reference
		ref = line[2]
		# all colon delimited count columns
		site_counts = line[3:]
		# split all count columns and record counts as integers
		# order in SYNC files is A, T, C, G, N, del
		# here, I ignore Ns
		count_list = get_sync_counts(site_counts)

		sync_key = ["A", "T", "C", "G", "-"]

		for i, count in enumerate(count_list):
			m = max(count)
			which_max = [i for i, j in enumerate(count) if j == m]
			if len(which_max) > 1:
				# where multiple sites have equal depth, shuffle and pick the first
				shuffle(which_max)
			consensus = sync_key[which_max[0]]
			pop_seq_dict[pop_names[i]].append(consensus)

for pop in pop_names:
	print(">" + str(pop))
	print("".join(pop_seq_dict[pop]))
