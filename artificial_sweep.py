import sys
import random

in_col = sys.argv[1]

for i in open(in_col).readlines():
	i = i.strip("\n")
	counts = i.split(":")
	counts = [int(count) for count in counts]
	# ignore N and del
	n_count = counts[4]
	del_count = counts[5]
	counts = counts[0:4]
	if counts.count(0) < 3:
		tot_count = 0
		n_list = []
		for i, base in enumerate(counts):
			if base > 0:
				tot_count += base
				n_list.append(str(i)*base)
		n_list = "".join(n_list)
		base_to_fix = int(random.choice(n_list))
		new_counts_list = [0,0,0,0,n_count,del_count]
		new_counts_list[base_to_fix] = tot_count
		new_counts_list = [str(count) for count in new_counts_list]
		print(":".join(new_counts_list))
	else:
		print(i)


