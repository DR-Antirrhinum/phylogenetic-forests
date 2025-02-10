import itertools
from argparse import ArgumentParser
import treeXY_funcs as tf


##########################
# command line arguments #
##########################
parser = ArgumentParser(prog="treeXY",
                        description="Calculate window-averaged pi-statistics from mapped pool-seq data")

parser.add_argument("-f", "--file",
                    required=True,
                    help="Input SYNC file to be processed.",
                    metavar="sync_file")

parser.add_argument("-m", "--min_depth",
                    type=int,
                    default=15,
                    help="Minimum depth, across all populations, for a site to be included")

parser.add_argument("-M", "--max_depth",
                    type=int,
                    default=200,
                    help="Maximum depth, across all populations, for a site to be included")

parser.add_argument("-a", "--min_allele_depth",
                    type=int,
                    default=2,
                    help="Minimum depth for an allele call")

parser.add_argument("-A", "--min_allele_pops",
                    type=int,
                    default=2,
                    help="Minimum number of populations required for an allele call")

parser.add_argument("--ignore_multiallelic",
                    action="store_true",
                    help="Ignore sites with >2 alleles, rather than removing least common allele")

args = parser.parse_args()

##########
# treeXY #
##########

# generate file names from args
sync_name = args.file
sync_name = sync_name.split(".")[0]
sync_name = sync_name.split("/")[-1]

h_args = args.min_depth, args.max_depth, args.min_allele_depth, args.min_allele_pops
h_args = list(map(str, h_args))
w_file_name = sync_name + "_" + "_".join(h_args) + "_treeXY.csv"

pop_names = tf.read_pop_names(args.file)

# open window output file
with open(w_file_name, "w") as out_file:
    piw_headers = ["piw_" + name for name in pop_names]
    pit_headers = ["_".join(pair) for pair in list(itertools.combinations(pop_names, 2))]
    pit_headers = ["piT_" + name for name in pit_headers]
    dxy_headers = ["_".join(pair) for pair in list(itertools.combinations(pop_names, 2))]
    dxy_headers = ["dXY_" + name for name in dxy_headers]
    D_headers = ["_".join(pair) for pair in list(itertools.combinations(pop_names, 2))]
    D_headers = ["D_" + name for name in D_headers]
    # write header
    out_file.write("scaffold" + "," + "window_start" + "," + "window_end" + "," + "n_window_sites" + "," +
                   ",".join(piw_headers) + "," + ",".join(pit_headers) + "," + ",".join(dxy_headers) + "," +
                   ",".join(D_headers) + "\n")

# messy way to get max_pos for output
window_details = tf.initialise_windows(args.file, 50000, 25000)
max_pos = window_details[1]

# initialise lists to store summed stats for averaging
pos_piw_sums = [0] * len(pop_names)
pop_comps = [i for i in itertools.combinations(pop_names, 2)]
pos_pit_sums = [0] * len(pop_comps)
pos_dxy_sums = [0] * len(pop_comps)
pos_D_sums = [0] * len(pop_comps)

n_valid = 0
with open(args.file) as file:
    for line in file:
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
        count_list = tf.get_sync_counts(site_counts)
        pop_names = [str(i) for i in range(1, len(count_list) + 1)]

        # read depth checks
        pop_dpth = tf.check_read_depth(count_list, args.min_allele_depth, args.min_depth, args.max_depth)
        if sum(pop_dpth) == len(pop_names):
            # get all possible alleles
            allele_stats = tf.check_allele_num(count_list, pop_dpth, args.min_allele_depth)
            pop_dpth = allele_stats[0]
            alleles = allele_stats[1]
            valid_allele_num = True

            if len(alleles) > 2 and args.ignore_multiallelic:
                valid_allele_num = False
            else:
                while len(alleles) > 2:
                    line = tf.filter_triallelic(count_list, scaff, pos, ref)
                    # strip trailing newline
                    line = line.strip("\n")
                    line = line.split("\t")
                    # all colon delimited count columns
                    site_counts = line[3:]
                    count_list = tf.get_sync_counts(site_counts)
                    # get alleles
                    allele_stats = tf.check_allele_num(count_list, pop_dpth, args.min_allele_depth)
                    pop_dpth = allele_stats[0]
                    alleles = allele_stats[1]

            # proceed with piT / dXY / tree calculations for biallelic and monoallelic sites
            if valid_allele_num and len(alleles) > 0 and all([i > 0 for i in pop_dpth]):
                n_valid += 1
                # get piw, piT, and dXY
                pos_piw_vals = tf.get_site_stats(alleles, count_list, pop_names, pop_dpth)[0]
                pos_pit_vals = tf.get_site_stats(alleles, count_list, pop_names, pop_dpth)[1]
                pos_dxy_vals = tf.get_site_stats(alleles, count_list, pop_names, pop_dpth)[2]
                pos_D_vals = tf.get_site_stats(alleles, count_list, pop_names, pop_dpth)[3]

                pos_piw_sums = [x + y for x, y in zip(pos_piw_sums, pos_piw_vals)]
                pos_pit_sums = [x + y for x, y in zip(pos_pit_sums, pos_pit_vals)]
                pos_dxy_sums = [x + y for x, y in zip(pos_dxy_sums, pos_dxy_vals)]
                pos_D_sums = [x + y for x, y in zip(pos_D_sums, pos_D_vals)]

piw_average = [str(i / n_valid) for i in pos_piw_sums]
pit_average = [str(i / n_valid) for i in pos_pit_sums]
dxy_average = [str(i / n_valid) for i in pos_dxy_sums]
da_average = [str(i / n_valid) for i in pos_D_sums]
key = range(1, max_pos + 1)
w_average = [str(n_valid), piw_average, pit_average, dxy_average, da_average]

tf.write_window(scaff, key, w_average, w_file_name)
