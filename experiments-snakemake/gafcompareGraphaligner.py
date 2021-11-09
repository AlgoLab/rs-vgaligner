import argparse
import pandas as pd
import re

# CLI arguments parsing
parser = argparse.ArgumentParser(description='Compare two GAFs given as input.')
parser.add_argument('GAF1', help='Path to the first GAF file')
parser.add_argument('REF', help='Path to the reference GAF file')
args = vars(parser.parse_args())

# Read GAF files
my_gaf = pd.read_csv(args["GAF1"], 
                        sep='\t', 
                        names=["name", "qlen", "qstart", "qend", 
                                "strand", "path", "plen", "pstart", "pend", 
                                "residue", "alblock", "quality", "extra",
                                "extra1", "extra2", "extra3", "extra4"])
ref_gaf = pd.read_csv(args["REF"], 
                        sep='\t', 
                        names=["name", "qlen", "qstart", "qend", 
                                "strand", "path", "plen", "pstart", "pend", 
                                "residue", "alblock", "quality", "extra",
                                "extra1", "extra2"])    # final fields are a bit different, 
                                                        # this should not matter too much
# Compare nodes
#assert len(my_gaf.index) == len(ref_gaf.index), "GAFs have different column numbers"
jaccard_list = []
reads_found = 0
total_ref_reads = len(ref_gaf.index)

for i in range(total_ref_reads):
        ref_gaf_row = ref_gaf.iloc[i]
        ref_gaf_read_name = ref_gaf_row["name"]

        if ref_gaf_read_name in my_gaf["name"].values:

                # The same read has been found in my_gaf
                reads_found += 1
                my_gaf_row = my_gaf[my_gaf["name"] == ref_gaf_read_name].iloc[0]
                
                # Get the string representing the path
                my_gaf_nodes_str = my_gaf_row["path"]
                ref_gaf_nodes_str = ref_gaf_row["path"]

                # Find tuples (orient, nodeid) and perform the comparison 
                my_gaf_tuples = re.findall("(>|<)([0-9]+)", my_gaf_nodes_str)
                ref_gaf_tuples = re.findall("(>|<)([0-9]+)", ref_gaf_nodes_str)

                # Convert ids to integers
                my_gaf_int = list(map(lambda x: +int(x[1]) if x[0]=='>' else -int(x[1]), my_gaf_tuples))
                ref_gaf_int = list(map(lambda x: +int(x[1]) if x[0]=='>' else -int(x[1]), ref_gaf_tuples))

                if my_gaf_int == ref_gaf_int:
                        jaccard = 1.0
                else:
                        # Find intersection and union
                        my_gaf_min = min(my_gaf_int)
                        my_gaf_max = max(my_gaf_int)

                        ref_gaf_min = min(ref_gaf_int)
                        ref_gaf_max = max(ref_gaf_int)

                        intersec = range(max(my_gaf_min, ref_gaf_min), min(my_gaf_max, ref_gaf_max))
                        union = range(min(my_gaf_min, ref_gaf_min), max(my_gaf_max, ref_gaf_max))

                        # Compute jaccard
                        jaccard = len(intersec)/len(union) if len(union) else 0

                print("jaccard for {} is: {}".format(ref_gaf_read_name, jaccard))
                jaccard_list.append(jaccard)

print("Matching reads: {}/{}".format(reads_found, total_ref_reads))
print("AVG Jaccard is: {}".format(sum(jaccard_list)/len(jaccard_list) if jaccard_list else 0))

# Print jaccard as comma-separated string
jaccard_list_string = [str(val) for val in jaccard_list]
print("Jaccard list is: \n {}".format(','.join(jaccard_list_string)))
