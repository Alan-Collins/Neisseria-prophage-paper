"""
usage: python3 ./blast_average_link_hier_clust_output_clusters.py blast_file

purpose: to process and all vs all blast of MGEs to produce a distance
matrix (consisting of 1 - percent length aligned between MGE pairs).
scipy is then used to perform hierarchical agglomerative clustering of 
MGEs using the average linkage method. the resulting dendrogram's branch
color information is used to output a list of cluster members for each 
cluster.

blast_file is all vs all MGEs blastn with -outfmt "6 std slen qlen"

Outputs:

pla.txt: the % length aligned of phage1 against phage2
    phage1\tphage2\tpla
clusters.txt: Clusters produced by scipy. One cluster per line
    cluster_number\tlist of cluster members
dendrogram.(png|svg): plot of dendrogram in png and svg format
tree.nwk: tree in newick format for analysis in software such as iToL

"""
import sys
import numpy
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree, fcluster, set_link_color_palette
from matplotlib import pyplot as plt
from matplotlib.colors import rgb2hex
from collections import defaultdict
import seaborn as sns

def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    

    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick


blast_file = sys.argv[1]


MGE_size = {}
aligned_regions = defaultdict(list)
# {qid: {sid: [list, of, aligned, regions]}}

# For each MGE, make a list of the regions of its genome that align with
# each other MGE
with open(blast_file, 'r') as fin:
# with open("test_blast.txt", 'r') as fin:
    for line in fin:
        bits = line.split()
        qid = bits[0]
        sid = bits[1]
        qstart = int(bits[6])
        qend = int(bits[7])
        qlen = int(bits[13])

        if qstart > qend:
            qstart, qend = qend, qstart

        # If not done yet, store MGE size
        if qid not in MGE_size:
            MGE_size[qid] = qlen

        aligned_regions[(qid,sid)].append([qstart,qend])
      
# Process aligned_regions to identify total length that aligned
# Divide that length by length of phage and store in new dict

pla_dict = {}
# {qid: [(sid, pla), ...]}

for (qid, sid), loc_ls in aligned_regions.items():
    loc_ls.sort(key=lambda x: x[0]) # Work through list in order
    uniq_regions = [loc_ls[0]]
    for loc in loc_ls:
        # Check if locs are overlapping existing region
        extend_right_check = [(
        # Starts within existing region
        loc[0] >= existing[0] and loc[0] <= existing[1]
        # End outside existing region
        and loc[1] > existing[1]
        ) for existing in uniq_regions]
        if any(extend_right_check):
            idx = extend_right_check.index(True)
            uniq_regions[idx][1] = loc[1]
            continue

        # Else check if this start is outside all existing regions
        if loc[0] > uniq_regions[-1][1]:
            uniq_regions.append(loc)

    # Sum uniq_region lengths
    aligned_len = sum([i[1]-i[0]+1 for i in uniq_regions])
    pla_dict[(qid,sid)] = aligned_len/MGE_size[qid]

with open("pla.txt", 'w') as fout:
    for (qid, sid), pla in pla_dict.items():
        fout.write("{}\t{}\t{}\n".format(qid, sid, pla))

dist_dict = {}
# Just keep the highest pla score for each pair as clustering wants
# score to be the same in both directions.
for (qid, sid), pla in pla_dict.items():
    if (sid, qid) in pla_dict.keys():
        if pla > pla_dict[(sid, qid)]:
            dist_dict[(sid, qid)] = 1-pla
    else:
      dist_dict[(qid, sid)] = 1-pla

phage_tuples = list(dist_dict.keys())

phages = []

for i in phage_tuples:
    phages.append(i[0])
    phages.append(i[1])

unique_phages = sorted(list(set(phages)))

dist_list = []

n_samples = len(unique_phages)

for i in range(n_samples):
    for j in range(i + 1, n_samples):
        if (unique_phages[i], unique_phages[j]) in dist_dict.keys():
            dist_list.append(dist_dict[(unique_phages[i], unique_phages[j])])
        elif (unique_phages[j], unique_phages[i]) in dist_dict.keys():
            dist_list.append(dist_dict[(unique_phages[j], unique_phages[i])])
        else:
            dist_list.append(1)

numpyarray = numpy.array(dist_list)

# to make small dend:
#   dn = dendrogram(Z, leaf_rotation = 90, labels=unique_phages)
#   comment out: fig.set_size_inches(20, 100) 

# to make big dend: 
#   dn = dendrogram(Z, orientation='left', labels=unique_phages)
#   fig.set_size_inches(20, 100)

Z = linkage(numpyarray, 'average') # average linkage clustering method; returns clustering encoded as a linkage matrix Z


sns.set_palette('flare', 160)
palette = sns.color_palette()
hex_palette = []
for i in palette:
  hex_palette.append(rgb2hex(i))
set_link_color_palette(hex_palette)



# fig = plt.figure(figsize=(100, 10))
fig = plt.figure()
dn = dendrogram(Z, orientation='left', labels=unique_phages, get_leaves=True) #  orientation='left', 
# fig.tight_layout()
fig.set_size_inches(20, 100) #100, 10)
plt.savefig("dendrogram.png", dpi=600)
plt.savefig("dendrogram.svg")

cluster_idxs = defaultdict(list)
cluster_leaves = defaultdict(list)

terminal_branch_count = 0

for c, pi, y in zip(dn['color_list'], dn['icoord'], dn['dcoord']):
    for x, leg in enumerate(pi[1:3]):
        if x ==0:
            if y[0] == 0 or y[1] == 0:
                i = (leg - 5.0) / 10.0
                if abs(i - int(i)) < 1e-5:
                    cluster_idxs[c].append(int(i))
                    terminal_branch_count+=1
        elif x == 1:
            if y[2] == 0 or y[3] == 0:
                i = (leg - 5.0) / 10.0
                if abs(i - int(i)) < 1e-5:
                    cluster_idxs[c].append(int(i))
                    terminal_branch_count+=1

idxs_all =[]
for i in cluster_idxs.values():
    idxs_all += i


for c, l in cluster_idxs.items():
    i_l = [dn['ivl'][i] for i in l]
    cluster_leaves[c] = i_l

# print(len(cluster_leaves.keys())) 

cluster_file_contents = []
clus_num = 1
for colour, cluster in cluster_leaves.items():
    if colour != 'C0':
        cluster_file_contents.append("{}\t{}".format(clus_num, " ".join(cluster)))
        clus_num += 1
    else:
        cluster_file_contents.append("Unclustered\t{}".format(" ".join(cluster)))

with open("clusters.txt", 'w') as fout:
    fout.write('\n'.join(cluster_file_contents) + '\n')


tree = to_tree(Z,False)
newick = getNewick(tree, "", tree.dist, unique_phages)

with open("tree.nwk", 'w') as outfile:
    outfile.write(str(newick))
