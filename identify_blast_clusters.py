#!/usr/bin/env python3

"""
usage: .identify_blast_clusters.py processed_blast_output.txt outfile

infile: blast output that you have processed so that only hits that
you consider good enough to warrant collapsing remain.
outfile: table summarizing identified clusters and their representative
e.g.:
representative_phage\tcluster_member1 cluster_member2 etc

Given blast results that have been thresholded somehow (e.g. only hits
covering 95% or more of subject length), collapses any completely
connected clusters to a single representative and represented members.

The cluster representative is chosen by looking for the largest member
of the cluster and then the member with the highest average %id for its
matches, then the first in the list if multiple are equal.

For incompletely connected clusters, identifies the most highly
connected member node and makes that the representative of all connected
nodes. This is repeated for remaining nodes until no connected nodes
remain.

Outputs a table with a representative of the cluster in one column and
the members of the cluster (space-delimited) in the other.
"""

import sys

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile, 'r') as blastfile:
	blastfilelines = blastfile.readlines()
blastfile.close()

dict_list = []
for line in blastfilelines:
	x_index = -1
	y_index = -1
	x = line.split()[0]
	y = line.split()[1]
	xlen = int(line.split()[4])
	ylen = int(line.split()[3])
	xyid = float(line.split()[2])

	if len(dict_list) == 0:
		dict_list.append({x:[1,xlen,xyid],y:[1,ylen,xyid]})
	else:
		for i in range(len(dict_list)):
			
			if x in dict_list[i].keys():
				x_index = int(i)
			
			if y in dict_list[i].keys():
				y_index = int(i)
			
			if x_index != -1 and y_index != -1:
				if x_index == y_index:
					dict_list[x_index][x][0] += 1
					dict_list[x_index][y][0] += 1
					dict_list[x_index][x][2] += xyid
					dict_list[x_index][y][2] += xyid
				
				else:
					for k,v in dict_list[y_index].items():
						dict_list[x_index][k] = v
					
					dict_list[x_index][x][0] += 1
					dict_list[x_index][y][0] += 1
					dict_list[x_index][x][2] += xyid
					dict_list[x_index][y][2] += xyid
					del dict_list[y_index]
				break
			
			if i == len(dict_list)-1:
				if x_index == -1 and y_index == -1:
					dict_list.append({x:[1,xlen,xyid],y:[1,ylen,xyid]})
				
				elif x_index == -1 and y_index != -1:
					dict_list[y_index][x] = [1,xlen,xyid]
					dict_list[y_index][y][0] += 1
					dict_list[y_index][y][2] += xyid
				
				elif x_index != -1 and y_index == -1:
					dict_list[x_index][x][0] += 1
					dict_list[x_index][x][2] += xyid
					dict_list[x_index][y] = [1,ylen,xyid]
				
				else:
					print("somethings gone wrong looking for entries in the dicts.")

# Find fully connected groups in network

connected_dict_list = []
partially_connected_dict_list = []

for i in dict_list:
	exp_cons = 2*(len(i.keys())-1) # A fully connected group will have 2(n-1) connections
	keep = True
	for k,v in i.items():
		if v[0] != exp_cons:
			keep = False
	if keep:
		connected_dict_list.append(i)
	else:
		partially_connected_dict_list.append(i)


# Pick representative of groups to keep (prefer biggest, then highest average identity with partners, then first in list)

reps_list = []

for i in connected_dict_list:
	best = 0
	longest = 0
	most_id = 0
	for k,v in i.items():
		if v[1] > longest:
			best = k
			longest = v[1]
			most_id = v[2]
		if v[1] == longest:
			if v[2] > most_id:
				best = k
				longest = v[1]
				most_id = v[2]
	to_collapse = [x for x in i.keys() if x != best]

	reps_list.append({best:to_collapse})


# If you want to collapse mostly-connected groups as well. This picks the best genome that is above the threshold with every other in the cluster, but when some members of the cluster are lower than the threshold with one another.

for i in partially_connected_dict_list:
	exp_cons = 2*(len(i.keys())-1) # A fully connected group will have 2(n-1) connections
	best = 0
	longest = 0
	most_id = 0
	if not any([v[0] == exp_cons for v in i.values()]): # If no 1 genome is connected to all the others then process further
		blast_subset = []
		group_of_interest = list(i.keys())
		
		for line in blastfilelines:
			if any([x in line for x in group_of_interest]):
				blast_subset.append(line)
		
		subset_con_dict = {k:[] for k in i.keys()}
		
		for line in blast_subset:
			a, b = line.split()[:2]
			subset_con_dict[a].append(b)
			subset_con_dict[b].append(a)
			subset_con_dict[a] = list(set(subset_con_dict[a]))
			subset_con_dict[b] = list(set(subset_con_dict[b]))

		finished_genomes = []
		rep = 1
		while any([len(v) > 0 for v in subset_con_dict.values()]):
			max_cons = max([len(v) for v in subset_con_dict.values()])
			for k,v in subset_con_dict.items():
				if len(v) == max_cons:
					reps_list.append({k:list(set(v))})
					finished_genomes += [k] + v
					subset_con_dict = {k:[] for k in i.keys()}
					fin_gen_set = set(finished_genomes)
					for line in blast_subset:
						a, b = line.split()[:2]
						if a not in fin_gen_set and b not in fin_gen_set:
							subset_con_dict[a].append(b)
							subset_con_dict[b].append(a)
							subset_con_dict[a] = list(set(subset_con_dict[a]))
							subset_con_dict[b] = list(set(subset_con_dict[b]))
					break

			
	else:
		for k,v in i.items():
			if v[0] == exp_cons:
				if v[1] > longest:
					best = k
					longest = v[1]
					most_id = v[2]
				if v[1] == longest:
					if v[2] > most_id:
						best = k
						longest = v[1]
						most_id = v[2]
		to_collapse = [x for x in i.keys() if x != best]

		reps_list.append({best:to_collapse})


with open(outfile, 'w+') as outf:
	for i in reps_list:
		for k, v in i.items():
			outf.write(k + '\t' + " ".join(v) + '\n')
outf.close()


