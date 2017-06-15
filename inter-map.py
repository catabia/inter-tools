#inter-map.py

"""
# -----------------------------------------------------------------------
#
# inter-map.py
#
# by Hannah Catabia
# Last Modified: 2017-05-30
#
# This script creates GEXF files showing a user-input gene set
# mapped onto an interactome network.
# 
# 
# Inter-Tools: A Command-line Toolkit for Interactome Research
#
# by Hannah Catabia, Caren Smith, and Jose Ordovas
# 
# 
# -----------------------------------------------------------------------
# 
# 
# This program will produce a CSV or TSV file containing
# protein-protein interactions curated from user-input
# MITAB files.
# 
# 
# * Required input:
# 
#	-g  A CSV or TSV file with a list of genes in the first
#		column, identified with Entrez IDs or Uniprot KBs
# 	
#	-e 	Indicates the the genes are identified with Entrez IDs
#	
#	-u 	Indicates the the genes are identified with UniProt KBs
#
#		ENTREZ ID USAGE EXAMPLE:
#		> python inter-map.py -g gene_set.csv -e
#
#		UNIPROT KB USAGE EXAMPLE:
#		> python inter-map.py -g gene_set.csv -u
# 
# * Optional input:  
# 
#	-i 	A CSV or TSV file containing an interactome dataset with
#		the following column headers:
#			entrez_A, entrez_B, uniprot_A, uniprot_B, symbol_A, symbol_B,
#			interaction_detection_methods, publication_identifiers,
#			taxid_A, taxid_B, interaction_type, source_database, from_files
#		If none is given, the script will search for 'interactome.csv'
#		in your directory.
#		USAGE EXAMPLE:
#		> python inter-map.py -g gene_set.csv -e -i mouse_interactome.tsv
#
#	-s 	Include self-interactions in the interactome network.
#		USAGE EXAMPLE:
#		> python inter-map.py -g gene_set.csv -u -s
# 
#
# -----------------------------------------------------------------------
"""

import sys
import argparse
import os.path
import csv
import pandas as pd
import networkx as nx
import datetime
from matplotlib import pyplot as plt
import random
import numpy as np

'''
# =============================================================================
           S T A R T   D E F I N I T I O N S 
# =============================================================================
'''

# =============================================================================
# Error checks all user input from the command line
def command_line_error_check(args):

	#Check to make sure the gene_file and interactome_files are a CSV or a TSV
	file_error= '''
	\n\n
	ERROR: Both the interactome and the gene set must be either a CSV or TSV file.
		%s or %s does not end with the file extension .csv or .tsv.

	EXAMPLE:
		./inter-map.py -g gene_set.tsv -i interactome.csv -e
	\n\n
	'''%(args.gene_file, args.interactome_file)
	
	if args.gene_file[-4:].lower() not in ['.csv', '.tsv']:
		print file_error
		sys.exit(0)
	if args.interactome_file[-4:].lower() not in ['.csv', '.tsv']:
		print file_error
		sys.exit(0)

	#Make sure the gene_file and interactome_file are in the same directory
	path_error = '''
	\n\n
	ERROR: Both the interactome and the gene set files must be in
		the same directory as this program.
	\n\n
	'''
	if not os.path.isfile(args.gene_file) or not os.path.isfile(args.interactome_file):
		print path_error
		sys.exit(0)

	#Make sure we have an Entrez or UniProt identifier designation
	id_error = '''
	\n\n
	ERROR: You must designate the gene set as a list of either Entrez IDs
		or UniProt IDs.  If no command is given, the program will assume
		you have provided Entrez IDs.

	EXAMPLE FOR ENTREZ IDs:
		./inter-map.py -g gene_set.tsv -i interactome.csv -e
	
	EXAMPLE FOR UNIPROT KBs:
		./inter-map.py -g gene_set.tsv -i interactome.csv -u
	\n\n
	'''
	if args.entrez==args.uniprot:
		print id_error
		sys.exit(0)


# =============================================================================
# Open the gene set file and read in the genes
def read_gene_set(gene_file):
	#Is it a CSV or a TSV?
	if gene_file[-4:].lower() == '.tsv':
		delim = '\t'
	else:
		delim = ','
	#open the file and put the genes into a list
	with open(gene_file, 'rb') as gene_file:
		gene_reader = csv.reader(gene_file, delimiter = delim)
		genes = [str(row[0]) for row in gene_reader]
	#check for a header that is not a number (entrez_id)
	for item in ['entrez', 'gene', 'id', 'symbol']:
		if item in genes[0].lower():
			genes.pop(0)
	return genes

# =============================================================================
# Open and read the interactome dataset file
def read_interactome(interactome_file):
	#Is it a CSV or a TSV?
	if interactome_file[-4:].lower() == '.tsv':
		delim = '\t'
	else:
		delim = ','
	#read it using pandas
	interactome = pd.read_csv(interactome_file, sep = delim, usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12])
	#enforce that the interactome file have this specific header
	header = ['entrez_A', 'entrez_B', 'uniprot_A', 'uniprot_B', 'symbol_A', 'symbol_B', 'interaction_detections_methods',
				'publication_identifiers', 'taxid_A', 'taxid_B', 'interaction_type', 'source_database', 'from_files']
	header_error = '''
	\n
	%s must have the following 13 column headers replicated EXACTLY:

	entrez_A, entrez_B, uniprot_A, uniprot_B, symbol_A, symbol_B,
	interaction_detection_methods, publication_identifiers,
	taxid_A, taxid_B, interaction_type, source_database, from_files
	\n
	''' %(interactome_file,)
	if list(interactome) != header:
		print header_error
		sys.exit(0)
	return interactome

# =============================================================================
# Make an graph in NetworkX interactome dataset
def graph_interactome(interactome_df, entrez, uniprot, self):
	G = nx.Graph()
	#get rid of missing id entries
	if entrez:
		interactome_df.dropna(subset=['entrez_A', 'entrez_B'])
	else:
		interactome_df.dropna(subset=['uniprot_A', 'uniprot_B'])
	#make a list
	interactome_list = interactome_df.values.tolist()
	#list the self-edges
	self_edges = []
	for row in interactome_list:
		if entrez:
			node1 = str(row[0])
			node2 = str(row[1])
		else:
			node1 = str(row[2])
			node2 = str(row[3])
		#adding interactions while dealing with self-edges
		if node1==node2 and not self:
			self_edges.append(node1)
		else:
			G.add_node(node1, entrez=row[0], uniprot=row[2], symbol=row[4], taxid=row[8], gene_set=False, weight = 0, interactome_degree = 0, gene_set_degree = 0)
			G.add_node(node2, entrez=row[1], uniprot=row[3], symbol=row[5], taxid=row[9], gene_set=False, weight = 0, interactome_degree = 0, gene_set_degree = 0)
			G.add_edge(node1, node2, interaction_detection_method=row[6], publication_identifiers=row[7], interaction_type=row[10], source_database=row[11], from_files=row[12])
			if node1==node2:
				self_edges.append(node1)
	#make the interactome into one connected component
	components = sorted([list(x) for x in nx.connected_components(G)], key=len)
	largest_component = components.pop()
	G = nx.subgraph(G, largest_component)
	#make the interactome degree an attribute of each gene
	for gene in G.nodes():
		G.node[gene]['interactome_degree'] = nx.degree(G, gene)
	# a list of eliminated genes
	eliminated = [item for sublist in components for item in sublist]
	return G, self_edges, eliminated

# =============================================================================
# Map the gene set into the interactome network
def map_gene_set(G, genes):
	#lists of found and not found genes
	found = []
	not_found = []
	for gene in genes:
		if gene in G:
			G.node[gene]['gene_set'] = True
			found.append(gene)
		else: 
			not_found.append(gene)
	return G, found, not_found


# =============================================================================
# Create GEXF and PDF files
def create_visual_files(G, found, entrez, uniprot):
	#the gene set subgraph
	J = G.subgraph(found)
	#make an attribute for the gene set degree
	#weight = gene_set_degree/interactome_degree
	for gene in found:
		G.node[gene]['gene_set_degree'] = nx.degree(J, gene)
		G.node[gene]['weight'] = (float(G.node[gene]['gene_set_degree'])/float(G.node[gene]['interactome_degree']))*100
	#connected and unconnected subsets of gene set subgraph
	connected_nodes = []
	unconnected_nodes = []
	for node in found:
		if J.degree(node) > 0:
			connected_nodes.append(node)
		else:
			unconnected_nodes.append(node)
	#create a PDF picture of the connected nodes
	H = J.subgraph(connected_nodes)
	#GEXF files
	print '\nWriting the GEXF files...\n'
	nx.write_gexf(J, 'gene_set_network.gexf')
	nx.write_gexf(G, 'full_interactome_network.gexf')
	#Done!
	return H, connected_nodes, unconnected_nodes

# =============================================================================
# Find the largest connected component of the gene set in the interactome network
def get_lcc(G, H):
	#get the connected components
	connected_components = sorted([list(x) for x in nx.connected_components(H)], key=len)
	lcc = connected_components.pop()
	lcc_graph = H.subgraph(lcc)
	return lcc

# =============================================================================
# Find the mean shortest distance of the gene set in the interactome network
def get_msd(G, H, connected_nodes, unconnected_nodes):
	print '\nCalculating the mean shortest distance (MSD)...\n'
	#connected nodes have distance of 1
	shortest_distances = [1] * len(connected_nodes)
	#find all distances between unconnected node pairs:
	pair_lengths = dict((node, {}) for node in unconnected_nodes)
	for nodeA in unconnected_nodes:
		for nodeB in unconnected_nodes:
			if nodeA != nodeB:
				if nodeB not in pair_lengths[nodeA]:
					path_len= nx.shortest_path_length(G, source=nodeA, target=nodeB)
					pair_lengths[nodeA][nodeB] = path_len
					pair_lengths[nodeB][nodeA] = path_len
	#now add in distances for connected nodes
	for nodeA in unconnected_nodes:
		for nodeB in connected_nodes:
			pair_lengths[nodeA][nodeB] = nx.shortest_path_length(G, source=nodeA, target=nodeB)
	#now add shortest distances
	for node in unconnected_nodes:
		shortest_distances.append(min(pair_lengths[node].values()))
	#take the mean
	msd = float(sum(shortest_distances))/len(shortest_distances)
	return msd

'''
# =============================================================================
# PROGRAM INITIALIZES HERE
# =============================================================================
'''

if __name__ == '__main__':
	#A few blank lines to begin with
	print '\n\n\n'

	#-------------------------------------------------------------------------------------------------
	# Parsing and error checking the command line
	#-------------------------------------------------------------------------------------------------
	parser = argparse.ArgumentParser(description='Map a gene set to the interactome')
	parser.add_argument('-g', '--genes', type=str, dest='gene_file', action='store', default= '', help = 'List of genes by Entrez ID or UniProt KB (single column in a CSV or TSV file).')
	parser.add_argument('-i', '--interactome', type=str, dest='interactome_file', action='store', default= 'interactome.csv', help='An interactome dataset created with inter-build.py')
	parser.add_argument('-e', '--entrez', dest='entrez', action='store_true', default=False, help='Use if the gene set uses Entrez IDs')
	parser.add_argument('-u', '--uniprot', dest='uniprot', action='store_true', default=False, help='Use if the gene set uses UniProt KBs')
	parser.add_argument('-s', '--self', dest='self', action='store_true', default=False, help= 'Use if you wish to include self-interactions in the interactome')
	args = parser.parse_args()
	# Error-checking the command line
	command_line_error_check(args)
	#upack the variables
	gene_file = args.gene_file
	interactome_file = args.interactome_file
	entrez = args.entrez
	uniprot = args.uniprot
	self = args.self

	
	#-------------------------------------------------------------------------------------------------
	# Create network and add to the results file
	#-------------------------------------------------------------------------------------------------
	results =  '''
MAP-RESULTS.txt
----------------------------------------------------------------------------
%s
Gene set: %s
Interactome dataset: %s

	''' % ('{:%Y-%b-%d_%H:%M:%S}'.format(datetime.datetime.now()), gene_file, interactome_file)

	#read in genes and interactome interactome
	genes = read_gene_set(gene_file)
	interactome = read_interactome(interactome_file)
	#graph the interactome
	G, self_edges, eliminated = graph_interactome(interactome, entrez, uniprot, self)
	#update results
	results = results + '''
INTERACTOME INFO:
%s contains %i unique genes and %i unique interactions.
	''' % (interactome_file, len(G.nodes()), len(G.edges()) )

	if self:
		results = results + '''
%i of interactions in the dataset are self-interactions. \n
		''' % len(self_edges)
	else:
		results = results + '''
%i self-interactions were eliminated from the dataset. \n
		''' % len(self_edges)

	results = results + '''
The following %i genes were eliminated from the interactome because they 
were not connected to the largest connected component of the interactome:

%s

	'''%(len(eliminated), ', '.join(eliminated))


	#map the gene set to the itneractome
	G, found, not_found = map_gene_set(G, genes)
	#update the results file
	total = len(found) + len(not_found)

	results = results + '''
GENE SET INFO:
%s contains %i unique genes.
%i genes (%f %%) were found in the interactome.
The folowing %i genes were not found in the interactome dataset:

%s
	''' % (interactome_file, total, len(found), float(len(found))/float(total)*100, len(not_found), ', '.join(not_found) )

	#create the visual files and the connected subgraph
	H, connected_nodes, unconnected_nodes = create_visual_files(G, found, entrez, uniprot)

	#update results
	results = results + '''
Of the %i genes found, %i (%f %%) were connected to at least one other gene in the set.
The following %i (%f %%) genes were found in the interactome, but were not directly connected another gene from the gene set:

%s

	''' % (len(found), len(connected_nodes), float(len(connected_nodes))/float(total)*100, 
		len(unconnected_nodes), float(len(unconnected_nodes))/float(total)*100, ', '.join(unconnected_nodes) )

	#find the largest connected component of the gene set
	lcc = get_lcc(G, H)
	#update results
	results = results + '''
LARGEST CONNECTED COMPONENT:
The LCC of the gene set contains %i genes:
%s

	''' % (len(lcc), ', '.join(lcc))

	#find the mean shortest distance of the gene set
	msd = get_msd(G, H, connected_nodes, unconnected_nodes)
	#update results
	results = results + '''
MEAN SHORTEST DISTANCE (MSD):
The MSD of the gene set is %f.

	''' % (msd,)

	#Write the results file
	with open('MAP-RESULTS.txt', 'w') as results_file:
		results_file.write(results)
	print results

	#Finished
	print '\n\nDone!\n\n'