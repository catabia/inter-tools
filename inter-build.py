#inter-build.py

"""
# -----------------------------------------------------------------------
#
# inter-build.py
#
# by Hannah Catabia
# Last Modified: 2017-05-30
#
# This script creates a single, species-specific, customized 
# interactome dataset from an unlimited number of MITAB files.
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
#   -m	One or several mitab files
#		USAGE EXAMPLE:
#		> python inter-map.py -m mitab1.txt mitab2.txt
# 
# * Optional input:  
# 
#   -t 	An NCBI taxonomy ID in integer form.
#		If none is given, 9606 (homo sapiens) will be used as default.
#		USAGE EXAMPLE:
#		> python inter-map.py -m mitab1.txt -t 10090
# 
#   -i 	A PSI-MI interaction detection method ancestor in integer form.
#		If none is given, 0045 (experimental detection method) will be used as default.
#		USAGE EXAMPLE:
#		> python inter-map.py -m mitab1.txt -d 0045
#
#   -f 	The format of the output file, either a CSV or TSV.
#		If none is given, CSV will be used as a default.
#		USAGE EXAMPLE:
#		> python inter-map.py -m mitab1.txt -f tsv
#
#   -o 	The name of the output file.
#		If none is given, 'interactome' will be used as a default.
#		USAGE EXAMPLE:
#		> python inter-map.py -m mitab1.txt -o human_interactome
#
#   -d 	Calculate the diameter of the interactome.
#		The default is to not calculate it, as it is time-intensive
#		USAGE EXAMPLE:
#		> python inter-map.py -m mitab1.txt -s human_interactome
# 
#
# -----------------------------------------------------------------------
"""

#various imports
import sys
import argparse
import os.path
import csv
import math
import datetime
import re
import pandas as pd
import networkx as nx
import collections
import numpy as np
#mygene
import mygene
mg = mygene.MyGeneInfo()
#rpy2
import rpy2
from rpy2.rinterface import RRuntimeError
from rpy2.robjects.packages import importr
from rpy2.robjects import StrVector, DataFrame, r
#matplotlib
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles


# =============================================================================
# Double-checks the MITAB file to make sure the columns are in line
def mitab_format_check(header, filename):
	#change to lowercase
	header = [x.lower() for x in header]
	error_message = '\n The MITAB file %s is not properly formatted and will not be included in the interactome.' %filename
	#check each block
	if 'interactor a' not in header[0]:
		print error_message
		return False
	elif 'interactor b' not in header[1]:
		print error_message
		return False
	elif 'detection' not in header[2]:
		print error_message
	elif 'publication' not in header[3]:
		print error_message
		return False
	elif 'taxid' not in header[4] or 'taxid' not in header[5]:
		print error_message
		return False
	elif 'type' not in header[6]:
		print error_message
	elif 'source' not in header[7]:
		print error_message
		return False
	else:
		return True


# =============================================================================
# Finds either an entrez of uniprot ID in the interactor column of the MITABs
def find_label(labels, filename):
	if 'entrez' in labels:
		label = 'entrez'
	elif 'uniprot' in labels:
		label = 'uniprot'
	else:
		print '''\n The MITAB file %s does not contain entrezid or uniprotkb, 
		or is formatted incorrectly. 
		\n Please check file and try running the program again. \n''' %filename
		sys.exit(0)
	return label


# =============================================================================
# Takes the entrez or uniprot ID out of a string
def extract_id(id_list, label):
	id_list = id_list.split('|')
	for item in id_list:
		if label in item:
			item = item.partition(':')[2].strip()
			#get rid of unwanted characters
			item = re.split(',.|;: ', item)[0]
			return item
	return None


# =============================================================================
# Convert MyGene output into a dictionary of dictionaries
def to_dict_of_dict(dict_list):
		#create a dictionary of dictionaries
		dict_of_dict = {}
		for item in dict_list:
			#make sure all gene objects have 'notfound' attribute
			if 'notfound' not in item.keys():
				item['notfound'] = False
			dict_of_dict[item['query']] = item
		return dict_of_dict


# =============================================================================
# Adds Uniprot KB and gene symbol columns to interactions identified by Entrez IDs
def entrez_id_convert(entrez_mitabs, taxids):
	if entrez_mitabs:
		#concatenate the entrez files and make a list of A and B identifier columns
		entrez_concat = pd.concat(entrez_mitabs)
		entrez_concat = entrez_concat.dropna(subset=['idA', 'idB'])
		#create a list of ids to look up with MyGene
		entrez_ids = list(set(entrez_concat['idA'].tolist() + entrez_concat['idB'].tolist()))
		#get translations from MyGene API
		print '\nQuerying MyGene API for %i Entrez IDs...\n' %len(entrez_ids)
		entrez_id_list = mg.querymany(entrez_ids, scope='entrezgene', fields='uniprot,symbol', species=taxids, entrezonly=True)
		print '\n\n\n'
		#create a dictionary of dictionaries from MyGene query
		entrez_dict = to_dict_of_dict(entrez_id_list)
		#mitab as list of lists
		listed = entrez_concat.values.tolist()
		#remove the 'notfound' entries
		listed = [row for row in listed if not entrez_dict[row[0]]['notfound'] and not entrez_dict[row[1]]['notfound']]
		#add columns for uniprot and symbol
		new_listed = []
		for row in listed:
			new_row = []
			#entrez IDs
			new_row.extend(row[0:2])
			#uniprot A
			try:
				new_row.append(str(entrez_dict[row[0]]['uniprot']['Swiss-Prot']))
			except:
				new_row.append(None)
			#uniprot B
			try:
				new_row.append(str(entrez_dict[row[1]]['uniprot']['Swiss-Prot']))
			except:
				new_row.append(None)
			#symbol A
			try:
				new_row.append(str(entrez_dict[row[0]]['symbol']))
			except:
				new_row.append(None)
			#symbol B
			try:
				new_row.append(str(entrez_dict[row[1]]['symbol']))
			except:
				new_row.append(None)
			#rest of row
			new_row.extend(row[2:])
			new_listed.append(new_row)
		return new_listed
	else:
		return []


# =============================================================================
# Adds Entrezt IDs and gene symbol columns to interactions identified by UniProt KBs
def uniprot_id_convert(uniprot_mitabs, taxids):
	#create a dictionary of uniprot ids
	if uniprot_mitabs:
		#concatenate the uniprot files and make a list of A and B identifier columns
		uniprot_concat = pd.concat(uniprot_mitabs)
		uniprot_concat = uniprot_concat.dropna(subset=['idA', 'idB'])
		#create a list of uniprot ids to look up in MyGene
		uniprot_ids = list(set(uniprot_concat['idA'].tolist() + uniprot_concat['idB'].tolist()))
		#get translations from MyGene API
		print '\nQuerying MyGene API for %i UniProt IDs...\n' %len(uniprot_ids)
		uniprot_id_list = mg.querymany(uniprot_ids, scope='uniprot', species=taxids, fields='entrez,symbol', entrezonly=True)
		print '\n\n\n'
		#create a dictionary of dictionaries from MyGene query
		uniprot_dict = to_dict_of_dict(uniprot_id_list)
		#mitab as list of lists
		listed = uniprot_concat.values.tolist()
		#remove the 'notfound' entries
		listed = [row for row in listed if not uniprot_dict[row[0]]['notfound'] and not uniprot_dict[row[1]]['notfound'] ]
		#add columns for entrez and symbol
		new_listed = []
		for row in listed:
			new_row = []
			#entrez A
			try:
				new_row.append(str(uniprot_dict[row[0]]['_id']))
			except:
				pass
			#uniprot B
			try:
				new_row.append(str(uniprot_dict[row[1]]['_id']))
			except:
				pass
			#uniprot IDs
			new_row.extend(row[0:2])
			#symbol A
			try:
				new_row.append(str(uniprot_dict[row[0]]['symbol']))
			except:
				new_row.append(None)
			#symbol B
			try:
				new_row.append(str(uniprot_dict[row[1]]['symbol']))
			except:
				new_row.append(None)
			#rest of row
			new_row.extend(row[2:])
			new_listed.append(new_row)
		return new_listed
	else:
		return []


# =============================================================================
# Use dictionary structure to eliminate duplicate interactions
def to_inter_dict(interactions):
	inter_dict = {}
	for row in interactions:
		try:
			#dictionary key is entrez id, min-to-max
			key = str(min( int(row[0]), int(row[1]) )) + '-' + str(max(int(row[0]), int(row[1])))
		except:
			pass
		if key not in inter_dict:
			inter_dict[key] = row
		else:
			if row[12] not in inter_dict[key][12]:
				#if the dictionary key is already there, add on the file name to the from_files column
				inter_dict[key][12] = inter_dict[key][12] + '|' + row[12]
	#turn the dictionary back into a list of lists, and return it
	interactions = []	
	for key in inter_dict.keys():
		interactions.append(inter_dict[key])
	return interactions

	 # columns:
	 # 0 entrezA, 1 entrezB, 2 uniprotA, 3 uniprotB, 4 symbolA, 5 symbolB, 
	 # 6 interaction_detection_methods, 7 publication_identifiers, 8 taxidA , 9 taxidB, 
	 # 10 interaction_type, 11 source_database, 12 from_files


# =============================================================================
# Extract the approprite information from the MITABs and create a single list of lists
def unify_mitabs(mitab_files, taxids):
	print '\nFinding Entrez and UniProt IDs within the MITAB files...\n'
	#lists of mitab files and dicts of entrez and uniprot id conversions
	entrez_mitabs = []
	uniprot_mitabs = []
	# open and create dataframes from the mitab files
	for mitab_file in mitab_files:
		mitab = []
		with open(mitab_file, 'rb') as tsvfile:
			reader = csv.reader(tsvfile, delimiter = '\t')
			for row in reader:
				mitab.append([row[0], row[1], row[6], row[8], row[9], row[10], row[11], row[12]])
		headers = mitab.pop(0)
		mitab = pd.DataFrame(mitab, columns=headers)
		#checking MITAB is correctly formated:
		if mitab_format_check(list(mitab), mitab_file):
			#rename the columns for query simplicity
			mitab.columns = ['idA', 'idB', 'detection_method', 'publication', 'taxidA', 'taxidB', 'interaction_type', 'source_database']
			mitab['from_files'] = mitab_file
			#figure out if the mitab uses entrez or uniprot for main identifiers
			label = find_label(mitab['idA'][0], mitab_file)
			mitab['idA'] = mitab['idA'].apply(extract_id, args = (label,))
			mitab['idB'] = mitab['idB'].apply(extract_id, args = (label,))
			#place the mitab into the correct list based on whether it is entrez or uniprot
			if label == 'entrez':
				entrez_mitabs.append(mitab)
			elif label == 'uniprot':
				uniprot_mitabs.append(mitab)
	#create a list of lists for
	interactions = entrez_id_convert(entrez_mitabs, taxids) + uniprot_id_convert(uniprot_mitabs, taxids)
	#create a dictionary of non-repeating interactions
	print '\nSearching for duplicate interactions...\n'
	return to_inter_dict(interactions)


# =============================================================================
# Checks interaction detection method ontology (using ontoCAT in R)
def check_detection_methods(interactions, detection):
# keeps only interactions with PSI-MI interaction dection method in ancestry
# detection should be a 4-digit number in string format
	new_data = []
	#To given percentages
	length = len(interactions)
	#dictionary for interactions that have already been inspected
	is_approved = {}
	#Use ontoCAT to check interaction ontology
	for row in interactions:
		#detection method number
		if 'I:' in row[6]:
			method = row[6].split('I:')[1][0:4]
		else:
			method = '0'
		#detection method is the term
		if method == detection:
			new_data.append(row)
		# check the dictionary
		elif method in is_approved:
			if is_approved[method]:
				new_data.append(row)
			else:
				pass
		#check ontology
		else:
			parents = ontoCAT.getAllTermParentsById(MI, "MI_" + method)
			# Boolean for stating whether ancestor appears in ontology
			appears = False
			for term in parents:
				if detection in str(term):
					new_data.append(row)
					appears = True
					break
			#place into approved dictionary
			is_approved[method] = appears
	return new_data


# =============================================================================
# Creates a PDF of a Venn diagram of the 3 largest MITABs as represented in the final dataset
def draw_venn(interactions, mitab_files, output):
	mitab_sets = []
	mitab_dict = {}
	#count the number of interactions from each file
	for mitab in mitab_files:
		mitab_set = []
		for i in range(len(interactions)):
			if mitab in interactions[i][12]:
				mitab_set.append(i)
		mitab_sets.append(set(mitab_set))
		mitab_dict[mitab] = mitab_set
	if len(mitab_files) == 1:
		pass
	elif len(mitab_files) == 2:
		v = venn2(mitab_sets, mitab_files)
	elif len(mitab_files) == 3:
		v = venn3(mitab_sets, mitab_files)
	else:
		lengths = [len(s) for s in mitab_sets]
		zipped = sorted(zip(lengths, mitab_sets, mitab_files), reverse=True)[:3]
		mitab_sets3 = [zipped[0][1], zipped[1][1], zipped[2][1]]
		mitab_files3 = [zipped[0][2], zipped[1][2], zipped[2][2]]
		v = venn3(mitab_sets3, mitab_files3)
	plt.savefig('venn_' + output + '.pdf', bbox_inches='tight')
	return mitab_dict


# =============================================================================
# Create a CSV or TSV file with the final interactome dataset
def create_interactome_file(interactions, output_file, file_format):
	header = ['entrez_A', 'entrez_B', 'uniprot_A', 'uniprot_B', 'symbol_A', 'symbol_B', 'interaction_detections_methods',
				'publication_identifiers', 'taxid_A', 'taxid_B', 'interaction_type', 'source_database', 'from_files']
	interactions.insert(0, header)
	delimiter = {'tsv': '\t', 'csv': ','}[file_format]
	wr = csv.writer(open(output_file + '.' + file_format, 'wb'), delimiter=delimiter)
	wr.writerows(interactions)

# =============================================================================
# Writes a text file with information about the final interactome dataset
# Also creates the histogram file
def create_text_file(interactions, taxids, detection, mitab_dict, output, diameter):
	interactions.pop(0)
	#info will be appended with all of the text
	info =  '''
%s
----------------------------------------------------------------------------
%s
MITAB Files: %s

	''' % (output, '{:%Y-%b-%d_%H:%M:%S}'.format(datetime.datetime.now()), ', '.join(mitab_files))

	info = info + '''	
INPUTS:
Taxonomy IDs: %s
Interaction detection methods: %s

	'''%(taxids, detection)

	# put interactome into a graph
	G = nx.Graph()
	for row in interactions:
		G.add_edge(row[0], row[1])
	#numbers of genes and interactions
	nodes = len(G.nodes())
	edges = len(G.edges())
	#count self-edges
	self_edges = 0
	for edge in G.edges():
		if edge[0] == edge[1]:
			self_edges += 1
	info = info + '''	
FINAL INTERACTOME:
This interactome contains:
%i unique genes
%i unique interactions, of which %i (%f %%) are interactions

	'''%(nodes, edges, self_edges, float(self_edges)/float(edges)*100 )

	for mitab in mitab_dict.keys():
		this_set = set(mitab_dict[mitab])
		other_sets = [value for key, value in mitab_dict.items() if key != mitab]
		other_sets = set([item for sublist in other_sets for item in sublist])
		difference = this_set.difference(other_sets)
		info = info + '''
MITAB FILE: %s
%i (%f %%) of the %i interactions in the final interactome came from %s.
%i  (%f %%) of the %i interactions in the final interactome may ONLY be found in %s.

		'''%(mitab, len(this_set), float(len(this_set))/float(edges), edges, mitab,
			len(difference), float(len(difference))/float(edges), edges, mitab)

	#get the largest connected component (lcc)
	connected_components = sorted([list(x) for x in nx.connected_components(G)], key=len, reverse = True)
	LCC = G.subgraph(connected_components[0])
	OCC = G.subgraph([item for sublist in connected_components[1:] for item in sublist])
	info = info + '''
LARGEST CONNECTED COMPONENT
The final interactome network contains %i components.
The LCC of the final interactome contains:
%i genes
%i interactions
Inter-Tools will use the LCC of the final interactome for its analysis.
All together, the other connected components contain:
%i genes
%i interactions
Inter-Tools will ignore these in its analysis.

	'''%(len(connected_components), len(LCC.nodes()), len(LCC.edges()), len(OCC.nodes()), len(OCC.edges()))

	#opt-in for diameter
	if diameter:
		print '\n\n Calculating the diameter of the interactome...\n\n'
		interactome_diameter = nx.diameter(LCC)
		info = info + '''
DIAMETER
The interactome has a diameter of %i.

		''' %(interactome_diameter,)

	#create this histogram
	#source: https://networkx.github.io/documentation/development/examples/drawing/degree_histogram.html
	degree_sequence = sorted(nx.degree(G).values(),reverse=True) # degree sequence
	degreeCount=collections.Counter(degree_sequence)
	deg, cnt = zip(*degreeCount.items())
	print cnt
	#log base 10 scale
	cnt = np.log10(cnt)
	print cnt
	fig, ax = plt.subplots()
	plt.bar(deg, cnt, width=0.80, color='b')
	plt.title("Histogram of Interactions per Gene")
	plt.ylabel("Number of Genes (Log 10 Scale)")
	plt.xlabel("Interactions per Gene (Degree)")
	ax.set_xticks([d+0.4 for d in deg])
	ax.set_xticklabels(deg)
	plt.savefig("histogram_" + output + ".pdf")
	#plt.show()
	#percentiles
	perc100 = np.percentile(degree_sequence, 100)
	perc75 = np.percentile(degree_sequence, 75)
	perc50 = np.percentile(degree_sequence, 50)
	perc25 = np.percentile(degree_sequence, 25)
	perc1 = np.percentile(degree_sequence, 1)	
	info = info + '''
DISTRIBUTION OF INTERACTIONS PER GENE (NODE DEGREE)
100th percentile:  %i
75th percentile:  %i
50th percentile:  %i
25th percentile:  %i
1st percentile:  %i

		''' %(perc100, perc75, perc50, perc25, perc1)

	with open('BUILD-RESULTS_'+ output + '.txt', 'w') as text_file:
		text_file.write(info)
	print info



'''
# =============================================================================
# PROGRAM INITIALIZES HERE
# =============================================================================
'''

if __name__ == '__main__':
	print #blank line to begin

	#-------------------------------------------------------------------------------------------------
	# Parsing the command line
	#-------------------------------------------------------------------------------------------------
	parser = argparse.ArgumentParser()
	parser.add_argument('-m', '--mitab', type=str, nargs='+', dest='mitab_files', action='store', default= None, help='List of MITAB files separated by spaces, ie: -m mitab1.txt mitab2.txt')
	parser.add_argument('-t', '--taxid', type=str, nargs='+', dest='taxids', action='store', default='9606', help='NCBI Taxonomy ID in integer form, ie: -t 9606')
	parser.add_argument('-i', '--detection', type=str, nargs='+', dest='detection', action='store', default='0045', help= 'PSI-MI Interaction Detection Method ancestor in integer form, ie: -d 0045')
	parser.add_argument('-f', '--format', type=str, dest='format', action='store', default='csv', help='Output file format in CSV or TSV, ie: -f csv')
	parser.add_argument('-o', '--outfile', type=str, dest='output', action='store', default='interactome_' + '{:%Y-%b-%d-%H:%M:%S}'.format(datetime.datetime.now()), help='Name for interactome file, ie: -o human_interactome' )
	parser.add_argument('-d', '--diameter', dest='diameter', action='store_true', default=False, help= 'Calculate the diameter of the interactome, ie: -s')
	# Unpacking the command line inputs
	args = parser.parse_args()
	mitab_files = args.mitab_files
	taxids = args.taxids
	detection = args.detection
	format = args.format
	output = args.output
	diameter = args.diameter

	#-------------------------------------------------------------------------------------------------
	# Checking/installing/loading ontoCAT (BioconductoR) and PSI-MI Ontology
	#-------------------------------------------------------------------------------------------------
	# ontoCAT installed here as variable for later use
	try:
		ontoCAT = importr('ontoCAT')
	except RRuntimeError:
		print '''

        This program requires ontoCAT (BioconductoR package) to run.

        Currently installing ontoCat to your machine...

        '''
		base = importr('base')
		base.source("http://www.bioconductor.org/biocLite.R")
		biocinstaller = importr("BiocInstaller")
		biocinstaller.biocLite("ontoCAT")
		ontoCAT = importr('ontoCAT')
	# PSI-MI stored here as variable MI for later use
	print '\nLoading PSI-MI Ontology...\n'
	MI = ontoCAT.getOntology('https://raw.githubusercontent.com/MICommunity/psidev/master/psi/mi/rel25/data/psi-mi25.obo')

	#-------------------------------------------------------------------------------------------------
	# Error checking the command line
	#-------------------------------------------------------------------------------------------------
	#checking for MITAB file input
	mitab_error = '''

	ERROR: You must specify at least one MITAB file to process.
	This file should be in the same directory as the program.

	EXAMPLE:
	./mitab_preprocess.py -m file1.mitab.txt file2.mitab.txt file3.mitab.txt

	'''
	if mitab_files == None:
		print mitab_error
		sys.exit(0)
	for mitab_file in mitab_files:
		if not os.path.exists(mitab_file):
			print mitab_error
			sys.exit(0)

	#checking for error in detection method input
	detection_error = '''
	ERROR: You must enter a 4-digit PSI-MI interaction detection method code.

	PSI-MI ontology can be explored at this site:
	http://www.ebi.ac.uk/ontology-lookup/browse.do?ontName=MI

	'''
	# make sure it's a PSI-MI detection method and get it's name
	try:
		method = ontoCAT.getTermById(MI, 'MI_' + detection)
	except:
		print detection_error
		sys.exit(0)

	#checking for error in file format input
	format_error = '''
	ERROR: You may specify only CSV (comma-separated value) 
	or TSV (tab-delimited) for the output file. For example:
	./mitab_preprocess.py -m file.mitab.txt -f TSV

	Defaults to CSV if no input is made.


	'''
	
	if format.lower() not in ['csv', 'tsv']:
		print format_error
		sys.exit(0)


	#------------------------------------------------
	# Assembling the interactome
	#------------------------------------------------

	#getting the mitabs into a unified list-of-lists
	interactions = unify_mitabs(mitab_files, taxids)
	#checking the interaction deteftion method
	print 'Extracting interactons with the following PSI-MI interaction detection method in its ancestry:'
	print method, '\n'
	interactions = check_detection_methods(interactions, detection)
	mitab_dict = draw_venn(interactions, mitab_files, output)
	#save interactome
	print '\nCreating an interactome database file named %s.%s...\n' %(output, format.lower())
	create_interactome_file(interactions, output, format.lower())
	#save info text file
	print '\nCreating an informational text file named %s.txt...\n' %(output, )
	create_text_file(interactions, taxids, detection, mitab_dict, output, diameter)
	# And we're...
	print '\nDone!\n\n'








