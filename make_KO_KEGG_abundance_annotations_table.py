#!/usr/bin/env python

"""
Author: Jennah Dharamshi
Date: 210217

Make KEGG KO abundance table accross a set of taxa from GhostKOALA output files
and the KEGG KO annotation table (ko00001.keg).

"""
#General modules
import sys
import os
import re
import argparse
import pandas as pd
import scipy
import pylab
import scipy.cluster.hierarchy as sch

################################################################################

#Command-line usage with argparser

parser = argparse.ArgumentParser(prog='make_KO_KEGG_abundance_annotations_table.py',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	description='Outputs the abundance of KEGG KOs accross a set of taxa based on \
		output from GhostKOALA and the KEGG annotation table.')

parser.add_argument('-in', '--input_file', required=True, help = 'Path to the \
	input GhostKOALA file.')

parser.add_argument('-a', '--annotations', required=True, help = 'Path to the \
	file containing the KEGG annotations file ko00001.keg.')

parser.add_argument('-out', '--output_file', required=True, help = 'Path to the \
	output file which will be an abundance table of KOs accross the given taxa \
	alongside KEGG annotations.')

parser.add_argument('-taxa', '--taxa_order', required=True, help = 'Path to a \
	file listing the order of taxa for the final output.')

args = parser.parse_args()

if not os.path.exists(args.input_file):
	raise IOError("This file does not exist")

if not os.path.exists(args.annotations):
	raise IOError("This file does not exist")

################################################################################
#FUNCTIONS
def make_KEGG_dict(annotations):
	KEGG_dict = {}
	with open(annotations, 'r') as anno:
		for line in anno:
			if line.startswith("D"):
				line = " ".join(line.strip("\n").split())
				line = line.replace("D ","D\t").replace("; ",";\t").replace(" [EC:","\t[EC:")
				line = re.sub('(^D\tK[0-9]*) ', '\\1\t', line)
				KO = line.split("\t")[1]
				gene_name = line.split("\t")[2]
				if ";" in line:
					annotation = line.split("\t")[3].strip(";")
				else:
					annotation = ""
				if "\t[EC:" in line:
					EC_number = line.split("\t")[4]
				else:
					EC_number = ""
				if KO not in KEGG_dict.keys():
					KEGG_dict[KO] = [EC_number, gene_name, annotation]
	return KEGG_dict

def count_KOs_per_genome(KEGG_dict, input_file):
	taxa_KO_counts = {}
	with open(args.input_file, 'r') as input:
		for line in input:
			try:
				KO = line.split("\t")[1].strip("\n")
			except:
				pass
			else:
				taxon = ("@".join(line.split("@")[0:2]).replace("=","_"))
				if taxon not in taxa_KO_counts.keys():
					taxa_KO_counts[taxon] = []
				taxa_KO_counts[taxon].append(KO)
	return taxa_KO_counts

def make_KO_table_dict(KEGG_dict, taxa_KO_counts):
	KO_table_dict = {}
	for taxon in taxa_KO_counts.keys():
		KO_table_dict[taxon] = {}
		for KO in KEGG_dict.keys():
			KO_table_dict[taxon][KO] = {}
			KO_table_dict[taxon][KO] = taxa_KO_counts[taxon].count(KO)
			#print(KO_table_dict)
	return KO_table_dict


################################################################################
#IMPLEMENTATION

#Make dictionary of KEGG annotations
KEGG_dict = make_KEGG_dict(args.annotations)
KEGG_table = pd.DataFrame(KEGG_dict).transpose()
KEGG_table.columns = ['EC Number', 'Gene Name', 'KEGG KO Annotation']

#Make dictionary of  KOs for each taxon
taxa_KO_counts = count_KOs_per_genome(KEGG_dict, args.input_file)

#Make dictionary of number of each KO for each taxon and annotation
KO_table_dict = make_KO_table_dict(KEGG_dict, taxa_KO_counts)

abundance_table = pd.DataFrame(KO_table_dict)

#Concatenate dataframes
result = pd.concat([KEGG_table, abundance_table], axis=1)

#Reorder columns by tree order
column_names = ["EC Number","Gene Name", "KEGG KO Annotation"]
with open(args.taxa_order, 'r') as input:
	for line in input:
		taxon = line.strip("\n")
		column_names.append(taxon)

result = result.reindex(columns=column_names)

print(result)
result.to_csv(args.output_file, sep='\t')
