#!/usr/bin/env python

"""
Author: Jennah Dharamshi
Date: 220802

Script that gives an overview of KO abundance accross taxa based on
emapper v2 bestOG assignments with the eggNOG v.5 database. The KOs are
transferred from the closest orthologs by the tool. This script will provide
the sum of each KO accross all NOGs.

All taxa together in one file.
"""

#General modules
import sys
import os
import re
import argparse
import pandas as pd

################################################################################

#Command-line usage with argparser

parser = argparse.ArgumentParser(prog='overview_of_NOGs_abundance',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	description='Gives an overview of gene content based on emapper \
	NOG assignments present in a given number of taxa')

parser.add_argument('-in', '--input_file', required=True, help = 'Path to a \
	file listing the paths of all emapper annotations files to include.')

parser.add_argument('-taxa', '--taxa_order', required=True, help = 'Path to a \
	file listing taxa and their order.')

parser.add_argument('-a', '--annotations', required=True, help = 'Path to the \
	file containing the KEGG annotations file ko00001.keg.')

parser.add_argument('-out', '--output_file', required=True, help = 'Path for \
	the output abundance table of NOGs.')

args = parser.parse_args()

################################################################################

#Check input file exists
if not os.path.exists(args.input_file):
	raise IOError("This file does not exist")

################################################################################
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

def count_KOs_per_genome(input_file):
	taxa_KO_counts = {}
	with open(args.input_file, 'r') as input:
		for line in input:
			if "#" not in line:
				KOs = line.replace("ko:","").split("\t")[11].split(",")
				if KOs[0] != '-':
					taxon = "@".join(line.split("\t")[0].split("@")[0:2]).replace("=","_")
					if taxon not in taxa_KO_counts.keys():
						taxa_KO_counts[taxon] = []
					taxa_KO_counts[taxon].extend(KOs)
	return taxa_KO_counts

def make_KO_table_dict(KEGG_dict, taxa_KO_counts):
	KO_table_dict = {}
	for taxon in taxa_KO_counts.keys():
		KO_table_dict[taxon] = {}
		for KO in KEGG_dict.keys():
			KO_table_dict[taxon][KO] = {}
			KO_table_dict[taxon][KO] = taxa_KO_counts[taxon].count(KO)
	return KO_table_dict

###############################################################################

#IMPLEMENTATION

#Make dictionary of KEGG annotations
KEGG_dict = make_KEGG_dict(args.annotations)
KEGG_table = pd.DataFrame(KEGG_dict).transpose()
KEGG_table.columns = ['EC Number', 'Gene Name', 'KEGG KO Annotation']

#Get list of all taxa
taxa = []
with open(args.taxa_order, 'r') as f:
	for line in f:
		taxa.append(line.strip("\n"))

number_taxa = len(taxa)
print("Total number of taxa: %s" %number_taxa)

#Make dictionary of  KOs for each taxon
taxa_KO_counts = count_KOs_per_genome(args.input_file)

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
