#!/usr/bin/env python

"""
Author: Jennah Dharamshi
Date: 220617

Script that gives an overview of NOG abundance accross taxa based on emapper v2
bestOG assignments with the eggNOG v.5 database.

Must assign the number of lineages in which the NOG must be present in.

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

parser.add_argument('-taxa', '--taxa_list', required=True, help = 'Path to a \
	file listing taxa and their order.')

parser.add_argument('-num', '--number', required=True, help = 'Number \
	of taxa that the NOG must be in.')

parser.add_argument('-out', '--output_file', required=True, help = 'Path for the \
	output file containing NOGs')

parser.add_argument('-tsv', '--absence_presence', required=True, help = 'Path for \
	the output abundance table of NOGs.')

args = parser.parse_args()

################################################################################

#Check input file exists
if not os.path.exists(args.input_file):
	raise IOError("This file does not exist")

#Output files
NOGs_out = open(args.output_file,'w')
presence_out = open(args.absence_presence,'w')

#Assign percentage cutoff to variable
taxa_to_pass = int(args.number)

################################################################################

def get_NOG_list(input_file):
	"Get list of all NOGs found across taxa."
	all_NOGs = []
	taxon_NOGs = {}
	with open(input_file, 'r') as f:
		for line in f:
			if "#" not in line:
				acc_taxon = "@".join(line.split("\t")[0].split("@")[0:2]).replace("=","_")
				NOG = line.split("\t")[4].rsplit(',', 1)[0].split("@")[0]
				if acc_taxon not in taxon_NOGs.keys():
					taxon_NOGs[acc_taxon] = []
				if NOG not in taxon_NOGs[acc_taxon]:
					taxon_NOGs[acc_taxon].append(NOG)
	for taxon in taxon_NOGs.keys():
		for NOG in taxon_NOGs[taxon]:
			all_NOGs.append(NOG)
	return all_NOGs

def get_passed_NOG_list(all_NOGs, taxa_to_pass):
	"Get list of NOGs present in at the least the number of taxa necessary \
	to be further considered"
	all_NOGs_set = set(all_NOGs)
	passed_NOGs = []
	print("There are a total of %s UNIQUE NOGs found accross all taxa" %(len(all_NOGs_set)))
	for NOG in all_NOGs_set:
		if all_NOGs.count(NOG) >= taxa_to_pass:
			passed_NOGs.append(NOG)
	return passed_NOGs

def make_abundance_table(passed_NOGs, input_file, taxa):
	taxa_NOG_dict = {}
	with open(args.input_file, 'r') as f:
		for line in f:
			if "#" not in line:
				acc_taxon = "@".join(line.split("\t")[0].split("@")[0:2]).replace("=","_")
				NOG = line.split("\t")[4].rsplit(',', 1)[0].split("@")[0]
				if acc_taxon not in taxa_NOG_dict.keys():
					taxa_NOG_dict[acc_taxon] = {}
					for NOG in passed_NOGs:
						taxa_NOG_dict[acc_taxon][NOG] = 0
				if NOG in passed_NOGs:
					taxa_NOG_dict[acc_taxon][NOG] += 1
	return taxa_NOG_dict

def NOG_info(passed_NOGs, input_file):
	NOG_annotations = {}
	with open(input_file, 'r') as f:
		for line in f:
			if "#" not in line:
				NOG = line.split("\t")[4].rsplit(',', 1)[0].split("@")[0]
				cat = line.split("\t")[6]
				desc = line.split("\t")[7]
				gene = line.split("\t")[8]
				if NOG in passed_NOGs:
					if NOG not in NOG_annotations.keys():
						NOG_annotations[NOG] = cat + "\t" + desc + "\t" + gene
	return NOG_annotations

def output_NOGs(passed_NOGs, NOGs_out, NOG_annotations):
	for NOG in passed_NOGs:
		NOGs_out.write("%s\t%s\n" % (NOG, NOG_annotations[NOG]))

def output_abundance_table(taxa_NOG_dict, presence_out, taxa):
	absence_presence_table = pd.DataFrame(taxa_NOG_dict)
	print(absence_presence_table)
	absence_presence_table = absence_presence_table.reindex(columns=taxa)
	print(absence_presence_table)
	absence_presence_table.to_csv(presence_out,sep='\t')
################################################################################

#IMPLEMENTATION

#Get list of all NOGs accross all lineages
#taxa, number_taxa = get_taxa_list(args.input_file)
taxa = []
with open(args.taxa_list, 'r') as f:
	for line in f:
		taxa.append(line.strip("\n"))

number_taxa = len(taxa)
print("Total number of taxa: %s" %number_taxa)
print("Number of taxa NOG must be found in: %s" %taxa_to_pass)

all_NOGs = get_NOG_list(args.input_file)
print("There are a total of %s NOGs found accross all taxa" %(len(all_NOGs)))

passed_NOGs = get_passed_NOG_list(all_NOGs, taxa_to_pass)
print("There are a total of %s NOGs found in at least %s taxa" %((len(passed_NOGs)),taxa_to_pass))

#Make abundance table
taxa_NOG_dict = make_abundance_table(passed_NOGs, args.input_file, taxa)

#Get NOG information
NOG_annotations = NOG_info(passed_NOGs, args.input_file)
print(NOG_annotations)

#Outputs
output_NOGs(passed_NOGs, NOGs_out, NOG_annotations)
output_abundance_table(taxa_NOG_dict, presence_out, taxa)
