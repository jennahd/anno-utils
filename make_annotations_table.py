#!/usr/bin/env python

"""
Author: Jennah Dharamshi
Date: 221006

Script that puts together a table of annotations for each protein ORF in a genome.

Inputs:

Protein fasta file:
File with protein sequences from the given genome in fasta format.
The header ID should include genome_accession@species_name@protein_accession@protein_annotation.
The protein annotation could have been derived from the sequence in NCBI, or from prokka when calling proteins.

Annotation files:
Tools should have been run using the same protein sequences and header IDs as in the protein fasta file.
- Top diamond blastp hit output (--outfmt 6 qseqid qlen sseqid slen stitle length qcovhsp pident evalue bitscore)
- GhostKOALA KEGG KO output
- eggNOG emapper annotations output (-m diamond)
- Interproscan output (-appl Pfam, TIGRFAM)

Output:
The result will be a tsv file with all annotations per protein sequences
"""

#General modules
import sys
import os
import re
import argparse
import pandas as pd
from Bio import SeqIO

################################################################################

#Command-line usage with argparser

parser = argparse.ArgumentParser(prog='make_annotations_table.py',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	description='Will output an overview table of annotations per protein sequence \
	from the provided blastp, ghostKOALA, eggNOG emapper, and interproscan output.')

parser.add_argument('-faa', '--faa', required=True, help = 'Path to the \
	fasta file with protein sequences and information in headers.')

parser.add_argument('-blastp', '--blastp', required=True, help = 'Path to the \
	blastp output file.')

parser.add_argument('-KOALA', '--KOALA', required=True, help = 'Path to the \
	ghost/blastKOALA output file.')

parser.add_argument('-emapper', '--emapper', required=True, help = 'Path to the \
	eggNOG emapper annotations output file.')

parser.add_argument('-interpro', '--interpro', required=True, help = 'Path to the \
	Interproscan output file.')

parser.add_argument('-KEGGdesc', '--KEGGdesc', required=True, help = 'Path to the \
	file containing the KEGG KO descriptions file ko00001.keg.')

parser.add_argument('-out', '--output_file', required=True, help = 'Path for \
	the output table of annotations per protein.')

args = parser.parse_args()

################################################################################

#Check input files exist
if not os.path.exists(args.faa):
	raise IOError("This file does not exist")
if not os.path.exists(args.blastp):
	raise IOError("This file does not exist")
if not os.path.exists(args.KOALA):
	raise IOError("This file does not exist")
if not os.path.exists(args.emapper):
	raise IOError("This file does not exist")
if not os.path.exists(args.interpro):
	raise IOError("This file does not exist")

################################################################################
def make_KEGG_dict(KEGGdesc):
	KEGG_dict = {}
	with open(KEGGdesc, 'r') as anno:
		for line in anno:
			if line.startswith("D"):
				line = " ".join(line.strip("\n").split())
				line = line.replace("D ","D\t").replace("; ",";\t").replace(" [EC:","\t[EC:")
				line = re.sub('(^D\tK[0-9]*) ', '\\1\t', line)
				KO = line.split("\t")[1]
				gene_name = line.split("\t")[2].strip(";").replace(", ",",")
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

def make_protein_dict(faa):
	prot_list = []
	protein_dict = {}
	for seq_record in SeqIO.parse(faa, "fasta"):
		ID = seq_record.description
		prot_list.append(ID.split(" ")[0])
		genome, species, protein, prot_desc = ID.split("@",3)
		protein_dict[ID.split(" ")[0]] = genome, species, protein, prot_desc
	return prot_list, protein_dict

def make_blastp_dict(blastp, prot_list):
	blastp_dict = {}
	for hit in blastp:
		hit = hit.split("\t")[0:9]
		ID, length, acc, qcovhsp, pident, evalue = hit[0], hit[1], hit[2], hit[6], hit[7], hit[8]
		if "[" in hit[4]:
			species, desc = hit[4].split("[")[1].split("]")[0], hit[4].split("[")[0].split(" ",1)[1:][0]
		else:
			species, desc = "", hit[4].split(" ",1)[1:][0]
		blastp_dict[ID] = length, acc, species, desc, qcovhsp, pident, evalue
	for ID in prot_list:
		if ID not in blastp_dict.keys():
			blastp_dict[ID] = "","","","","","",""
	return blastp_dict

def make_KOALA_dict(KOALA, KEGG_dict):
	KOALA_dict = {}
	for hit in KOALA:
		try:
			ID = hit.split("\t")[0]
			KO = hit.split("\t")[1].strip("\n")
		except:
			ID = hit.strip("\n")
			KO, EC_number, gene_name, annotation = "", "", "", ""
		else:
			try:
				EC_number, gene_name, annotation = KEGG_dict[KO][0], KEGG_dict[KO][1], KEGG_dict[KO][2]
			except:
				KO, EC_number, gene_name, annotation = KO, "", "", ""
		KOALA_dict[ID] = KO, EC_number, gene_name, annotation
	return KOALA_dict

def make_emapper_dict(emapper, prot_list, KEGG_dict):
	emapper_dict = {}
	for hit in emapper:
		hit = hit.split("\t")
		ID, NOG_list, min_annot_lvl, COG_cat, Desc, Name, EC = hit[0], hit[4], hit[5], hit[6], hit[7], hit[8], hit[10]
		KOs = hit[11].replace("ko:","").split(",")
		ECs, genes, annos = [], [], []
		if KOs[0] != "-":
			for KO in KOs:
				try:
					ECs.append(KEGG_dict[KO][0]), genes.append(KEGG_dict[KO][1]), annos.append(KEGG_dict[KO][2])
				except:
					ECs.append(""), genes.append(""), annos.append("")
			ECs, genes, annos = ';'.join(ECs), ';'.join(genes).replace(", ",","), ';'.join(annos)
		else:
			ECs, genes, annos = "", "", ""
		KOs = ';'.join(KOs)
		emapper_dict[ID] = COG_cat, Name, EC, min_annot_lvl, Desc, NOG_list, KOs, ECs, genes, annos
	for ID in prot_list:
		if ID not in emapper_dict.keys():
			emapper_dict[ID] = "","","","","","","","","",""
	return emapper_dict

def make_interpro_dict(interpro, prot_list):
	interpro_dict = {}
	for hit in interpro:
		hit = hit.split("\t")[0:6]
		ID = hit[0]
		if ID not in interpro_dict.keys():
			interpro_dict[ID] = [], [], [], []
		if hit[3] == "Pfam":
			interpro_dict[ID][0].append(hit[4]), interpro_dict[ID][1].append(hit[5])
		elif hit[3] == "TIGRFAM":
			interpro_dict[ID][2].append(hit[4]), interpro_dict[ID][3].append(hit[5])
	for ID in prot_list:
		if ID in interpro_dict.keys():
			interpro_dict[ID] = ';'.join(interpro_dict[ID][0]),';'.join(interpro_dict[ID][1]),';'.join(interpro_dict[ID][2]),';'.join(interpro_dict[ID][3])
		else:
			interpro_dict[ID] = "","","",""
	return interpro_dict

###############################################################################

#IMPLEMENTATION

#Open files
faa = open(args.faa, 'r')
blastp = open(args.blastp, 'r')
KOALA = open(args.KOALA, 'r')
emapper = open(args.emapper, 'r')
interpro = open(args.interpro, 'r')

#Make dictionary of KEGG annotations
KEGG_dict = make_KEGG_dict(args.KEGGdesc)

#Make list of protein IDs and dictionary of info
prot_list, protein_dict = make_protein_dict(faa)
protein_dict = pd.DataFrame(protein_dict).transpose()

#Make dictionary of blastp info:
blastp_dict = make_blastp_dict(blastp, prot_list)
blastp_dict = pd.DataFrame(blastp_dict).transpose()

#Make dictionary of KOALA info:
KOALA_dict = make_KOALA_dict(KOALA, KEGG_dict)
KOALA_dict = pd.DataFrame(KOALA_dict).transpose()

#Make dictionary of emapper info:
emapper_dict = make_emapper_dict(emapper, prot_list, KEGG_dict)
emapper_dict = pd.DataFrame(emapper_dict).transpose()

#Make dictionary of interpro info:
interpro_dict = make_interpro_dict(interpro, prot_list)
interpro_dict = pd.DataFrame(interpro_dict).transpose()

#Concatenate dataframes
annotation_table = pd.concat([protein_dict, blastp_dict, KOALA_dict, emapper_dict, interpro_dict], axis=1)
annotation_table = annotation_table.replace("-", "")
annotation_table = annotation_table.replace(";", "")

column_names =["Genome_ID","Species_name","Protein_ID","NCBI_prokka_annotation",
				"protein_length","blastp_hit_acc","blastp_hit_species","blastp_hit_desc","blastp_hit_qcovhsp","blastp_hit_pident","blastp_hit_evalue",
				"ghostKOALA_KO","ghostKOALA_KO_EC_number","ghostKOALA_KO_gene_name","ghostKOALA_KO_description",
				"emapper_COG_cat","emapper_gene_name","emapper_EC","emapper_min_annot_lvl","emapper_description","emapper_all_NOGs","emapper_NOG_KO","emapper_NOG_KO_EC","emapper_NOG_KO_gene_name","emapper_NOG_KO_description",
				"interproscan_Pfam_domains","interproscan_Pfam_descriptions","interproscan_TIGRFAM_domains","interproscan_TIGRFAM_descriptions"]
annotation_table.columns = column_names

#Output annotation table
annotation_table.to_csv(args.output_file, sep='\t', index_label='faa_ID')
