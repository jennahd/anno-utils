# anno-utils
Annotation utility scripts

Parse annotations from various annotation tools.

### "make_KO_KEGG_abundance_annotations_table.py"

For making a KEGG KO abundance table accross a set of taxa from GhostKOALA (https://www.kegg.jp/ghostkoala/, should also work with BlastKOALA, and KofamKOALA) output files
and the KEGG KO annotation table (ko00001.keg).

```
python make_KO_KEGG_abundance_annotations_table.py \
  -in "ghostKOALA-KOs.tsv" \
  -a "ko00001.keg" \
  -out "ghostKOALA-KOs.OVERVIEW.tsv" \
  -taxa "taxa_order.list"
```

### "make_NOG_abundance_table.py"

For making an overview of NOG abundances and a list of NOGs accross taxa based on emapper v2 (https://github.com/eggnogdb/eggnog-mapper) bestOG assignments with the eggNOG v.5 database (http://eggnog5.embl.de). Uses "universal-level" NOGs.

```
python make_NOG_abundance_table.py \
  -in "emapper.annotations" \
  -num "number" \
  -out "NOGs.list" \
  -tsv "NOGs.OVERVIEW.tsv" \
  -taxa "taxa_order.list"
```
- "num " assign the number of lineages in which the NOG must be present in order to consider it

### "make_NOG_abundance_table_v2.py"

For making an overview of NOG abundances and a list of NOGs accross taxa based on emapper v2 (https://github.com/eggnogdb/eggnog-mapper) bestOG assignments with the eggNOG v.5 database (http://eggnog5.embl.de). Uses NOGs assigned at the lowest taxonomic level.

```
python make_NOG_abundance_table.py \
  -in "emapper.annotations" \
  -num "number" \
  -out "NOGs.list" \
  -tsv "NOGs.OVERVIEW.tsv" \
  -taxa "taxa_order.list"
```
- "num " assign the number of lineages in which the NOG must be present in order to consider it

### "make_KO_abundance_table_emapper_annotations.py"

For making an overview of KO abundance accross taxa based on emapper v2 (https://github.com/eggnogdb/eggnog-mapper) bestOG assignments with the eggNOG v.5 database (http://eggnog5.embl.de). The KOs are transferred from the closest orthologs by the tool. This script will provide the sum of each KO accross all NOGs. Requires the KEGG KO annotation table (ko00001.keg).

```
python make_KO_abundance_table_emapper_annotations.py \
  -in "emapper.annotations" \
  -taxa "taxa_order.list" \
  -a "ko00001.keg" \
  -out "eggNOG-emapper_KOs.OVERVIEW.tsv"
```

### "make_annotations_table.py"

For making a table of annotations for each protein in a genome.

Inputs:

Protein fasta file:
File with protein sequences (".faa") from the given genome in fasta format.
The header ID should include genome_accession@species_name@protein_accession@protein_annotation.
The protein annotation could have been derived from the sequence in NCBI, or from prokka when calling proteins.

Annotation files:
Tools should have been run using the same protein sequences and header IDs as in the protein fasta file.
- Top DIAMOND blastp (https://github.com/bbuchfink/diamond) hit output (--outfmt 6 qseqid qlen sseqid slen stitle length qcovhsp pident evalue bitscore)
- GhostKOALA (https://www.kegg.jp/ghostkoala) KEGG KO output
- eggNOG emapper (https://github.com/eggnogdb/eggnog-mapper) annotations output (-m diamond)
- InterProScan (https://github.com/ebi-pf-team/interproscan) output (-appl Pfam, TIGRFAM)
- Requires the KEGG KO annotation table (ko00001.keg).

Output:
The result will be a tsv file with all annotations per protein sequence

```
  python make_annotations_table.py \
    -faa "faa" \
    -blastp "diamond-blastp.tsv" \
    -KOALA "ghostKOALA-KOs.tsv" \
    -emapper "emapper.annotations" \
    -interpro "interproscan.tsv" \
    -KEGGdesc "ko00001.keg" \
    -out "annotations.tsv"
```
