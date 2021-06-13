import pandas as pd
import argparse
import os
import sys

def main(gene_list_file, refseq_file, refseq_to_gene_symbol_file, reference_genome_index, output_directory):

    """
        Parses a refseq bed to ouput a bed that include only the exons from genes in a gene list

        Arguments :

            gene_list_file : The file path to the text file containing gene symbol to parse with
            refseq_file : The file path to the refseq bed containing all the exons
            refseq_to_gene_symbol : The file path to the file containing refseq IDs and there corresponding gene symbol
            reference_genome_index : The file path to the fasta index file of the reference genome 
            output_directory : Path to the ouput directory

        Ouput files :

            A bed containing all the exons from the refseq file of the genes in gene list file.
            filename : refseq_coding_exons_filtered.bed 

            A bed containing all the exons from the refseq file of the genes in gene list file without any alternative contigs.
            filename : refseq_coding_exons_filtered_no_alt_contigs.bed

    """

    if not output_directory.endswith("/"):
        output_directory = output_directory + "/"

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    #Reading Gene List
    with open(gene_list_file) as f:
        gene_list = f.read().splitlines()

    #Reading gene ID conversion
    refeseq_id_to_gene_symbol = pd.read_csv(refseq_to_gene_symbol_file, sep="\t", names=["REFSEQ_ID", "GENE_SYMBOL"])
    refeseq_id_to_gene_symbol = refeseq_id_to_gene_symbol.drop_duplicates(["REFSEQ_ID"])

    #Reading Refseq with coding exons
    refseq = pd.read_csv(refseq_file, sep = "\t", names=["CHROM", "START", "END", "NAME", "SCORE", "STRAND"])

    #Testing if all genes in gene_list has a coresponding refseq ID
    gene_list_has_refseq_id = all(gene_symbol in refeseq_id_to_gene_symbol["GENE_SYMBOL"].to_list() for gene_symbol in gene_list)

    if not gene_list_has_refseq_id:
        gene_symbol_without_refseq_id = [gene_symbol for gene_symbol in gene_list if gene_symbol not in refeseq_id_to_gene_symbol["GENE_SYMBOL"].to_list()]
        print("The following genes doesnt have a corresponding refseq ID in the refseq ID to gene symbol file : {}".format(gene_symbol_without_refseq_id), file=sys.stderr)
        sys.exit(1)

    #Getting refseq IDs from name of coding exon (name needs to be in format refseq_id_cds_xxxxxx)
    refseq["REFSEQ_ID"] = refseq["NAME"].apply(lambda x: x.split("_cds")[0])

    #Merging gene symbol on refseq_id
    refseq_merged = refseq.merge(refeseq_id_to_gene_symbol, left_on="REFSEQ_ID", right_on="REFSEQ_ID", how="left")

    #Selecting coding exon of genes in gene_list
    refseq_filtered = refseq_merged[refseq_merged["GENE_SYMBOL"].isin(gene_list)]

    #Testing if all genes in gene list have coding exons in refseq
    genes_have_coding_exon_in_refseq = all(gene_symbol in refseq_filtered["GENE_SYMBOL"].to_list() for gene_symbol in gene_list)

    if not genes_have_coding_exon_in_refseq:
        genes_without_coding_exon_in_refseq = [gene_symbol for gene_symbol in gene_list if gene_symbol not in refseq_filtered["GENE_SYMBOL"].to_list()]
        print("The following genes doesnt have a coding exon in refseq : {}".format(genes_without_coding_exon_in_refseq), file=sys.stderr)
        sys.exit(1)
    
    #Testing if all coding exons of a gene in refseq is in filtered dataframe
    refeseq_id_to_gene_symbol_in_refseq = refeseq_id_to_gene_symbol[refeseq_id_to_gene_symbol["REFSEQ_ID"].isin(refseq["REFSEQ_ID"])]
    refseq_id_to_gene_symbol_in_gene_list = refeseq_id_to_gene_symbol_in_refseq[refeseq_id_to_gene_symbol_in_refseq["GENE_SYMBOL"].isin(gene_list)]
    genes_have_all_coding_exon_from_refseq = all(refseq_id in refseq_filtered["REFSEQ_ID"].to_list() for refseq_id in refseq_id_to_gene_symbol_in_gene_list["REFSEQ_ID"])

    if not genes_have_all_coding_exon_from_refseq:
        genes_without_refseq_coding_exon = [refseq_id for refseq_id in refseq_id_to_gene_symbol_in_gene_list["REFSEQ_ID"] if refseq_id not in refseq_filtered["REFSEQ_ID"].to_list()]
        print("The following exons where not filtered in and are coding exons for genes in gene list found in : {}".format(genes_without_refseq_coding_exon), file=sys.stderr)
        sys.exit(1)


    #Keeping only contigs that are in the reference genome index
    reference_genome_index = pd.read_csv(reference_genome_index, sep="\t", names=["CHROM", "LENGTH", "OFFSET", "LINEBASES", "LINEWIDTH"])
    refseq_filtered = refseq_filtered[refseq_filtered["CHROM"].isin(reference_genome_index["CHROM"])]

    #Creating a dataframe without alternative contigs
    refseq_filtered_no_alternativ_contigs = refseq_filtered[~refseq_filtered["CHROM"].str.endswith("_alt")]

    #Writing parsed refseq bed
    refseq_filtered.to_csv(output_directory + "refseq_coding_exons_filtered.bed", sep="\t", columns=["CHROM", "START", "END", "NAME"], index=False, header=False)
    refseq_filtered_no_alternativ_contigs.to_csv(output_directory + "refseq_coding_exons_filtered_no_alt_contigs.bed", sep="\t", columns=["CHROM", "START", "END", "NAME"], index=False, header=False)

if __name__ == "__main__":
    # enable command line arguments parsing. Easy to use and to add arguments. See : https://docs.python.org/fr/2/howto/argparse.html
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene_list_file", "-g", help="File path to the gene list")
    parser.add_argument("--refseq_bed_file", "-r", help="File path to the refseq bed file")
    parser.add_argument("--refeseq_id_to_gene_symbol", "-id", help="File path to the file with refseq IDs and its corresponding gene symbol")
    parser.add_argument("--reference_genome_index", "-fai", help="File path to reference genome fai file")
    parser.add_argument("--output_directory", "-o", help="File path to the output directory")
    args = parser.parse_args()

    main(args.gene_list_file, args.refseq_bed_file, args.refeseq_id_to_gene_symbol, args.reference_genome_index, args.output_directory)