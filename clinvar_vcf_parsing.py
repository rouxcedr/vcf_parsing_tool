import pandas as pd
import cyvcf2
import argparse
import gzip
import os
import sys


def vcf_parse_with_info(vfc, filters):
    """
        Parses a vcf with filters from INFO fields

        Arguments : 
            vcf : a cyvcf2 VCF object initiated with the clinvar vcf to parse
            filters : a python dictionary with the info fields to parse as keys and a list of values to parse with as values
        
        Return:
            variants_filtered_list : the list of cyvcf2 variant objects that match the criterias of the filters
            gene_variant_count : a dictionary with genes symbol as keys and variant count for the gene as values
    """

    variants_filtered_list = list()
    gene_variant_count = {key: 0 for key in filters["GENEINFO"]}

    for variant in vcf:

        variant_passed_filtering = True
        variant_genes_list = list()

        for info_field, field_filter in filters.items():
            
            #Parsing GENEINFO Field, need to split with "|" and then with ":" to get (gene symbol, gene id)
            if info_field == "GENEINFO" and variant.INFO.get(info_field) is not None:
                variant_genes = variant.INFO.get(info_field).split("|")
                info_filter_passed = any(gene.split(":")[0] in field_filter for gene in variant_genes)

                #To count number of variants per gene
                for gene in variant_genes:
                    if gene.split(":")[0] in field_filter:
                        variant_genes_list.append(gene.split(":")[0])
          
                    
            #Parsing CLNREVSTAT field
            if info_field == "CLNREVSTAT" and variant.INFO.get(info_field) is not None:
                info_filter_passed = variant.INFO.get(info_field) in field_filter
            

            #Parsing CLNSIG field
            if info_field == "CLNSIG" and variant.INFO.get(info_field) is not None:
                clnsig_list = variant.INFO.get(info_field).split(",_")
                info_filter_passed = any(clnsig_value in field_filter for clnsig_value in clnsig_list)

            variant_passed_filtering = variant_passed_filtering and info_filter_passed

        if variant_passed_filtering:
            variants_filtered_list.append(variant)

            #Updating variant count for coresponding gene 
            for gene in variant_genes_list:
                gene_variant_count[gene] +=1

    return (variants_filtered_list, gene_variant_count)


def main(gene_list_file, clinvar_file, output_directory):
    """
        Arguments :
            gene_list_file : the file path to the file containing gene symbols
            clinvar_file : the file path to the clinvar vcf to parse
            output_directory : the file path to the output directory
    """
    if not output_directory.endswith("/"):
        output_directory = output_directory + "/"

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    #Reading Gene List
    with open(gene_list_file) as f:
        gene_list = f.read().splitlines()

    #Reading clinvar vcf
    clinvar_vcf = cyvcf2.VCF(clinvar_file)

    #Creating filters
    filters = {
        "GENEINFO": gene_list,
        "CLNREVSTAT": ["criteria_provided,_multiple_submitters,_no_conflicts", "reviewed_by_expert_panel", "practice_guideline"],
        "CLNSIG": ["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"]
    }

    #Parsing clinvar vcf
    clinvar_variants_filtered, gene_variant_count = vcf_parse_with_info(clinvar_vcf, filters)
    print("Number of variants : ", len(clinvar_variants_filtered))

    with open(output_directory + "clinvar_variant_count_per_genes.txt", "w") as f:
        for gene in sorted(gene_variant_count):
            f.write(gene + "  # " + str(gene_variant_count[gene]) + "\n")

    #Creating new vcf with filtered variants
    output_vcf = cyvcf2.Writer(output_directory + "clinvar_filtered.vcf", clinvar_vcf)
    for variant in clinvar_variants_filtered:
        output_vcf.write_record(variant)
    output_vcf.close()

if __name__ == "__main__":
    # enable command line arguments parsing. Easy to use and to add arguments. See : https://docs.python.org/fr/2/howto/argparse.html
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene_list_file", "-g", help="File path to the gene list")
    parser.add_argument("--clinvar_file", "-c", help="File path to the clinvar cvf file")
    parser.add_argument("--output_directory", "-o", help="File path to the output directory")
    args = parser.parse_args()

    main(args.gene_list_file, args.clinvar_file, args.output_directory)
