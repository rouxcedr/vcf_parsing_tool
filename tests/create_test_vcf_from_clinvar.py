import cyvcf2
import os
import sys
import argparse
from collections import defaultdict

def main(clinvar_file, output_directory):

    clinvar_vcf=cyvcf2.VCF(clinvar_file)

    info_fields =  [
        "AF_ESP",
        "AF_EXAC",
        "AF_TGP",
        "ALLELEID",
        "CLNDN",
        "CLNDNINCL",
        "CLNDISDB",
        "CLNDISDBINCL",
        "CLNHGVS",
        "CLNREVSTAT",
        "CLNSIG",
        "CLNSIGCONF",
        "CLNSIGINCL",
        "CLNVC",
        "CLNVCSO",
        "CLNVI",
        "DBVARID",
        "GENEINFO",
        "MC",
        "ORIGIN",
        "RS",
        "SSR",
    ]

    all_possible_values = defaultdict(set)

    for variant in clinvar_vcf:
        for field in info_fields:
            all_possible_values[field].add(variant.INFO.get(field))

    for field in info_fields:
        index = 1
        vcf_test = open(output_directory + field + "_vcf_test.vcf", "a")
        vcf_test.write(clinvar_vcf.raw_header)
        # print(all_possible_values[field])
        for value in all_possible_values[field]:
            vcf_test.write("chr1\t{pos}\t{id}\tA\tT\t.\t.\t{field}={value};\n".format(
                pos = index,
                id = index,
                field = field,
                value = value
            ))

            index +=1
        vcf_test.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("clinvar_file", help="File path to the clinvar cvf file")
    parser.add_argument("output_directory", help="Path to ouput directory")
    args = parser.parse_args()

    main(args.clinvar_file, args.output_directory)

