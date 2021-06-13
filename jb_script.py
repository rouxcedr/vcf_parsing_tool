# Read list of genes:
genes_file = "gene_list.txt"
genes_list = [line.strip() for line in open(genes_file)]
# Define ClinVar functions (using cyvcf2):

def clinvar_genes(variant):
    genes = []
    if variant.INFO.get("GENEINFO") is not None:
        geneinfo = variant.INFO.get("GENEINFO").split("|")
        for info in geneinfo:
            info = info.split(":")
            gene = info[0]
            genes.append(gene)
    return genes

def clinvar_stars(variant):
    stars = 0
    if variant.INFO.get("CLNREVSTAT") is not None:
        clnrevstat = variant.INFO.get("CLNREVSTAT")
        if clnrevstat == "criteria_provided,_single_submitter":
            stars = 1
        if clnrevstat == "criteria_provided,_conflicting_interpretations":
            stars = 1
        if clnrevstat == "criteria_provided,_multiple_submitters,_no_conflicts":
            stars = 2
        if clnrevstat == "reviewed_by_expert_panel":
            stars = 3
        if clnrevstat == "practice_guideline":
            stars = 4
    return stars

def clinvar_clnsig(variant):
    rec = {"Affects":0, "Benign":0, "Benign/Likely_benign":0, "Conflicting_interpretations_of_pathogenicity":0, "Likely_benign":0, "Likely_pathogenic":0, "Pathogenic":0, "Pathogenic/Likely_pathogenic":0, "Uncertain_significance":0, "association":0, "association_not_found":0, "confers_sensitivity":0, "drug_response":0, "not_provided":0, "other":0, "protective":0,  "risk_factor":0}
    if variant.INFO.get("CLNSIG") is not None:
        clnsig = variant.INFO.get("CLNSIG")
        clnsig = clnsig.split(",")
        for sig in clnsig:
            if sig.startswith("_"):
                sig = sig[1:]
            rec[sig] += 1
    return rec

def clinvar_sigconf(variant):
    rec = {"Benign":0, "Likely_benign":0, "Likely_pathogenic":0, "Pathogenic":0, "Uncertain_significance":0}
    if variant.INFO.get("CLNSIGCONF") is not None:
        clnsigconf = variant.INFO.get("CLNSIGCONF")
        clnsigconf = clnsigconf.split(",")
        for sig in clnsigconf:
            sig = sig.split("(")[0]
            rec[sig] += 1
    return rec
    
# Get BED file of all ClinVar variants (clinvar_variants.bed):
from cyvcf2 import VCF
clinvar_vcf = "clinvar_20210529.vcf.gz"
outfile = open("clinvar_variants.bed", "w")
for variant in VCF(clinvar_vcf):
    in_gene = False
    one_star = False
    # Check if variant in genes_list provided:
    for gene in clinvar_genes(variant):
        if gene in genes_list:
            in_gene = True
    # Check if variant has at least 1 star:
    if in_gene is True:
        if clinvar_stars(variant) >= 1:
            one_star = True
    # Check variant classification:
    if in_gene is True and one_star is True:
        clnsig = clinvar_clnsig(variant)
        sigconf = clinvar_sigconf(variant)
        clinsig_count = clnsig["Pathogenic"] + clnsig["Pathogenic/Likely_pathogenic"] + clnsig["Likely_pathogenic"] + clnsig["drug_response"] + clnsig["risk_factor"]
        sigconf_count = sigconf["Likely_pathogenic"] + sigconf["Pathogenic"]
        if clinsig_count + sigconf_count >= 1:
            chrom = "chr" + variant.CHROM
            start = str(variant.start)
            end = str(variant.end)
            name = "clinvar_id:" + variant.ID
            line = [chrom, start, end, name]
            # Print variants < 500 bp to bed file:
            if variant.end - variant.start < 500:
                print("\t".join(line), file=outfile)
            # Calculate intervals for variants > 500 bp tp bed file (only keep 100 bp after start and before end):
            else:
                start_prox = str(variant.start)
                end_prox = str(variant.start + 100)
                start_dist = str(variant.end - 100)
                end_dist = str(variant.end)
                # Print variants >= 500 bp to bed file:
                line_prox = [chrom, start_prox, end_prox, name]
                line_dist = [chrom, start_dist, end_dist, name]
                print("\t".join(line_prox), file=outfile)
                print("\t".join(line_dist), file=outfile)
outfile.close()