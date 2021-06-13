import unittest
import sys
import pandas as pd

sys.path.append("../")
import clinvar_vcf_parsing

class TestVCFFiltering(unittest.TestCase):
    vcf_test = pd.read_csv("test_vcf.vcf")

    def test_vcf_geneinfo(self):
        result_tardbp = db_filtering.vcf_gene_info(vcf_test, ["TARDBP"])
        result_tp53_mthfr = db_filtering.vcf_gene_info(vcf_test, ["TP53", "MTHFR"])

    def test_vcf_clnrevstat(self):
        result_4_stars = db_filtering.vcf_clnrvstat(vcf_test, ["practice_guideline"])
        result_1_star = db_filtering.vcf_clnrvstat(vcf_test, ["criteria_provided,_conflicting_interpretations", "criteria_provided,_single_submitter"])
        result_more_than_2_stars = db_filtering.vcf_clnrvstat(vcf_test, ["practice_guideline", "reviewed_by_expert_panel", 	"criteria_provided,_multiple_submitters,_no_conflicts"])

    def test_vcf_clnsig(self):
        result_pathogenic = db_filtering.clnsig(vcf_test, ["Pathogenic"])
        result_pathogenic_and_likely_pathogenic = db_filtering.clnsig(vcf_test, ["Pathogenic", "Likely_pathogenic"])