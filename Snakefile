trait_tissue_nonmolecular = []
trait_tissue_molecular = []


with open("data/hukku_et_al/trait_tissue_nonmolecular.txt") as f:
    for line in f:
        trait_tissue_nonmolecular.append(line.strip())

with open("data/sinnott_armstrong_et_al/twas_coloc_sumstats/trait_tissue_molecular.txt") as f:
    for line in f:
        trait_tissue_molecular.append(line.strip())






rule all:
    input:
        "output/simulation/tables/table1.tab",
        "output/simulation/tables/tableS3.tab",
        "output/simulation/tables/tableS1.tab",
        "output/simulation/tables/tableS2.tab",
        "output/simulation/figures/figure3.pdf",
        "output/simulation/figures/figureS2.pdf",
        "output/simulation/figures/figureS3.pdf",
        "output/real_data_hukku/tables/table2.tab",
        "output/real_data_hukku/tables/tableS4.tab",
        expand("output/real_data_hukku/gene_set_enrichment_rst/{tt}_gene_set_enrichment.tab",tt = trait_tissue_nonmolecular),
        "output/real_data_hukku/tables/tableS5.tab",
        "output/real_data_hukku/supplement/Supplemental_Excel_File_1.tab",
        "output/real_data_sinnott_armstrong/bedtools/Urate_all_bedtools_output.bed",
        "output/real_data_sinnott_armstrong/bedtools/IGF_1_all_bedtools_output.bed",
        "output/real_data_sinnott_armstrong/bedtools/Testosterone_male_bedtools_output.bed",
        "output/real_data_sinnott_armstrong/bedtools/Testosterone_female_bedtools_output.bed",
        "output/real_data_sinnott_armstrong/intact_genes/intact_sig_genes.tab",
        "output/real_data_sinnott_armstrong/tables/table3.tab",
        "output/real_data_sinnott_armstrong/supplement/Supplemental_Excel_File_2.tab",
        "output/real_data_sinnott_armstrong/tables/table4.tab",
        "output/real_data_sinnott_armstrong/supplement/Supplemental_Excel_File_3.tab",
        "output/real_data_sinnott_armstrong/tables/table5.tab",
        "output/real_data_sinnott_armstrong/tables/tableS6.tab"





rule simulation_study:
    input:
        expand("data/simulated/alpha1_0/sim.phi_0_6.{sim_num}.summary",sim_num = list(range(1,101))),
        expand("data/simulated/alpha1_0_25/sim.phi_0_6.{sim_num}.summary",sim_num = list(range(1,101))),
        expand("data/simulated/alpha1_0_5/sim.phi_0_6.{sim_num}.summary",sim_num = list(range(1,101))),
        expand("data/simulated/alpha1_1/sim.phi_0_6.{sim_num}.summary",sim_num = list(range(1,101))),
        expand("data/simulated/alpha1_1_5/sim.phi_0_6.{sim_num}.summary",sim_num = list(range(1,101))),
        expand("data/simulated/with_alt_model/alpha1_0/sim.phi_0_6.{sim_num}.summary",sim_num = list(range(1,101))),
    output:
        "output/simulation/tables/table1.tab",
       	"output/simulation/tables/tableS3.tab",
        "output/simulation/tables/tableS1.tab",
        "output/simulation/tables/tableS2.tab",
        "output/simulation/figures/figure3.pdf",
        "output/simulation/figures/figureS2.pdf",
        "output/simulation/figures/figureS3.pdf"
    script:
        "scripts/simulation_analysis.R"



rule real_data_hukku:
    input:
        expand("data/hukku_et_al/{tt}.summary.out",tt = trait_tissue_nonmolecular),
        "data/gene_set/CTD_genes_pathways.csv",
        "data/gene_set/ensembl_bp_go_term_ref.tab",
        "data/biomart/biomart_ensembl_to_hgnc.tab"
    output:
        "output/real_data_hukku/tables/table2.tab",
        "output/real_data_hukku/tables/tableS4.tab",
    script:
        "scripts/hukku_real_data_analysis.R"	




rule go_gse:
    input:
        expand("data/hukku_et_al/{tt}.summary.out",tt = trait_tissue_nonmolecular),
        "data/gene_set/ensembl_bp_go_term_ref.tab",
    output:
        expand("output/real_data_hukku/gene_set_enrichment_rst/{tt}_gene_set_enrichment.tab",tt = trait_tissue_nonmolecular)
    threads: 8
    script:
        "scripts/gene_set_enrichment.R"




rule summarize_go_gse:
    input:
        expand("output/real_data_hukku/gene_set_enrichment_rst/{tt}_gene_set_enrichment.tab", tt= trait_tissue_nonmolecular)
    output:
        "output/real_data_hukku/tables/tableS5.tab",
        "output/real_data_hukku/supplement/Supplemental_Excel_File_1.tab"
    script:
        "scripts/summarize_go_gse.R"

         


rule real_data_sinnott_armstrong_bedtools_intersect:
    input:
        expand("data/sinnott_armstrong_et_al/twas_coloc_sumstats/{ttc}.summary.out",ttc = trait_tissue_molecular),
        "data/sinnott_armstrong_et_al/all_hits_msiggenes_ext_hg38.bed",
        "data/sinnott_armstrong_et_al/gwas_loci/Urate_all_hits_hg38.bed",
        "data/sinnott_armstrong_et_al/gwas_loci/IGF_1_all_hits_hg38.bed",
        "data/sinnott_armstrong_et_al/gwas_loci/Testosterone_male_hits_hg38.bed",
        "data/sinnott_armstrong_et_al/gwas_loci/Testosterone_female_hits_hg38.bed",
        "data/biomart/biomart_ensembl_to_hgnc.tab"
    output:
        "output/real_data_sinnott_armstrong/bedtools/Urate_all_bedtools_output.bed",
        "output/real_data_sinnott_armstrong/bedtools/IGF_1_all_bedtools_output.bed",
        "output/real_data_sinnott_armstrong/bedtools/Testosterone_male_bedtools_output.bed",
        "output/real_data_sinnott_armstrong/bedtools/Testosterone_female_bedtools_output.bed",
        "output/real_data_sinnott_armstrong/intact_genes/intact_sig_genes.tab" 
    script:
        "scripts/sinnott_armstrong_bedtools_intersect.R"




rule real_data_sinnott_armstrong:
    input:
        expand("data/sinnott_armstrong_et_al/twas_coloc_sumstats/{ttc}.summary.out",ttc = trait_tissue_molecular),
        "data/sinnott_armstrong_et_al/all_hits_msiggenes_ext_hg38.bed",
        "data/sinnott_armstrong_et_al/Supplementary_File_8.txt",
        "data/sinnott_armstrong_et_al/Supplementary_File_9.txt",
        "data/sinnott_armstrong_et_al/Supplementary_File_10.txt",
        "output/real_data_sinnott_armstrong/bedtools/Urate_all_bedtools_output.bed",
        "output/real_data_sinnott_armstrong/bedtools/IGF_1_all_bedtools_output.bed",
        "output/real_data_sinnott_armstrong/bedtools/Testosterone_male_bedtools_output.bed",
        "output/real_data_sinnott_armstrong/bedtools/Testosterone_female_bedtools_output.bed",
        "output/real_data_sinnott_armstrong/intact_genes/intact_sig_genes.tab",
        "data/sinnott_armstrong_et_al/gwas_loci/Supplementary_File_1.txt",
        "data/sinnott_armstrong_et_al/gwas_loci/Supplementary_File_2.txt",
        "data/sinnott_armstrong_et_al/gwas_loci/Supplementary_File_3.txt",
        "data/sinnott_armstrong_et_al/gwas_loci/Supplementary_File_4.txt",
        "data/biomart/biomart_ensembl_to_hgnc.tab"
    output:
        "output/real_data_sinnott_armstrong/tables/table3.tab",
        "output/real_data_sinnott_armstrong/supplement/Supplemental_Excel_File_2.tab",        
        "output/real_data_sinnott_armstrong/tables/table4.tab",
        "output/real_data_sinnott_armstrong/supplement/Supplemental_Excel_File_3.tab",
        "output/real_data_sinnott_armstrong/tables/table5.tab",
        "output/real_data_sinnott_armstrong/tables/tableS6.tab"
    script:
        "scripts/sinnott_armstrong_real_data_analysis.R"





