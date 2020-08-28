# rubisco_co_occurrence_gh
Pipeline to dowload, align, and analyze rubisco isoform sequences for phylogeny and metabolic marker co-occurrence analysis. 

Pipeline includes R scripts, Python notebooks, the CIPRES portal, and code previously published by external developers: BOOSTER, Tree2Fasta, and OxyPhen.

PIPELINE OUTLINE:

- Python notebook to download UniprotKB sequences:	"2019_4_5_Download_uniprot_API_with_phylum_and_class_and_fragment_rubisco.ipynb"
- R script to process RuBisCO seq download: 		"2019_2_25_Process_Avi_uniprotKB_downloads.R"
- R script to screen out incomplete sequences, etc:	"2019_3_21_scraping_out_environmentals_from_rubisco_lg_csv.R"
- CIPRES portal to creat alignment and tree of RuBisCO sequences: http://www.phylo.org/
- Calculate transfer bootstrapping of tree using BOOSTER: https://booster.pasteur.fr/
- FigTree to manually annotate clades according to reference RuBisCO isoforms: http://tree.bio.ed.ac.uk/software/
- Tree2Fasta to recover lists of annotated sequences: https://pubmed.ncbi.nlm.nih.gov/29506565/
- R scripts to process Tree2Fasta output: 			"2019_3_1_sort_uniprot_by_tree2fasta output.R"
- Python notebook for co-occurrence, taxonomy barchart: "2019_10_19_conditional_probability_and_updated_markers_and_phylogeny_bar_charts.ipynb"

Additional analysis: OxyPhen, to predict O2-utilizing markers per proteome of each Rubisco-containing taxon.
- Python notebook to download proteome files: "2019_4_16_OxyPhen_Download_Uniprot_API.ipynb"
- R script to generate formatted FASTA files: "2019_4_16_loop_for_formatting_oxyphen_proteome_csv_files_to_fasta.R"
- Python notebook to perform Oxyphen analysis: "2019_5_2_Troubleshoot_Oxyphen_multiome.ipynb"
- R script to merge Oxyphen output with taxonomy:  "2019_5_7_oxyphen_output_to_phylogeny.R"
- R scrpt to graph oxyphen CDF by isoform:		"2019_7_15_oxyphen_graphing_by_cdf_ggplot.R"
