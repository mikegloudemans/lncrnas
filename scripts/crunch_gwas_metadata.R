require(dplyr)

### Get number of GWAS hits at various thresholds, for each study

gwas_data = read.table("/oak/stanford/groups/smontgom/mgloud/projects/lncrnas/output/snp-lists/lncrnas-tier1_all-gwas_source-pval1e-05_source-window500000_snps-considered.txt", header=FALSE, skip=1)
colnames(gwas_data) = c("chr", "snp_pos", "pvalue", "source_trait", "source_file")

gwas_data$short_trait = substring(gwas_data$source_trait, 33)

num_gwas_hits = gwas_data %>% group_by(short_trait) %>% summarize(p1e_5 = sum(pvalue < 1e-5), p1e_6 = sum(pvalue < 1e-6), p1e_7 = sum(pvalue < 1e-7), p5e_8 = sum(pvalue < 5e-8), p1e_10 = sum(pvalue < 1e-10), p1e_15 = sum(pvalue < 1e-15))
print(as.data.frame(num_gwas_hits))
print(colSums(num_gwas_hits[,2:7]))
write.table(num_gwas_hits, "/oak/stanford/groups/smontgom/mgloud/projects/lncrnas/output/gwas_tables/num_gwas_hits.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)


### Get number of GWAS-GWAS overlaps at various thresholds, for each study

#
# Make table, for both lookup window size 0 and 10,000
#
for (window in c("0", "10000"))
{
	overlap_data = read.table(sprintf("/oak/stanford/groups/smontgom/mgloud/projects/lncrnas/output/snp-lists/lncrnas-tier1_all-gwas_gtex-ir_source-pval1e-05_lookup-pval1e-05_source-window500000_lookup-window%s_coloc-tests.txt", window), header=TRUE)
	overlap_data$short_trait = substring(overlap_data$source_trait, 33)

	overlap_eqtl_data = overlap_data %>% filter(grepl("eQTL", lookup_file))
	overlap_sqtl_data = overlap_data %>% filter(grepl("sQTL", lookup_file))

	# For now, we're just going to assume using a constant threshold of 1e-5 in
	# the lookup set, but we could easily vary that if necessary. Would just be another loop

	# All overlaps
	num_overlaps = overlap_data %>% group_by(short_trait) %>% summarize(p1e_5 = sum(source_pvalue < 1e-5), p1e_6 = sum(source_pvalue < 1e-6), p1e_7 = sum(source_pvalue < 1e-7), p5e_8 = sum(source_pvalue < 5e-8), p1e_10 = sum(source_pvalue < 1e-10), p1e_15 = sum(source_pvalue < 1e-15))
	write.table(num_overlaps, sprintf("/oak/stanford/groups/smontgom/mgloud/projects/lncrnas/output/gwas_tables/num_window%s_overlaps.tsv", window), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
	print(num_overlaps)
	print(colSums(num_overlaps[,2:7]))

	# eQTL only
	num_eqtl_overlaps = overlap_eqtl_data %>% group_by(short_trait) %>% summarize(p1e_5 = sum(source_pvalue < 1e-5), p1e_6 = sum(source_pvalue < 1e-6), p1e_7 = sum(source_pvalue < 1e-7), p5e_8 = sum(source_pvalue < 5e-8), p1e_10 = sum(source_pvalue < 1e-10), p1e_15 = sum(source_pvalue < 1e-15))
	write.table(num_eqtl_overlaps, sprintf("/oak/stanford/groups/smontgom/mgloud/projects/lncrnas/output/gwas_tables/num_window%s_eqtl_overlaps.tsv", window), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
	print(num_eqtl_overlaps)
	print(colSums(num_eqtl_overlaps[,2:7]))

	# sQTLs only
	num_sqtl_overlaps = overlap_sqtl_data %>% group_by(short_trait) %>% summarize(p1e_5 = sum(source_pvalue < 1e-5), p1e_6 = sum(source_pvalue < 1e-6), p1e_7 = sum(source_pvalue < 1e-7), p5e_8 = sum(source_pvalue < 5e-8), p1e_10 = sum(source_pvalue < 1e-10), p1e_15 = sum(source_pvalue < 1e-15))
	write.table(num_sqtl_overlaps, sprintf("/oak/stanford/groups/smontgom/mgloud/projects/lncrnas/output/gwas_tables/num_window%s_sqtl_overlaps.tsv", window), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
	print(num_sqtl_overlaps)
	print(colSums(num_sqtl_overlaps[,2:7]))

	# Make summary of totals
	summary = do.call(rbind, list(colSums(num_overlaps[,2:7]), colSums(num_eqtl_overlaps[,2:7]), colSums(num_sqtl_overlaps[,2:7])))
	rownames(summary) = c("all-QTL", "eQTL-only", "sQTL-only")
	write.table(summary, sprintf("/oak/stanford/groups/smontgom/mgloud/projects/lncrnas/output/gwas_tables/window%s_overlap_summary.tsv", window), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
}

# ^ This will just give the number of tests though, including many that are the same
# locus but just with different tissues. One could also ask how many 
# unique loci, and how many unique locus-target_trait pairs
