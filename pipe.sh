
# Deal with gnarly relative file paths that are inconsistent with the directory we're running them from
cat output/snp-lists/lncrnas-tier1_all-gwas_gtex-ir_source-pval1e-05_lookup-pval1e-05_source-window500000_lookup-window10000_coloc-tests.txt | sed s/munge\\/munged/data\\/gwas/g | sed s/\\/users\\/mgloud\\/projects\\/insulin_resistance\\///g | sed s/\\/users\\/mgloud\\/projects\\/brain_gwas\\///g > output/snp-lists/lncrnas-tier1_all-gwas_gtex-ir_source-pval1e-05_lookup-pval1e-05_source-window500000_lookup-window10000_coloc-tests-modified-paths.txt
