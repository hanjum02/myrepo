# loads and libraries the sleuth package into R's environment, used for DE analysis of 2 conditions
library(sleuth)

# variable that reads in the data stored in the sample_info.txt file, containing kallisto output and sample descriptions
stab = read.table("sample_info.txt",header=TRUE)

# establishes a sleuth object which will be used to run the DE analysis
so = sleuth_prep(stab)

#fits a model comparing the two conditions, 2dpi to 6dpi
so = sleuth_fit(so, ~condition, 'full')

#fits a reduced model to compare the likelihood ratio test
so = sleuth_fit(so, ~1, 'reduced')

# using the previously formed models, computes likelihood ratio test for DE analysis between 2dpi and 6pi conditions
so = sleuth_lrt(so, 'reduced', 'full')

# loads and libraries the dplyr package into R's environment, used for filtering data.frames
library(dplyr)

# extracts the results from the sleuth object and analysis and stores them in variable 'sleuth_table"'
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

# checks 'sleuth_table' for significant results based on our set threshold, or those with FDR/q-val < .05, then orders the results by pval
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval)

# writes the top results of the DE analysis, entries with FDR < .05, to file 'fdr05_results.txt'
write.table(sleuth_significant, file="fdr05_results.txt",quote = FALSE,row.names = FALSE)
