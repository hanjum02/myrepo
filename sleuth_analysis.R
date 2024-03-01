library(sleuth)
stab = read.table("sample_info.txt",header=TRUE)
so = sleuth_prep(stab)
#fit a model comparing the two conditions 
so = sleuth_fit(so, ~condition, 'full')
#fit the reduced model to compare in the likelihood ratio test 
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')

library(dplyr)
#extract the test results from the sleuth object 
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval)
write.table(sleuth_significant, file="fdr05_results.txt",quote = FALSE,row.names = FALSE)



