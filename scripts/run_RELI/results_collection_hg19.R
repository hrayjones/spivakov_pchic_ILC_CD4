library("writexl")
phenotypes<-readLines("/Volumes/miraldiNB/Zi/RELI/ILC3-PC-HiC/hg19/012025/Disease_RELI/run/EUR/phenotypes.list")
dir_out<-'/Volumes/miraldiNB/Zi/RELI/ILC3-PC-HiC/hg19/012025/Disease_RELI'
final_results<-NULL
for (ix in phenotypes){
  file_in<-file.path(dir_out,paste0('run/EUR/',ix,'/LDexpansion/snp_expanded/AllStats.tab'))
  result<-read.delim(file_in,header=TRUE,sep="\t")
  result<-result[,-c(1,12)]
  result$Phenotype<-paste0(ix)
  colnames(result)<-c('ATAC_library','TF','Overlap','Total','Ratio','Mean','Std_Dev','Z_score','Enrichment',
                      'P_val','Null_Model','Species','Phenotype')
  final_results<-rbind(final_results,result)
}
order<-c('Phenotype','ATAC_library','TF','Overlap','Total','Ratio','Mean','Std_Dev','Z_score','Enrichment',
         'P_val','Null_Model','Species')
final_results<-final_results[,order]
final_results<-final_results[order(final_results$P_val,decreasing=FALSE),]
file_out<-file.path(dir_out,paste0('Disease_RELI_results_raw_Pval.xlsx'))
write_xlsx(final_results,file_out)
#filter for phenotypes >=10 independent loci and perform Padj
final_results<-final_results[final_results$Total>=10,]
order_celltype<-unique(final_results$ATAC_library)
Padj_final_results<-NULL
for (ix in order_celltype){
  results_sub<-final_results[final_results$ATAC_library %in% ix,]
  results_sub$Padj<-p.adjust(results_sub$P_val,method='bonferroni',n=nrow(results_sub))
  Padj_final_results<-rbind(Padj_final_results,results_sub)
}
Padj_final_results$`-log10(Padj)`= -log10(Padj_final_results$`Padj`)
order<-c('Phenotype','ATAC_library','TF','Overlap','Total','Ratio','Mean','Std_Dev','Z_score','Enrichment',
         'P_val','Padj','Null_Model','Species','-log10(Padj)')
Padj_final_results<-Padj_final_results[,order]
Padj_final_results<-Padj_final_results[order(Padj_final_results$`-log10(Padj)`,decreasing=TRUE),]
file_out<-file.path(dir_out,paste0('Disease_RELI_results_min_10.xlsx'))
write_xlsx(Padj_final_results,file_out)
file_out<-file.path(dir_out,paste0('Disease_RELI_results_min_10.txt'))
write.table(Padj_final_results,file_out,sep="\t",row.names=FALSE,col.names = TRUE,quote=FALSE)