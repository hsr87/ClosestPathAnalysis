##################################################
###### This script calculates AUROC for inferring drug-disease associations using closest path analysis

###### This code uses 'ROCR' library to calculate AUROC values ######
#install.packages("ROCR")
library("ROCR")

###### AUC function ######
SingleDrugAUC <- function(ScoreMatrix, AnswerSet){
 FN = names(which(ScoreMatrix[AnswerSet] < 0))
 ScoreMatsorted = sort(ScoreMatrix,decreasing=T)
 AnsIndex = ScoreMatsorted
 AnsIndex[names(AnsIndex) %in% AnswerSet] = 1
 AnsIndex[!(names(AnsIndex) %in% AnswerSet)] = 0

 pred = prediction(as.numeric(ScoreMatsorted),as.numeric(AnsIndex))
 perf.tmp = performance(pred,"auc")

 perf = performance(pred,"tpr","fpr")
 AUC = as.numeric(perf.tmp@y.values)    
 return(AUC)
}

###### Known drug-target, drug-disease, and gene-disease associations from CTD ######
setwd("D:/INA/ClosestPathAnalysis/rcode") # the directory where your input files exist
dr_tg = read.table("../data/input_file/CTD_chem_gene_ixns_homo_sapiens.tsv",sep="\t",quote="",comment.char="",header=T)
dr_ds = read.table("../data/input_file/CTD_chemicals_diseases_direct.tsv",sep="\t",quote="",comment.char="",header=T)
dg_ds = read.table("../data/input_file/CTD_genes_diseases_direct.tsv",sep="\t",quote="",comment.char="",header=T)

###### Diseases (need to be changed) ######
ina_ds_lists = list.files("../data/input_file/used_drug_disease_pair")
ina_ds_names = tolower(sapply(strsplit(ina_ds_lists,"_images",fixed=T),"[",1))
temp_dg_names = tolower(dg_ds$DiseaseName); temp_dg_names = gsub(",","_",temp_dg_names,fixed=T); temp_dg_names = gsub(" ","_",temp_dg_names,fixed=T);
ina_ds_ids = vector();
ina_ds_check = vector();

for(i in 1:length(ina_ds_names)){
 temp_dg_name= as.vector(dg_ds$Disease_MeSHorOMIM_ID[grepl(ina_ds_names[i] , temp_dg_names, fixed=T)])
 ina_ds_check = c(ina_ds_check, length(temp_dg_name))
 ina_ds_ids = c(ina_ds_ids,paste(temp_dg_name,collapse="&&"))
}

#ina_ds_ids = dg_ds$Disease_MeSHorOMIM_ID[match(ina_ds_names,temp_dg_names)]


####### AUC #######
cla_tg_ds_files = list.files("../data/result_file/tr_ds_dist/")
cla_ds_ont = sapply(strsplit(cla_tg_ds_files,"_",fixed=T),"[",2)
cla_ds_ids = sapply(strsplit(cla_tg_ds_files,"_",fixed=T),"[",3)
cla_ds_ids = paste(cla_ds_ont,substr(cla_ds_ids,1,nchar(cla_ds_ids)-4),sep=":")

inter_ds_ids = 

Scores_all = vector(); AUC = vector(); DSname = vector();
for(dsnum in 1:length(cla_ds_ids))


setwd("../data/result_file/target_score")
Scores_all = vector();
filelists = list.files()
ds_ids = sapply(strsplit(filelists,"_",fixed=T),"[",2)
ds_ids = substr(ds_ids,1,nchar(ds_ids)-4)
filelists = filelists[ds_ids %in% DS_list_new]
ds_ids = ds_ids[ds_ids %in% DS_list_new]

AUC = vector()
DSname = vector()
for(dsnum in 1:length(filelists)){
 temp_score = read.table(filelists[dsnum],sep="\t",comment.char="",quote="",header=T)
 temp_result = data.frame(Scores=temp_score[match(dr_tg[,3],temp_score[,1]),2],drugNames=dr_tg[,1],targetIndex=rep(1,nrow(dr_tg)))
 drugScores = aggregate(. ~ drugNames, data = temp_result, FUN = sum)
 Scores = drugScores[,2] / drugScores[,3]
 Scores = 1 / Scores
 names(Scores) = drugScores[,1]
 dr_set = unique(as.vector(dr_ds[,1]))
 dr_ds_pos = unique(as.vector(dr_ds[dr_ds$Disease_MeSHorOMIM_ID == ds_ids[dsnum],1]))
 dr_ds_neg_whole = dr_set[!(dr_set %in% dr_ds_pos)]
 pos_num = length(dr_ds_pos)
 for(iter in 1:10){
  dr_ds_neg = dr_ds_neg_whole[sample.int(length(dr_ds_neg_whole),pos_num)]
  total_set = c(dr_ds_pos,dr_ds_neg)
  Scores_samp = Scores[names(Scores) %in% total_set]
  AUC_samp = SingleDrugAUC(Scores_samp, dr_ds_pos)
  AUC = c(AUC,AUC_samp)
  DSname = c(DSname,ds_ids[dsnum])
 }
}
 
results = data.frame(AUC,Disease_MeSHorOMIM_ID=DSname,DiseaseName=dr_ds[match(DSname,dr_ds[,4]),3])
results_mean = aggregate(. ~ DiseaseName, data=results, FUN=mean)

setwd("/Users/hsbr/Desktop/Research/INA/data/results/AUC")
write.table(results, "results_seven_diseases.txt",sep="\t", quote=F, col.names=T, row.names=F)
write.table(results_mean[,1:2], "results_mean_seven_diseases.txt",sep="\t", quote=F, col.names=T, row.names=F)

