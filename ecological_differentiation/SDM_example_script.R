#SCRIPT EXAMPLE FOR LESSER KESTREL SDMs


#example with the western clade during the breeding season - "breeding ovest"


library(raster)
library(dismo)
library(ggmap)
library(ENMeval)
library(rgdal)
library(SDMtune)
library(stringr)
library(fasterize)
library(xlsx)


#load the script to compute AICc
source("AICc_BC.R")
#taken from Brambilla et al. 2022
#https://dataverse.unimi.it/file.xhtml?persistentId=doi:10.13130/RD_UNIMI/ARAI8C/SLVDCC&version=1.1


current=stack("current_climate_stack.grd")


swd_breeding=readRDS("SWD_breeding_ovest_large_train_admixed_all.rds")
swd_breeding_test=readRDS("SWD_breeding_ovest_large_test_admixed_all.rds")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
############## select best regularization multiplier for basic model ##########
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


stack_prediction=current

rm_values=seq(0.5,5,0.5)
aic_rm=data.frame(rm=rm_values,AICc=NA,k=NA)

niter=1000



for (j in aic_rm$rm) {
  
  modellobase_breeding_rmj <- train(method = "Maxent", data = swd_breeding, fc = "lq", reg = j, iter = niter, addSamplovestoBg = F)
  nparam_rmj=nrow(modellobase_breeding_rmj@model@coeff)
  aic_i_rmj=AICc_BC(modellobase_breeding_rmj,swd_breeding)
  aic_rm[aic_rm$rm==j,"AICc"]<-aic_i_rmj
  aic_rm[aic_rm$rm==j,"k"]<-nparam_rmj
  
}


# select best regularization mutiplier
rm_ok_base_breeding=aic_rm$rm[which.min(aic_rm$AICc)]

print(rm_ok_base_breeding)

modellobase_breeding_rmsel <- train(method = "Maxent", data = swd_breeding, fc = "lq", reg = rm_ok_base_breeding, iter = niter, addSamplovestoBg = T)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
############## check correlation between variables #####################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

coord_bk=background_breeding[,1:2]
names(coord_bk)

bgcorr <- prepareSWD(species = "Bgs", a = coord_bk, env = current)

plotCor(bgcorr, method = "pearson", cor_th = 0.7)
#all <0.75, no issues

#not run as there are no variables so highly correlated:
#selected_variables_model <- varSel(modellobase_breeding_rmsel, metric = "auc", test = swd_breeding_test, bg4cor = bgcorr, method = "pearson", cor_th = 0.75, permut = 3)

selected_variables_model <- modellobase_breeding_rmsel 


swd_breeding_no_corr=SDMtune:::SWD(species = swd_breeding@species,
                                    coords = swd_breeding@coords,
                                    data = selected_variables_model@data@data,
                                    pa=swd_breeding@pa)


modello_no_var_corr <- train(method = "Maxent", data = swd_breeding_no_corr, fc = "lq", reg = rm_ok_base_breeding, iter = niter,  addSamplovestoBg = T)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
############## EXCLUDING VARIABLES WITH LAMBDA == 0 #####################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

df_variabili_lambda_0=data.frame(variabile=names(modello_no_var_corr@data@data),sel=NA)
df_variabili_lambda_0$variabile=as.character(df_variabili_lambda_0$variabile)

for (i in df_variabili_lambda_0$variabile) {
  sel_i=any(str_detect(modello_no_var_corr@model@coeff$feature,i))
  df_variabili_lambda_0[df_variabili_lambda_0$variabile==i,"sel"]<-sel_i
}

# BUILD A VECTOR WITH VARIABLES TO BE REMOVED
variabili_lambda_0=df_variabili_lambda_0[df_variabili_lambda_0$sel==F,"variabile"]

# PREPARE NEW SWD FOR VARIABLE SELECTION
swd_variable_selection=SDMtune:::SWD(species = swd_breeding_no_corr@species,
                                     coords = swd_breeding_no_corr@coords,
                                     data = swd_breeding_no_corr@data[,!(names(swd_breeding_no_corr@data)%in%variabili_lambda_0)]  ,
                                     pa=swd_breeding_no_corr@pa)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
############## LOOP FOR AICc-BASED MODEL SELECTION #####################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###




repeat {
  
  
  modello_loop_1 <- train(method = "Maxent", data = swd_variable_selection, fc = "lq", reg = rm_ok_base_breeding, iter = niter,  addSamplovestoBg = T)
  AICc_1=AICc_BC(modello_loop_1,swd_variable_selection)
  permutazione1=varImp(modello_loop_1,permut=10)
  variabile_peggiore=permutazione1$Variable[dim(permutazione1)[1]]
  swd_2=SDMtune:::SWD(species = swd_variable_selection@species,
                      coords = swd_variable_selection@coords,
                      data = swd_variable_selection@data[,!(names(swd_variable_selection@data)%in%variabile_peggiore)]  ,
                      pa=swd_variable_selection@pa)
  modello_loop_2 <- train(method = "Maxent", data = swd_2, fc = "lq", reg = rm_ok_base_breeding, iter = niter, addSamplovestoBg = T)
  AICc_2=AICc_BC(modello_loop_2,swd_2)
  delta_aicc=AICc_2-AICc_1
  
  
  
  if   (delta_aicc>0) {
    
    break
    
  }
  
  swd_variable_selection=swd_2
  
}

print(modello_loop_1@model@coeff)




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
############## 8. MODEL TUNING #####################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###







df_tuning=data.frame(rm=rep(rm_values,2),
                     features=rep(c("l","lq"),each=length(rm_values)),
                     AICc=NA,niter=NA)


for (i in c("l","lq")) {
  
  for (j in rm_values) {
    
    modello_ij <- train(method = "Maxent", data = swd_variable_selection, fc = i, reg = j, iter = niter,  addSamplovestoBg = T)
    aicc_ij=AICc_BC(modello_ij,swd_variable_selection)  
    niter_ij=modello_ij@model@results[rownames(modello_ij@model@results)=="Iterations"]
    df_tuning[df_tuning$rm==j&df_tuning$features==i,"AICc"]<-aicc_ij
    df_tuning[df_tuning$rm==j&df_tuning$features==i,"niter"]<-niter_ij
  }
  
}

tuned_reg=df_tuning$rm[which.min(df_tuning$AICc)]
tuned_feat=as.character(df_tuning$features[which.min(df_tuning$AICc)])

tuned_model<- train(method = "Maxent", data = swd_variable_selection, fc = tuned_feat, reg = tuned_reg, iter = niter,  addSamplovestoBg = T)




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
############## 9. MODEL EVELUATION #####################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


predict.test.i=predict(tuned_model,data=swd_breeding_test,type="cloglog")


auc_train_i=auc(tuned_model,test=swd_variable_selection)
auc_test_i=auc(tuned_model,test=swd_breeding_test)
tss_train_i=tss(tuned_model,test=swd_variable_selection)
tss_test_i=tss(tuned_model,test=swd_breeding_test)

mtss_threshold_tuned=tuned_model@model@results[rownames(tuned_model@model@results)=="Maximum.training.sensitivity.plus.specificity.Cloglog.threshold"]
ten_perc_threshold_tuned=tuned_model@model@results[rownames(tuned_model@model@results)=="X10.percentile.training.presence.Cloglog.threshold"]
mtp_threshold_tuned=tuned_model@model@results[rownames(tuned_model@model@results)=="Minimum.training.presence.Cloglog.threshold"]


pred_occ_test=predict.test.i[1:length(swd_breeding_test@pa[swd_breeding_test@pa==1])]
pred_binary_test_mtss=ifelse(pred_occ_test>mtss_threshold_tuned,1,0)
pred_binary_test_10th=ifelse(pred_occ_test>ten_perc_threshold_tuned,1,0)
pred_binary_test_mtp=ifelse(pred_occ_test>mtp_threshold_tuned,1,0)

OR_tuned=data.frame(threshold=c("mtss","mtp","10th_perc"),fa=NA,tp=NA,OR=NA)

OR_tuned[OR_tuned$threshold=="mtss","fa"]<-length(pred_binary_test_mtp[pred_binary_test_mtss==0])
OR_tuned[OR_tuned$threshold=="mtss","tp"]<-length(pred_binary_test_mtp[pred_binary_test_mtss==1])
OR_tuned[OR_tuned$threshold=="mtp","fa"]<-length(pred_binary_test_mtp[pred_binary_test_mtp==0])
OR_tuned[OR_tuned$threshold=="mtp","tp"]<-length(pred_binary_test_mtp[pred_binary_test_mtp==1])
OR_tuned[OR_tuned$threshold=="10th_perc","fa"]<-length(pred_binary_test_mtp[pred_binary_test_10th==0])
OR_tuned[OR_tuned$threshold=="10th_perc","tp"]<-length(pred_binary_test_mtp[pred_binary_test_10th==1])
OR_tuned$OR=OR_tuned$fa/(OR_tuned$fa+OR_tuned$tp)  

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
############## 10. get features #####################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# PLOTS FOR RESPONSES
for (y in names(tuned_model@data@data)) {
  
  plotResponse(tuned_model,var=y,type="cloglog",only_presence = F,marginal = T,fun=mean,rug = T)
  #  ggplot(data=a[1]$data, aes(x=x, y=y)) +
  #  geom_line()+
  #  geom_rug(mapping=aes(x=x),data=data.frame(x=tuned_model@data@data[which(tuned_model@data@pa==1),y],
  #                                            y=tuned_model@data@pa[which(tuned_model@data@pa==1)]*max(a[1]$data[,2])),sides="t")+
  #  geom_rug(mapping=aes(x=x),data=data.frame(x=tuned_model@data@data[which(tuned_model@data@pa==0),y],
  #                                          y=tuned_model@data@pa[which(tuned_model@data@pa==0)]+min(a[1]$data[,2])),sides="b")
  
  ggsave(paste("plot_response_","_","breeding",y,".jpg",sep=""), plot = last_plot(), device = "jpeg",width = 16, height = 12, units = c("cm"), dpi = 200)
}


#COMPUTE PERMUTATION IMPORTANCE
permutazione_tuned=varImp(tuned_model,permut=10)


#MAKE A TABLE (ONE ROW) WITH EVALUATION PARAMETERS
evaluation_df=data.frame(tss_train=tss_train_i,tss_test=tss_test_i,auc_train=auc_train_i,auc_test=auc_test_i,
                         mtss_fa=OR_tuned[OR_tuned$threshold=="mtss","fa"],
                         mtss_tp=OR_tuned[OR_tuned$threshold=="mtss","tp"],
                         mtss_OR=OR_tuned[OR_tuned$threshold=="mtss","OR"],
                         mtp_fa=OR_tuned[OR_tuned$threshold=="mtp","fa"],
                         mtp_tp=OR_tuned[OR_tuned$threshold=="mtp","tp"],
                         mtp_OR=OR_tuned[OR_tuned$threshold=="mtp","OR"],
                         tenth_perc_fa=OR_tuned[OR_tuned$threshold=="10th_perc","fa"],
                         tenth_perc_tp=OR_tuned[OR_tuned$threshold=="10th_perc","tp"],
                         tenth_perc_OR=OR_tuned[OR_tuned$threshold=="10th_perc","OR"])

write.csv(rbind(tuned_model@model@results,tuned_reg),paste("breeding","_results.csv",sep=""))

# WRITE OTHER RESULTS AND PARAMETERS
write.csv(tuned_model@model@coeff,paste("breeding","_coeff.csv",sep=""))
write.csv(df_tuning,paste("breeding","_df_tuning.csv",sep=""))
write.csv(permutazione_tuned,paste("breeding","_perm_imp.csv",sep=""))
write.csv(evaluation_df,paste("breeding","_evaluation.csv",sep=""))


# MAKE RASTER WITH MODEL PREDICTION
prediction_tuned=predict(tuned_model,data=stack_prediction,type="cloglog")
writeRaster(prediction_tuned, paste("breeding","_prediction.tif",sep=""),overwrite=T)



save.image(paste("modello_","_","breeding",".RData",sep=""))



# MAKE PROJECTIONS ####
#for any other scenario:

prediction_scenario=predict(tuned_model,data=scenario_stack,type="cloglog")

plot(prediction_scenario)
writeRaster(prediction_scenario, paste("scenario_breeding_ovest","_prediction.tif",sep=""))
