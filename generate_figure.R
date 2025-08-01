library(ggplot2)
library(MASS)
library(dplyr)
library(nnet)
library(glmnet)
library(twang)
library(pheatmap)
library(xtable)
require("ggrepel")

library(gridExtra)


# drug effect
betas=c(-2.513784,-2.513784,-2.50111,-2.528043,-2.522842,
        -1.174818,-1.144003,
        -0.5319548,-0.5141345,
        0.4884806,0.4888055,0.5057889)
		
save.path = ("~/ZWang-Lab")

pairwise_diff_unique <- function(a) {
  # Compute differences for all unique pairs
  differences <- combn(a, 2, FUN = function(x) x[1] - x[2])
  return(differences)
}

pairwise_diff_unique_name <- function(a) {
  # Compute differences for all unique pairs
  differences <- combn(a, 2, FUN = function(x) paste0(x[1], "v", x[2]))
  return(differences)
}

HR = pairwise_diff_unique_name(c(1:12))
true_CHR=pairwise_diff_unique(betas)

vs = c("1v1", "1v1", "1v1", "1v1", "1v2", "1v2", "1v3", "1v3", "1v4", "1v4", "1v4", 
       "1v1", "1v1", "1v1", "1v2", "1v2", "1v3", "1v3", "1v4", "1v4", "1v4",
       "1v1", "1v1", "1v2", "1v2", "1v3", "1v3", "1v4", "1v4", "1v4", 
       "1v1", "1v2", "1v2", "1v3", "1v3", "1v4", "1v4", "1v4", 
       "1v2", "1v2", "1v3", "1v3", "1v4", "1v4", "1v4", 
       "2v2", "2v3", "2v3", "2v4", "2v4", "2v4", 
       "2v3", "2v3", "2v4", "2v4", "2v4", 
       "3v3", "3v4", "3v4", "3v4", "3v4", "3v4", "3v4",
       "4v4", "4v4", "4v4")
	   
evaluation = data.frame(lambda=NA,
                        BIAS=NA,
                        BIAS_rel=NA,
                        MSE=NA,
                        Scenario=NA,
                        HR=NA,
                        vs=NA,
                        Prevalence=NA)
						
for(k_index in c(4,8,15)){
  for(s_index in c("a","b","c","d","e","f","g","h","i","j")){
    file.names = list.files(path=file.path(save.path,paste0("result_",k_index,s_index)), pattern="fit")

    result = list()
    for(s in 1:5){
      df = read.table(file.path(save.path,paste0("result_",k_index,s_index),file.names[s]),header = F)[,-c(1,2)]
      this.name = gsub("fit|.txt","",file.names[s])
      final_df=t(apply(df,1,pairwise_diff_unique))

      tmp = data.frame(lambda = this.name,
                       BIAS = colMeans(sapply(1:66, function(j) final_df[,j]-true_CHR[j])), 
                       BIAS_rel = colMeans(sapply(1:66, function(j) final_df[,j]-true_CHR[j]))/true_CHR, 
                       MSE = colMeans(sapply(1:66, function(j) (final_df[,j]-true_CHR[j])^2)),
                       Scenario=s_index,
                       HR=HR,
                       vs=vs,
                       Prevalence=paste0(k_index,"%"))
      evaluation = rbind(evaluation, tmp)
      print(paste0(k_index,s_index," - ",this.name,": ",nrow(df)))
    }
  }
}

evaluation = evaluation[-1,]
evaluation$lambda = factor(evaluation$lambda, 
                           levels = c("NoReg", 
                                      "SingleSelNoReg", "DoubleSelNoReg",
                                      "SingleSelLap","DoubleSelLap"),
                           labels = c("SLR", 
                                      "PSS", "PDS",
                                      "PSS-Lap","PDS-Lap"))

evaluation$Prevalence = factor(evaluation$Prevalence,
                         levels = c("4%","8%","15%"))

evaluation$Scenario = toupper(evaluation$Scenario)

evaluation$group_comp = ifelse(evaluation$vs%in%c("1v1","2v2","3v3","4v4"),
                               "Within group","Across group")

p_MSE = ggplot(evaluation, 
       aes(x=Scenario,y=MSE, fill = lambda)) +
  stat_boxplot(geom = "errorbar", width = 0.75) +
  geom_boxplot(outliers = F) +
  theme_classic()+
  theme(legend.position = "bottom",
        legend.box.spacing = unit(0, 'cm'),
        legend.text = element_text(size=10,color="black"),
        legend.title = element_blank(),
        axis.text = element_text(size=10,face="bold",color="black"),
        axis.title = element_text(size=14,color="black"),
        strip.text = element_text(color="black",family = "sans"))+
  ylab("MSE")+
  facet_grid(rows=vars(Prevalence), 
             scales="free",space="free_x",
             axes="all_x",
             labeller = labeller(Prevalence = label_both))+
  ggtitle(NULL)+
  scale_fill_manual(values = c("#FFF08C","#F28CA4","#A1B1EB", "#e6194b", "#4363d8"))

p_Bias = ggplot(evaluation, aes(x=Scenario,y=abs(BIAS),
                                fill = lambda)) +
  stat_boxplot(geom = "errorbar", width = 0.75) +
  geom_boxplot(outliers = F) +
  theme_classic()+
  theme(legend.position = "bottom",
        legend.box.spacing = unit(0, 'cm'),
        legend.text = element_text(size=10,color="black"),
        legend.title = element_blank(),
        axis.text = element_text(size=10,face="bold",color="black"),
        axis.title = element_text(size=14,color="black"),
        strip.text = element_text(color="black",family = "sans"))+
  ylab("Bias")+
  scale_y_log10()+
  facet_grid(rows=vars(Prevalence), 
             scales="free",space="free_x",
             axes="all_x",
             labeller = labeller(Prevalence = label_both))+
  ggtitle(NULL)+
  scale_fill_manual(values = c("#FFF08C","#F28CA4","#A1B1EB", "#e6194b", "#4363d8"))


evaluation$group_comp = ifelse(evaluation$group_comp=="Within group","Within drug class","Between drug class")
evaluation$group_comp = factor(evaluation$group_comp,
                               levels=c("Within drug class","Between drug class"))
							   
p_BIAS_group = ggplot(evaluation, aes(x=Scenario,y=abs(BIAS),
                                      fill = lambda)) +
  stat_boxplot(geom = "errorbar", width = 0.75) +
  geom_boxplot(outliers = F) +
  theme_classic()+
  theme(legend.position = "bottom",
        legend.box.spacing = unit(0, 'cm'),
        legend.text = element_text(size=10,color="black"),
        legend.title = element_blank(),
        axis.text = element_text(size=10,face="bold",color="black"),
        axis.title = element_text(size=14,color="black"),
        strip.text = element_text(color="black",family = "sans"))+
  ylab("Bias")+
  scale_y_log10()+
  facet_grid(rows=vars(Prevalence), 
             cols=vars(group_comp),
             scales="free",space="free_x",
             axes="all_x",
             labeller = labeller(Prevalence = label_both))+
  ggtitle(NULL)+
  scale_fill_manual(values = c("#FFF08C","#F28CA4","#A1B1EB", "#e6194b", "#4363d8"))
