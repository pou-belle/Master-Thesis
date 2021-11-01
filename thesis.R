#Characterization of immune infiltrate in ealry breast cancer based on a multiplex imaging approach 
#The following script is an example of the steps that I followed. 
#All the steps were applied to all the markers in the same way as shown below.


library(tidyverse)
library(sqldf)
library(PairedData)
library(ggplot2)
library(ggpubr)
library(gridExtra)

#Confidence 100% and columns from the raw data
columns_<- select(`SBG_2004_TMA56_Core[1,11,G]_[18377,61546]_cell_seg_data`, 3:5,174,103,78,83,93,88,47,108)
first_file<- columns_[c(`SBG_2004_TMA56_Core[1,11,G]_[18377,61546]_cell_seg_data`$Confidence=="100.00%" ),]
write.table(first_file, "/*/SBG_2004_TMA56_Core[1,11,G]_[18377,61546]_cell_seg_data.txt", sep="\t")

#Split cell_seg to PROMIX and SBG 
cell_seg <-Promix_SBG_merged_cell_seg_data
x <- sbg_all_proteins
sbg_all_proteins %>%
group_by(Phenotype)

for(Phenotype in unique(sbg_all_proteins$Phenotype)){
  write.csv(x[sbg_all_proteins$Phenotype==Phenotype,], file = paste0("newfile_", Phenotype, ".csv"))
}


#Plots of the values for the values/ example for CD4
plot(CD4$Cytoplasm.Opal.520.Mean)


#Thresolding from the minimal value at 5% of the range was selected/ Example of CD4 
CD4_cutoff<- min(CD4$Cytoplasm.Opal.520.Mean)+(max(CD4$Cytoplasm.Opal.520.Mean)-min(CD4$Cytoplasm.Opal.520.Mean))*0.05


#Applying threshold in order to found the cell expression/ Example of CD4 
col_sbg$Expression_520<-ifelse((col_sbg$Cytoplasm.Opal.520.Mean>= 2.88665),1,0)


#Remove Blank surface
sbg_seg<-test[!test$Tissue.Category == "Blank", ]


#Merge with the second raw db for the px 
columns<-sbg_tissue_seg_summary%>%select(2:3,10)
merged <- merge(columns, sbg_expression, by=c("Sample.Name","Tissue.Category"))


#Convert px to mm
merged_summary_sbg_exp$um=(merged_summary_sbg_exp$Region.Area..pixels.)*0.246016


#um to mm 
merged_summary_sbg_exp$mm<-(merged_summary_sbg_exp$um) * 0.000001


#Merging  with Patients_id
merged <- merge(merged_summary_sbg_exp, patients_id, by="Sample.Name",all.x=TRUE)


#Spliting for tumor and stroma for each patient/each patient one-line of information 

test <- function(x) { 
  dat <- read.table(x,sep = "\t",header = TRUE)
  reshaping1 <- reshape(dat, idvar=c("PATIENT.ID"), timevar = "Tissue.Category", direction="wide")
  #reshaping1 <- reshaping1[-c(2,3,5,6)]

  return (reshaping1)

}

filenames <- dir(getwd(), pattern =".txt")
for(i in 1:length(filenames)) {
 splited <- test(filenames[i])
}


#Create the condition for expressions and count the cells and divide them with the surface 

#TOTAL CD4,CD8,CD68,FOXP3
#Example of CD4 
a=unique(aggregate(TRUE_time_point$Expression_520==1,by=list(PATIENT.ID=TRUE_time_point$PATIENT.ID,Tissue.Category=TRUE_time_point$Tissue.Category,mm=TRUE_time_point$mm),FUN=sum))
write.table(a, "/dir/dir/cd4.txt", sep="\t")


cd4_st<-cd4_totals %>%
  group_by(PATIENT.ID) %>%
  mutate(ratio_cd4_totals = x /mm) %>%
  ungroup()
write.table(cd4_st, "/dir/dir/cd4_st.txt", sep="\t")


#Co-expressions of: 

#Example of PDL1-CD4
a=unique(aggregate(TRUE_time_point$Expression_520==1 & TRUE_time_point$Expression_540==1,by=list(PATIENT.ID=TRUE_time_point$PATIENT.ID,Tissue.Category=TRUE_time_point$Tissue.Category,mm=TRUE_time_point$mm),FUN=sum))
write.table(a, "/dir/dir//pdl1_cd4.txt", sep="\t")

pdl1_cd4<-pdl1_cd4_t %>% 
  group_by(PATIENT.ID) %>% 
  mutate(ratio_pdl1_cd4 = x/ mm) %>% 
  ungroup()

write.table(pdl1_cd4, "/dir/dir/pdl1_cd4_ratios.txt", sep="\t")



#PDL1-CK(different approach)
#I created a seperate dataframe for PDL1_CK// only in Tumor 
#I want one patient-one line 

sbg_seg<-pdl1_ck_dataset[!pdl1_ck_dataset$Tissue.Category == "Stroma", ]
write.table(sbg_seg, "/dir/dir/pdl1_ck_dataset_.txt", sep="\t")


cell_seg <-Promix_SBG_merged_cell_seg_data
x <- pdl1_ck_tumor
pdl1_ck_tumor %>%
group_by(PATIENT.ID)

for(PATIENT.ID in unique(pdl1_ck_tumor$PATIENT.ID)){
  write.csv(x[pdl1_ck_tumor$PATIENT.ID==PATIENT.ID,], file = paste0("newfile_", PATIENT.ID, ".csv"))
}


# Example of PD1-CD4 
a=unique(aggregate(TRUE_time_point$Expression_520==1 & TRUE_time_point$Expression_650==1,by=list(PATIENT.ID=TRUE_time_point$PATIENT.ID,Tissue.Category=TRUE_time_point$Tissue.Category,mm=TRUE_time_point$mm),FUN=sum))

write.table(a, "/dir/dir/pd1_cd4.txt", sep="\t")

#CD4 ST
pd1_cd4_st<-pd1_cd4_t%>% 
  group_by(PATIENT.ID) %>% 
  mutate(ratio_pd1_cd4 = x / mm) %>% 
  ungroup()

write.table(pd1_cd4_st, "/dir/dir/pd1_cd4_ratios.txt", sep="\t")


# Example of FOXP3-CD4 
a=unique(aggregate(TRUE_time_point$Expression_520==1 &  TRUE_time_point$Expression_nucleus==1,by=list(PATIENT.ID=TRUE_time_point$PATIENT.ID,Tissue.Category=TRUE_time_point$Tissue.Category,mm=TRUE_time_point$mm),FUN=sum))

write.table(a, "/dir/dir/foxp3_cd4.txt", sep="\t")

#CD4 ST
foxp3_cd4_st<- foxp3_cd4_t %>% 
  group_by(PATIENT.ID) %>% 
  mutate(ratio_foxp3_cd4 = x / mm) %>% 
  ungroup()

write.table(foxp3_cd4_st, "/dir/dir/foxp3_cd4_ratios.txt", sep="\t")



#Apply medians to find H and L expression to totals and pdl1-ck 


dataframe_sbg$cd4_tumor_ex<-ifelse(is.na(dataframe_sbg$ratio_cd4_totals.Tumor), "NA",
                             ifelse(dataframe_sbg$ratio_cd4_totals.Tumor >270.9524407 ,"H","L"))

dataframe_sbg$cd4_stroma_ex<-ifelse(is.na(dataframe_sbg$ratio_cd4_totals.Stroma), "NA",
                              ifelse((dataframe_sbg$ratio_cd4_totals.Stroma >136.5281),"H","L"))

dataframe_sbg$cd8_tumor_ex<-ifelse(is.na(dataframe_sbg$ratio_cd8_totals.Tumor), "NA",
                            ifelse((dataframe_sbg$ratio_cd8_totals.Tumor >23.87186051 ),"H","L"))

dataframe_sbg$cd8_stroma_ex<-ifelse(is.na(dataframe_sbg$ratio_cd8_totals.Stroma), "NA",
                              ifelse((dataframe_sbg$ratio_cd8_totals.Stroma >68.28577818),"H","L"))

dataframe_sbg$cd68_tumor_ex<-ifelse(is.na(dataframe_sbg$ratio_cd68_totals.Tumor), "NA",
                            ifelse((dataframe_sbg$ratio_cd68_totals.Tumor >59.3201796 ),"H","L"))

dataframe_sbg$cd68_stroma_ex<-ifelse(is.na(dataframe_sbg$ratio_cd68_totals.Stroma), "NA",
                                    ifelse((dataframe_sbg$ratio_cd68_totals.Stroma >142.3092889 ),"H","L"))

pdl1_ck_density_$pdl1_ck_ex<-ifelse(pdl1_ck_density_$`PDL1-CK density`>6.97596428 ,"H","L")

#apply >0 , positive negative
dataframe_sbg$foxp3_tumor_ex<-ifelse(is.na(dataframe_sbg$ratio_foxp3_totals.Tumor), "NA",
                            ifelse((dataframe_sbg$ratio_foxp3_totals.Tumor >0 ),"P","N"))
dataframe_sbg$foxp3_stroma_ex<-ifelse(is.na(dataframe_sbg$ratio_foxp3_totals.Stroma), "NA",
                                ifelse((dataframe_sbg$ratio_foxp3_totals.Stroma >0 ),"P","N"))


#PDL1, PD1 > 0 , p, n 
dataframe_sbg$pdl1_cd4_tumor_ex<-ifelse(is.na(dataframe_sbg$ratio_foxp3_totals.Tumor), "NA",
  ifelse((dataframe_sbg$ratio_pdl1_cd4.Tumor >0 ),"P","N"))




#Baseline Characteristics for each condition and each compartment, create contigency tables 

#ER
cd4_totals_tumor<-table(all_included_sbg$`er positive`, all_included_sbg$cd4_tumor_ex)
cd4_totals_stroma<-table(all_included_sbg$`er positive`, all_included_sbg$cd4_stroma_ex)


#Menopausal 

cd4_totals_tumor<-table(all_included_sbg$`menopausal status` , all_included_sbg$cd4_tumor_ex)
cd4_totals_stroma<-table(all_included_sbg$`menopausal status` , all_included_sbg$cd4_stroma_ex)



#Treatarm
cd4_totals_tumor<-table(all_included_sbg$`treatarm` , all_included_sbg$cd4_tumor_ex)
cd4_totals_stroma<-table(all_included_sbg$`treatarm` , all_included_sbg$cd4_stroma_ex)



#tumor_size
cd4_totals_tumor<-table(all_included_sbg$`tumour size`>=20 , all_included_sbg$cd4_tumor_ex)
cd4_totals_stroma<-table(all_included_sbg$`tumour size`>=20 , all_included_sbg$cd4_stroma_ex)




#Grade
cd4_totals_tumor<-table(all_included_sbg$`elston grade` , all_included_sbg$cd4_tumor_ex)
cd4_totals_stroma<-table(all_included_sbg$`elston grade` , all_included_sbg$cd4_stroma_ex)



# PGR Positive
cd4_totals_tumor<-table(all_included_sbg$`pgr positive` , all_included_sbg$cd4_tumor_ex)
cd4_totals_stroma<-table(all_included_sbg$`pgr positive` , all_included_sbg$cd4_stroma_ex)



#ER+ PGR 
cd4_totals_tumor<-table(all_included_sbg$`er and pgr negative` , all_included_sbg$cd4_tumor_ex)
cd4_totals_stroma<-table(all_included_sbg$`er and pgr negative` , all_included_sbg$cd4_stroma_ex)




#HER2 Overexxpression 
cd4_totals_tumor<-table(all_included_sbg$`her2 overexpression` , all_included_sbg$cd4_tumor_ex)
cd4_totals_stroma<-table(all_included_sbg$`her2 overexpression` , all_included_sbg$cd4_stroma_ex)




#PDL1-CK
pdl1_ck_1<-table(pdl1_ck_d$`elston grade` , pdl1_ck_d$pdl1_ck_ex)
write.table(pdl1_ck, "/dir/dir//pdl1_ck_grade.txt", sep="\t")

pdl1_ck_2<-table(pdl1_ck_d$`pgr positive` , pdl1_ck_d$pdl1_ck_ex)
write.table(pdl1_ck, "/dir/dir//pdl1_ck_pgr_positive.txt", sep="\t")

pdl1_ck_3<-table(pdl1_ck_d$`er and pgr negative` , pdl1_ck_d$pdl1_ck_ex)
write.table(pdl1_ck, "/dir/dir//pdl1_ck_er_pgr_negative.txt", sep="\t")

pdl1_ck_4<-table(pdl1_ck_d$`her2 overexpression` , pdl1_ck_d$pdl1_ck_ex)
write.table(pdl1_ck, "/dir/dir//pdl1_ck_her2_overexpressed.txt", sep="\t")

pdl1_ck_5<-table(pdl1_ck_d$`tumour size`>=20 , pdl1_ck_d$pdl1_ck_ex)
write.table(pdl1_ck, "/dir/dir//pdl1_ck_tumor_size.txt", sep="\t")

pdl1_ck_6<-table(pdl1_ck_d$`menopausal status`, pdl1_ck_d$pdl1_ck_ex)
write.table(pdl1_ck, "/dir/dir//pdl1_ck_menop.txt", sep="\t")

pdl1_ck_7<-table(pdl1_ck_d$`er positive`, pdl1_ck_d$pdl1_ck_ex)
write.table(pdl1_ck, "/dir/dir//pdl1_ck_ER.txt", sep="\t")


#Removing uncertain category from baseline characteristics 
cd4_totals_tumor<-table(pairwised_data$`menopausal status` , pairwised_data$cd4_tumor_ex)
cd4_totals_tumor<-cd4_totals_tumor[c(1,2,3), ]
cd4_totals_tumor<-table(pairwised_data$`her2 overexpression` , pairwised_data$cd4_tumor_ex)
cd4_totals_tumor<-cd4_totals_tumor[c(1,2), ]
cd4_totals_tumor_her2_overexpressed<- cd4_totals_tumor_her2_overexpressed[ -c(3) ]




#Statistical tests


#P-adjust
p.adjust(pvalues,method="fdr") 


#Wilcoxon 
wilcox.test(x= pairwised_data$`CD4 Single T`, y=pairwised_data$`CD4 Single S`, paired = TRUE, alternative = "two.sided")

#Wilcoxon plots
#Example for CD4,CD8, CD68, FOXP3
CD4.SINGLE.TUMOR<-pairwised_data$`CD4 Single T`
CD4.SINGLE.STROMA<- pairwised_data$`CD4 Single S`
d<-data.frame(Tumor= CD4.SINGLE.TUMOR,Stroma = CD4.SINGLE.STROMA)
p1<- ggpaired(d,cond1 = "Tumor", cond2 = "Stroma", line.color = "gray", line.size = 0.4,
fill = "condition" ,palette = "npg", title = "CD4 Single cells",subtitle = "p-value = 4.520000e-01",
            ylab = "Cell Density(positive cells/mm2)",
            xlab = "Patients N=79",
            legend.title = "Compartments") 

CD8.SINGLE.TUMOR<- pairwised_data$`CD8 Single T`
CD8.SINGLE.STROMA<- pairwised_data$`CD8 Single S`
d1<-data.frame(Tumor= CD8.SINGLE.TUMOR, Stroma=CD8.SINGLE.STROMA)
p2<- ggpaired(d1,cond1 = "Tumor", cond2 = "Stroma", line.color = "gray", line.size = 0.4,
fill = "condition" ,palette = "npg", ylab = "Cell Density(positive cells/mm2)", title = "CD8 Single cells",subtitle = "p-value = 9.856000e-07",
            xlab = "Patients N=79",
            legend.title = "Compartments")

CD68.SINGLE.TUMOR<- pairwised_data$`CD68 Single T`
CD68.SINGLE.STROMA<- pairwised_data$`CD68 Single S`
d2<-data.frame(Tumor=CD68.SINGLE.TUMOR,Stroma=CD68.SINGLE.STROMA)
p3<- ggpaired(d2,cond1 = "Tumor", cond2 = "Stroma", line.color = "gray", line.size = 0.4,
    fill = "condition" ,palette = "npg", ylab = "Cell Density(positive cells/mm2)",title = "CD68 Single cells",subtitle = "p-value = 1.158000e-06 ",
            xlab = "Patients N=79",
            legend.title = "Compartments")

FOXP3.SINGLE.TUMOR<- pairwised_data$`FOXP3 Single T`
FOXP3.SINGLE.STROMA<- pairwised_data$`FOXP3 Single S`
d3<- data.frame(Tumor=FOXP3.SINGLE.TUMOR,Stroma=FOXP3.SINGLE.STROMA)
p4<- ggpaired(d3,cond1 = "Tumor", cond2 = "Stroma", line.color = "gray", line.size = 0.4,
    fill = "condition" ,palette = "npg" ,ylab = "Cell Density(positive cells/mm2)",title = "FOXP3 Single cells",subtitle = "p-value = 1.638400e-08 ",
            xlab = "Patients N=79",
            legend.title = "Compartments")


postscript(file= "singles.eps",
           width = 8.72, height = 11.69,
           horizontal = FALSE, onefile = FALSE,
           paper = "special")
grid.arrange(p1, p2,p3,p4 ,nrow = 2)
dev.off()



#Spearman, pairs from stroma and tumor 
cor.test(x= pairwised_data$`CD4 Single T`, y=pairwised_data$`CD4 Single S`, method = "spearman")


data_rcorr <-as.matrix(pairwised_data[, 2:32])

mat_2 <-rcorr(data_rcorr)

write.table(p_value, "spear.txt", sep="\t")


upper<-spear
upper[upper.tri(spear)]<-""
upper<-as.data.frame(upper)

upper.tri(spear, diag = FALSE)

ggcorrplot(upper, hc.order = TRUE, type = "lower",
   outline.col = "white", order = "hclust",
   ggtheme = ggplot2::theme_gray,
   colors = c("#6D9EC1", "white", "#E46726"))
dev.off()

ggsave("correlation.png")
ggcorrplot(spear, nbreaks = 4,min_size = 0, max_size = 1, hjust = 0.75, size = 2,layout.exp = 10,type = "upper")
dev.off()



#Fisher test 
fisher.test(contigency tables)






























