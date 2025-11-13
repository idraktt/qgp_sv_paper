library(dplyr)

merged<-merge(mapping, najeeb.merged, by.x="PHENOTYPIC_ID", by.y="DummyID")

write.csv(merged, file = "phenotype_merged_ranking_mapping",row.names=FALSE,quote = FALSE)

demographics <- read.delim("~/QF-QBB-RES-ACC-0032/scripts/amy/demographics.tsv", stringsAsFactors=FALSE)
questionnaire <- read.delim("~/QF-QBB-RES-ACC-0032/scripts/amy/questionnaire.tsv", stringsAsFactors=FALSE)
BMI <- read.csv("~/QF-QBB-RES-ACC-0032/Phenotypes_6218/BMI.pheno")
`HBA.1C.%` <- read.csv("~/QF-QBB-RES-ACC-0032/Phenotypes_6218/HBA 1C %.pheno", stringsAsFactors=FALSE)
mapping <- read.delim("~/QF-QBB-RES-ACC-0032/mapping/mapping.tsv")
pop_structure <- read.delim("~/QF-QBB-RES-ACC-0032/pop_structure.tsv", header=FALSE)
vw_pivot_measurement <- read.delim("~/QF-QBB-RES-ACC-0032/scripts/amy/amy_subpop/vw_pivot_measurement.tsv")
vw_pivot_bioimpedance <- read.delim("~/QF-QBB-RES-ACC-0032/scripts/amy/vw_pivot_bioimpedance.tsv", stringsAsFactors=FALSE)



mashael<-merge(mapping,BMI,by.x = "PHENOTYPIC_ID",by.y = "DummyID")

mashael<-merge(mashael,`HBA.1C.%`,by.x = "PHENOTYPIC_ID",by.y = "DummyID")

mashael<-merge(mashael,questionnaire,by.x = "QBB_ID",by.y = "PATIENT_IDENTIFIER")

mashael<-merge(mashael,demographics,by.x = "QBB_ID",by.y = "PATIENT_IDENTIFIER")

mashael<-merge(mashael,vw_pivot_measurement,by.x = "QBB_ID",by.y = "PATIENT_IDENTIFIER")

mashael<-merge(mashael,vw_pivot_bioimpedance,by.x = "QBB_ID",by.y = "PATIENT_IDENTIFIER")

mashael<-merge(mashael,pop_structure,by.x = "SAMPLE_NAME",by.y = "V1")


mashael$obesity_class<-"NA"
mashael$surgery_status<-"None"
mashael$bariatric_status<-"FALSE"
mashael$gastric_resection<-"FALSE"
mashael$batch<-"1"



mashael[which(mashael$GENDER_NAME == "Male"),]$surgery_status<-mashael[which(mashael$GENDER_NAME == "Male"),]$NQ_A11_M
mashael[which(mashael$GENDER_NAME == "Female"),]$surgery_status<-mashael[which(mashael$GENDER_NAME == "Female"),]$NQ_A11_F

mashael[which(mashael$BMI <= 18.5),]$obesity_class<-"Underweight"
mashael[which(mashael$BMI > 18.5 & mashael$BMI <= 24.9),]$obesity_class<-"Normal"
mashael[which(mashael$BMI > 24.9 & mashael$BMI <= 29.9),]$obesity_class<-"Overweight"
mashael[which(mashael$BMI > 29.9 & mashael$BMI <= 34.9),]$obesity_class<-"Obesity I"
mashael[which(mashael$BMI > 34.9 & mashael$BMI <= 39.9),]$obesity_class<-"Obesity II"
mashael[which(mashael$BMI > 39.9),]$obesity_class<-"Obesity III"



mashael[grep("Bariatric", mashael$surgery_status), ]$bariatric_status<-"TRUE"
mashael[grep("Gastric resection", mashael$surgery_status), ]$gastric_resection<-"TRUE"
mashael[grep("SI", mashael$NAME), ]$batch<-"2"

surgeries_1st_batch<-table(mashael[which(mashael$batch ==  "1"),]$surgery_status)


output<-
  select(
    mashael,
    QBB_ID,
    PHENOTYPIC_ID,
    BAM_LOC,
    NAME,
    SAMPLE_NAME,
    AGE,
    GENDER_NAME,
    HW_OUT_STANDING_HEIGHT,
    HW_OUT_WEIGHT,
    BMI,
    BMI_Raw_without_outlier,
    obesity_class,
    surgery_status,
    bariatric_status,
    gastric_resection
)



write.table(output,"/home/ealiyev_qgp/QF-QGP-RES-PUB-003/MA/mashael.tsv",sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
write.table(surgeries,"/home/ealiyev_qgp/QF-QGP-RES-PUB-003/MA/surgeries.tsv",sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
write.table( surgeries_1st_batch,"/home/ealiyev_qgp/QF-QGP-RES-PUB-003/MA/surgeries_1st_batch.tsv",sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)




amy_test<-merge(amy_6000_counts,mashael,by.x = "X1",by.y = "NAME")

amy_test$BMI<-as.numeric(amy_test$BMI)
amy_test$AMY1_round<-as.numeric(amy_test$AMY1_round)
amy_test$AMY2A_round<-as.numeric(amy_test$AMY2A_round)
amy_test[which(amy_test$AMY1_round>14),]$AMY1_round<-14
amy_test$AMY2A_round[which(amy_test$AMY2A_round=="0")]<-NA

outlierMax = mean(amy_test$BMI,na.rm = TRUE) + 3*sd(amy_test$BMI,na.rm = TRUE)
outlierMin = mean(amy_test$BMI,na.rm = TRUE) - 3*sd(amy_test$BMI,na.rm = TRUE)

amy_test$BMI[amy_test$BMI> outlierMax] = NA
amy_test$BMI[amy_test$BMI< outlierMin] = NA




summary(lm(BMI ~ AMY2A_round+AGE+factor(GENDER_NAME),
           data=amy_test[which( amy_test$batch == "2" & (amy_test$V2 == "Blueberry" | amy_test$V2 == "Coral")),]))

boxplot(amy_test$BMI_CALC * 10000)

summary(lm(BMI ~ AMY1_round+AGE+factor(GENDER_NAME),
           data=amy_test[which(amy_test$batch == "2"),]))

amy_test$BMI_CALC= as.numeric(as.character(amy_test$HW_OUT_WEIGHT)) / (as.numeric(as.character(amy_test$HW_OUT_STANDING_HEIGHT))^2)


table(amy_test$NQ_A10_F)

table(amy_test$BMI_CALC)

table(amy_test[which(amy_test$HBA.1C.. > 6.5),]$batch)







