
library(data.table)
library(xlsx)
library(dplyr)
library(tidyr)
library(reshape2)
library(maftools)
library(ggplot2)
library(ggpubr)
library(stringr)
library(pheatmap)


wd <- "D:/Github-Projects/Pipeline_Delly-SV/"
setwd(wd)

output <- paste0("OUTPUT-SV_Filter_R.",Sys.Date())
dir.create(output)
#                                        ===================================
#                                                ANALISE DOS RESULTADOS
#                                        ===================================


#==============================================================================#
# STEP 0 -  Datasets reorganization ####
#==============================================================================#
all.sv<-fread("Result_Delly_SV.2023-10-28/step6_annovar/ROP-annovar.hg38_multianno.txt", header = T)

all.sv<- fread(file="Result_Delly_SV.2023-10-28/step6_annovar/ROP-annovar.hg38_multianno.txt", header=T, sep="\t", na.strings=c(".","NA"), fill=TRUE, check.names = FALSE)
samples_ID <-fread("Result_Delly_SV.2023-10-28/step6_annovar/sampleID.list", header = F)

# Rename columns
#==============================================================================#
colnames(all.sv) <-  c( colnames(all.sv[,c(1:11)]),"Otherinfo1","Otherinfo2","Otherinfo3","CHR","POS","ID","REF","ALT","QUAL","FILTER", "INFO", "FORMAT", 
                          samples_ID$V1)

# Creating INDEX (CHR,POS,REF,ALT) and select the main columns
#==============================================================================#
all.sv$index <- all.sv$ID
all.sv <-data.frame(all.sv)
exCol=colnames(all.sv[,-which(names(all.sv) %in% c("Otherinfo1","Otherinfo2","Otherinfo3","ID","index"))])
all.sv <- select(all.sv, index,all_of(exCol))
dim(all.sv)   #19142   105 
fwrite(all.sv, paste0(output,"/tab0.tsv"), quote = F, sep="\t")



#==============================================================================#
#                 STEP 1 - Clear dataset ####
#==============================================================================#
# tab1<-as.data.frame(fread(paste0(output,"/tab0.tsv"), header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE))
tab1<- all.sv


## Restructure tables, all samples in the same column ###
#==============================================================================#
tab1 <- reshape2::melt(tab1,id=colnames(tab1[,1:which(names(tab1)== "FORMAT")]))
colnames(tab1)[colnames(tab1) == "variable"] <- "SAMPLE"


unique(tab1$FORMAT) #GT:GL:GQ:FT:RCL:RC:RCR:RDCN:DR:DV:RR:RV

FORMAT<-strsplit(unique(tab1$FORMAT), ":")[[1]]
tab2<- separate(data = tab1, col = value, into =FORMAT, sep = ":")

table(tab2$GT)
# tab2 <- tab2[!grepl("/\\.", tab2$GT),]
tab2 <- tab2[!grepl("^0/0", tab2$GT),]
tab2 <- tab2[!grepl("^0\\|0", tab2$GT),]
tab2 <- tab2[!grepl("^\\./\\.", tab2$GT),]
table(tab2$GT)

#==============================================================================#
### SEPARAR OS ALELOS DE CADA INDIVIDUO
#==============================================================================#
setDT(tab2)[,paste0("GT", 1:2) := tstrsplit(GT, "/")]
#==============================================================================#

tab2$SV <- substr(tab2$index, 1, 3)
table(tab2$SV)
table(tab2$GT)

#==============================================================================#
### AMOSTRAS EXCLUIDAS POR DADOS CLINICOS
#==============================================================================#
tab2$SAMPLE <- gsub("\\.", "-", tab2$SAMPLE)

samples.excl <- c('ROP-116','ROP-118','ROP-120','ROP-121','ROP-122','ROP-123','ROP-125','ROP-126',
                  'ROP-127','ROP-130','ROP-131','ROP-132','ROP-23','ROP-84',
                  'ROP-83','ROP-107','ROP-81' ,'ROP-91' ,'ROP-129')
samples.excl <-paste(samples.excl, collapse = "|" )

tab2 <- tab2[!grepl(samples.excl, tab2$SAMPLE),]
tab.info <- as.data.frame(tab2$INFO)


table(tab2$SAMPLE,tab2$SV)


#==============================================================================#
#                 STEP 2 - FILTERS ####
#==============================================================================#


#==============================================================================#
### SELECIONAR OS PRECISE & PASS
#==============================================================================#
tab3 <- tab2[grepl("^PRECISE", tab2$INFO),]
table(tab3$SV)

tab3 <- subset(tab3, FT == "PASS") 
table(tab3$SV)


# Quantificando mutações compartilhadas >=2
#==============================================================================#
aux <- data.frame(table(tab3$index))
tab3 <- merge(tab3, aux, by.x = "index", by.y = "Var1")
colnames(tab3)[names(tab3)=="Freq"]<-"MUTATION_shared"
shared2exclude <- which((tab3$MUTATION_shared >= 2) )
nrow(unique(tab3[shared2exclude,"index"]))  #3048 
tab3 <- tab3[-shared2exclude,] 

# selecionando evento em regiões codificadoras
#==============================================================================#
data.table(table(tab3$Func.refGene))
# FALTA IMPLEMENTAR.  REMOVER APENAS OS EVENTOS (DEL, INV) QEU OCORREM EM REGIÃO INTRONICA

location2keep_mut.xg <- c(
  which(tab3$Func.refGene=="exonic"),
  which(tab3$Func.refGene=="exonic;splicing"),
  which(tab3$Func.refGene=="splicing"))
tab3 <- tab3[location2keep_mut.xg,]

tab3 <- subset(tab3, Func.refGene %in% c("splicing", "exonic", "exonic;splicing", "splicing;exonic"))



table(tab3$SV)
table(tab3$SAMPLE,tab3$SV)
tab.count <- data.frame(table(tab3$SAMPLE,tab3$SV))
tab.count<- as.data.frame(reshape2::dcast(tab.count, Var1 ~ Var2,  value.var = "Freq"))
tab.count$SVs <- 
  
fwrite(tab3, paste0(output,"/tabela_all.SVs.tsv"), quote = F, sep="\t")
fwrite(tab.count, paste0(output,"/tabela_count.SVs.tsv"), quote = F, sep="\t")




dgenes<-fread("input_R/driver-genes.txt", header = F, col.names = "gene")
dMMR<-fread("input_R/dMMR_genes.txt", header = F, col.names = "gene")

resp.trat <- read.xlsx("Dados Clinicos JP final  - Dez2022 - FILTRADO.xlsx", sheetIndex = 1, header = T)
resp.trat<- resp.trat[rowSums(is.na(resp.trat)) != ncol(resp.trat), -c(6)]
resp.trat <- resp.trat[resp.trat$WXS== "Y",]
resp.trat$Distance.from.anal.verge <- gsub("cm", "", resp.trat$Distance.from.anal.verge)
resp.trat$Tumor.size <- gsub("cm", "", resp.trat$Tumor.size)

all.mut<-fread("output_R/mutation_codRegion.tsv", header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
all.mut.flt<-fread("output_R/mutation_codRegion_Cov30Frq02.tsv", header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
somatic<-fread("output_R/somatic_mutation-ALL.tsv", header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
somatic<-subset(somatic, (sample %in% resp.trat$Sample) )
somatic_Samples<- fread("output_R/MATH_TMB-SAMPLES.tsv", header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
somatic_Samples<-subset(somatic_Samples, (sample %in% resp.trat$Sample) )

somatic_Samples.resp <- merge(somatic_Samples, resp.trat, by.x = "sample", by.y = "Sample")

somatic <- merge(somatic, resp.trat[,c("Sample","ID_Exoma")], by.x = "sample", by.y = "Sample", all.x = T)
somatic<- somatic %>% 
  mutate(Gene.refGene = strsplit(as.character(Gene.refGene), ";")) %>% 
  unnest(Gene.refGene)

somatic<-subset(somatic, ExonicFunc.refGene != "unknown" )
somatic<-subset(somatic, ExonicFunc.refGene != "" )

somatic_Samples.resp$dMMR <-"sem_dMMR"
mutados<- unique(subset(somatic, class == "dMMR", select = ID_Exoma))
somatic_Samples.resp[somatic_Samples.resp$ID_Exoma %in% mutados$ID_Exoma]$dMMR <-"com_dMMR"


# CARREGANDO INFORMAÇÕES DE SVs
#===================================#
sv<- fread(file="OUTPUT-SV_Filter_R.2023-11-08/tabela_count.SVs.tsv", header=T, sep="\t", fill=TRUE, check.names = FALSE)
sv<- as.data.frame(sv)
sv<-sv %>% 
  mutate_all(replace_na, 0)

somatic_Samples.resp.sv <- somatic_Samples.resp

somatic_Samples.resp <- merge( somatic_Samples.resp, sv, by.x = "ID_Exoma", by.y = "sample", all.x = T)
somatic_Samples.resp$SV[is.na(somatic_Samples.resp$SV)] <- 0
somatic_Samples.resp$BND[is.na(somatic_Samples.resp$BND)] <- 0
somatic_Samples.resp$DEL[is.na(somatic_Samples.resp$DEL)] <- 0
somatic_Samples.resp$DUP[is.na(somatic_Samples.resp$DUP)] <- 0
somatic_Samples.resp$INV[is.na(somatic_Samples.resp$INV)] <- 0
somatic_Samples.resp$status_SV <- "sem SV"
somatic_Samples.resp[somatic_Samples.resp$SV > 0, "status_SV"] <- "com SV"


aux <-somatic_Samples.resp

print(table(aux$Resposta))
print(as.matrix(table(aux[,c("status_SV", "Resposta")]), nr=2))
dt3<- matrix(table(aux[,c("status_SV", "Resposta")]), nr=2)
print(chisq.test(dt3))
print(fisher.test(dt3))

print(table(aux$Resposta2))
print(as.matrix(table(aux[,c("status_SV", "Resposta2")]), nr=2))
dt3<- matrix(table(aux[,c("status_SV", "Resposta2")]), nr=2)
print(chisq.test(dt3))
print(fisher.test(dt3))

print(table(aux$Resposta3))
print(as.matrix(table(aux[,c("status_SV", "Resposta3")]), nr=2))
dt3<- matrix(table(aux[,c("status_SV", "Resposta3")]), nr=2)
print(chisq.test(dt3))
print(fisher.test(dt3))


#==============================================================================
# PACIENTES COM INV
somatic_Samples.resp$status_INV <- "sem INV"
somatic_Samples.resp[somatic_Samples.resp$INV > 0, "status_INV"] <- "com INV"


aux <-somatic_Samples.resp
print(table(aux$Resposta))
print(as.matrix(table(aux[,c("status_INV", "Resposta")]), nr=2))
dt3<- matrix(table(aux[,c("status_INV", "Resposta")]), nr=2)
print(fisher.test(dt3))

print(table(aux$Resposta2))
print(as.matrix(table(aux[,c("status_INV", "Resposta2")]), nr=2))
dt3<- matrix(table(aux[,c("status_INV", "Resposta2")]), nr=2)
print(fisher.test(dt3))

print(table(aux$Resposta3))
print(as.matrix(table(aux[,c("status_INV", "Resposta3")]), nr=2))
dt3<- matrix(table(aux[,c("status_INV", "Resposta3")]), nr=2)
print(fisher.test(dt3))


#==============================================================================
# PACIENTES COM SV
#==============================================================================
sv<- sv %>% 
  mutate(Gene.refGene = strsplit(as.character(Gene.refGene), ";")) %>% 
  unnest(Gene.refGene)

sv.genes=unique(select(sv, c("Gene.refGene","variable", "type")))
sv.genes.heatmap= as.data.frame(table(unique(select(sv, c("Gene.refGene","variable")))))

tab<-sv[c("Gene.refGene","variable")]
tab<-as.data.frame(table(tab))
tab<-as.matrix(dcast(tab, Gene.refGene ~ variable, value.var = "Freq"))
rownames(tab) <- tab[,1]
tab<-type.convert(tab[,-1], as.is = T)

sv_type<-sv[c("Gene.refGene","type")]
sv_type<-as.data.frame(table(sv_type))
sv_type<-dcast(sv_type, Gene.refGene ~ type, value.var = "Freq")
row.names(sv_type) <- sv_type[,1]
sv_type<-subset(sv_type, select=-1)

my_sample_col <- data.frame(subset(somatic_Samples.resp, 
                                   select = c("sample","Resposta3")))
row.names(my_sample_col) <- my_sample_col[,1]
my_sample_col<-subset(my_sample_col, select=-1)

pheatmap((tab), color=colorRampPalette(c("white", "darkred"))(100),
         annotation_col = my_sample_col, annotation_row = sv_type, legend = F ,
         cluster_cols  = T, cluster_rows  = T)

dev.off()

#==============================================================================#
#         Chi.sq - Verificando a correlação de alguma assinatura
#==============================================================================#

aux <-somatic_Samples.resp

print(table(aux$cosmic_sig))
print(as.matrix(table(aux[,c("dMMR", "cosmic_sig")]), nr=2))
dt3<- matrix(table(aux[,c("dMMR", "cosmic_sig")]), nr=2)
print(fisher.test(dt3))
print(chisq.test(dt3))

aux<-as.data.frame(table(subset(final_Ass,  select=c("sample","ass_id","dMMR", "Resposta3"))), stringsAsFactors = F)
aux <-somatic_Samples.resp
List_Ass <-unique(final_Ass$ass_id)

for (id in List_Ass) {
  # id="6"
  print(id)
  
  aux <-somatic_Samples.resp
  aux$dMMR <-"sem_dMMR"
  mutados<- unique(subset(somatic, class == "dMMR", select = ID_Exoma))
  aux[aux$ID_Exoma %in% mutados$ID_Exoma]$dMMR <-"com_dMMR"
  aux$ass_id<-"0"
  ass_sample<-subset(final_Ass,  ass_id==id ,select=c("sample"))
  aux[aux$sample %in% ass_sample$sample]$ass_id <-id
  
  print(table(aux$Resposta))
  print(as.matrix(table(aux[,c("ass_id", "Resposta")]), nr=2))
  dt3<- matrix(table(aux[,c("ass_id", "Resposta")]), nr=2)
  print(fisher.test(dt3))
  print(chisq.test(dt3))
  print("=====================================================================")
}



