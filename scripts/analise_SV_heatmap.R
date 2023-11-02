
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


wd <- "/media/vande/HDD/PROJETOS-HSL_BP/resultados_Exoma/"
setwd(wd)

#                                        ===================================
#                                                ANALISE DOS RESULTADOS
#                                        ===================================

##################################################
### CARREGANDO TABELAS DE RESULTADOS
##################################################
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
sv<- fread(file="../resultados_SV/tabela_count.SVs.tsv", header=T, sep="\t", fill=TRUE, check.names = FALSE)
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



