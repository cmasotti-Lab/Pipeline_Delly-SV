
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
library(StructuralVariantAnnotation)

wd <- "D:/Github-Projects/Pipeline_Delly-SV/"
setwd(wd)
#hoje<-Sys.Date()
hoje<-"2023-11-10"
output <- paste0("OUTPUT-SV_Filter_R.",hoje)
dir.create(output)
#                                        ===================================
#                                                ANALISE DOS RESULTADOS
#                                        ===================================


#==============================================================================#
# STEP 0 -  Datasets reorganization ####
#==============================================================================#
all.sv<- fread(file="Result_Delly_SV.2023-10-28/step6_annovar/ROP-annovar.hg38_multianno.txt", header=T, sep="\t", na.strings=c(".","NA"), fill=TRUE, check.names = FALSE)
samples_ID <-fread("Result_Delly_SV.2023-10-28/step6_annovar/sampleID.list", header = F)

# Rename columns
#==============================================================================#
colnames(all.sv) <-  c( colnames(all.sv[,c(1:11)]),"Otherinfo1","Otherinfo2","Otherinfo3","CHR","POS","ID","REF","ALT","QUAL","FILTER", "INFO", "FORMAT", 
                          samples_ID$V1)

# Creating INDEX com ID and selecting the main columns
#==============================================================================#
all.sv$index <- all.sv$ID
all.sv <-data.frame(all.sv)
exCol=colnames(all.sv[,-which(names(all.sv) %in% c("Otherinfo1","Otherinfo2","Otherinfo3","ID","index"))])
all.sv <- dplyr::select(all.sv, index,all_of(exCol))
dim(all.sv)   #55882   105
#fwrite(all.sv, paste0(output,"/tab0.tsv"), quote = F, sep="\t")

#==============================================================================#
#                 STEP 1 - Clear dataset ####
#==============================================================================#
tab1<- all.sv

## Restructure tables, all samples in the same column ###
#==============================================================================#
tab1 <- reshape2::melt(tab1,id=colnames(tab1[,1:which(names(tab1)== "FORMAT")]))
colnames(tab1)[colnames(tab1) == "variable"] <- "SAMPLE"

unique(tab1$FORMAT) #GT:GL:GQ:FT:RCL:RC:RCR:RDCN:DR:DV:RR:RV
FORMAT<-strsplit(unique(tab1$FORMAT), ":")[[1]]
tab2<- separate(data = tab1, col = value, into =FORMAT, sep = ":")

table(tab2$GT)
tab2 <- tab2[!grepl("^0/0", tab2$GT),]
tab2 <- tab2[!grepl("^0\\|0", tab2$GT),]
tab2 <- tab2[!grepl("^\\./\\.", tab2$GT),]
table(tab2$GT)

remove(all.sv, tab1)
gc()

#==============================================================================#
### SEPARAR OS ALELOS DE CADA INDIVIDUO
#==============================================================================#
setDT(tab2)[,paste0("GT", 1:2) := tstrsplit(GT, "/")]
#==============================================================================#

tab2$ALL <- substr(tab2$index, 1, 3)
table(tab2$ALL)
table(tab2$GT)

#==============================================================================#
### AMOSTRAS EXCLUIDAS POR DADOS CLINICOS
#==============================================================================#
tab2$SAMPLE <- gsub("\\.", "-", tab2$SAMPLE)

samples.excl <- c('ROP-116','ROP-118','ROP-120','ROP-121','ROP-122','ROP-123','ROP-125','ROP-126',
                  'ROP-127','ROP-130','ROP-131','ROP-132','ROP-23','ROP-84', 'ROP-81','ROP-91','ROP-129')
#                 'ROP-83','ROP-107',)
samples.exclude <-paste(samples.excl, collapse = "|" )

tab2 <- tab2[!grepl(samples.exclude, tab2$SAMPLE),]
tab.info <- as.data.frame(tab2$INFO)

#==============================================================================#
#                 STEP 2 - FILTERS ####
#==============================================================================#

#==============================================================================#
### SELECIONAR OS PRECISE & PASS
#==============================================================================#
tab3 <- tab2[grepl("^PRECISE", tab2$INFO),]
table(tab3$ALL)

tab3 <- subset(tab3, FT == "PASS") 
table(tab3$ALL)

# Removendo mutações compartilhadas >=2
#==============================================================================#
aux <- data.frame(table(tab3$index))
tab3 <- merge(tab3, aux, by.x = "index", by.y = "Var1")
colnames(tab3)[names(tab3)=="Freq"]<-"MUTATION_shared"
shared2exclude <- which((tab3$MUTATION_shared >= 2) )
nrow(unique(tab3[shared2exclude,"index"]))  #8565
tab3excl <- tab3[shared2exclude,] 
tab3 <- tab3[-shared2exclude,] 
table(tab3$ALL)
# BND  DEL  DUP  INV 
# 48  241   26 2961

# Remover DEL em região: "intronic","intergenic","UTR3","UTR5","upstream;downstream","ncRNA_intronic","UTR5;UTR3","upstream" 
#==============================================================================#
data.table(table(tab3[tab3$ALL == "DEL",]$Func.refGene))
unique(tab3[tab3$ALL == "DEL",]$Func.refGene)

tab4 <- subset(tab3, !(ALL == "DEL" & Func.refGene %in% c("intronic","intergenic","UTR3","UTR5","upstream;downstream","ncRNA_intronic","UTR5;UTR3","upstream" )))
#==============================================================================#
#            IMPORTANTE           ####
# APLICAR FILTRO PARA DUP, INV e INS ?
tab4 <- subset(tab4, !(ALL == "DUP" & Func.refGene %in% c("intronic","intergenic","UTR3","UTR5","upstream;downstream","ncRNA_intronic","UTR5;UTR3","upstream" )))
tab4 <- subset(tab4, !(ALL == "INS" & Func.refGene %in% c("intronic","intergenic","UTR3","UTR5","upstream;downstream","ncRNA_intronic","UTR5;UTR3","upstream" )))
tab4 <- subset(tab4, !(ALL == "INV" & Func.refGene %in% c("intronic","intergenic","UTR3","UTR5","upstream;downstream","ncRNA_intronic","UTR5;UTR3","upstream" )))
#==============================================================================#

table(tab3$ALL,tab3$Func.refGene)
table(tab4$ALL,tab4$Func.refGene)

table(tab4$ALL)
# BND  DEL  DUP  INV 
# 48  117   25 2391
table(tab4$SAMPLE,tab4$ALL)

#==============================================================================#
#   Excluindo SV com anotação no DGV           ####
#==============================================================================#
tab5<-tab4[is.na(tab4$dgvMerged)]
table(tab5$ALL)
table(tab5$SAMPLE,tab5$ALL)

tab.count <- data.frame(table(tab5$SAMPLE,tab5$ALL))
tab.count<- as.data.frame(reshape2::dcast(tab.count, Var1 ~ Var2,  value.var = "Freq"))

tab.count$ALL <- rowSums(tab.count[,c(2:5)])
tab.count
sum(tab.count$ALL) #460


fwrite(tab4, paste0(output,"/tabela_all.SVs.tsv"), quote = F, sep="\t")
fwrite(tab5, paste0(output,"/tabela_all.SVs.exclDGV.tsv"), quote = F, sep="\t")
fwrite(tab.count, paste0(output,"/tabela_count.SVs.tsv"), quote = F, sep="\t")

#==============================================================================#
# CARREGANDO DADOS DO MUTECT2 E DADOS CLINICOS ####
#==============================================================================#
# carregando dados do Mutect2
#==============================================================================#
dgenes<-fread("D:/PROJETOS-HSL_BP/resultados_Mutect2-2023/input_Somatic-Filter/driver-genes.txt", 
              header = F, col.names = "gene")
dMMR <- fread("D:/PROJETOS-HSL_BP/resultados_Mutect2-2023/input_Somatic-Filter/dMMR_genes.txt", 
            header = F, col.names = "gene")
somatic.69<-fread("D:/PROJETOS-HSL_BP/resultados_Mutect2-2023/output-Mutect2_somaticFilters.2023-09-26/somatic_mutation.69samples.tsv",
                  header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
somatic.69<- somatic.69 %>% 
  mutate(Gene.refGene = strsplit(as.character(Gene.refGene), ";")) %>% 
  unnest(Gene.refGene)

somatic_Samples<- fread("D:/PROJETOS-HSL_BP/resultados_Mutect2-2023/output-Mutect2_SomaticFilters.2023-09-26/MATH_TMB-SAMPLES.tsv", 
                        header = T,  na.strings=c("NA"), fill=TRUE, check.names = FALSE)
somatic_Samples<-subset(somatic_Samples, (SAMPLE %in% somatic.69$SAMPLE) )


#==============================================================================#
# Carregando tabela de dados Clinicos ####
#==============================================================================#
Clinical <- read.xlsx("D:/PROJETOS-HSL_BP/resultados_Mutect2-2023/Dados Clinicos JP final  - Dez2022 - FILTRADO.xlsx", 
                      sheetIndex = 1, header = T)
Clinical<- Clinical[rowSums(is.na(Clinical)) != ncol(Clinical), -c(6)]
## Renomeando desfeichos clinicos ####
#==============================================================================#
Clinical$Response1 <- gsub("Completa", "nCRT-R" , Clinical$Resposta)
Clinical$Response1 <- gsub("Incompleta", "nCRT-NR" , Clinical$Response1)
Clinical$Response2 <- gsub("Completa", "nCRT-R" , Clinical$Resposta2)
Clinical$Response2 <- gsub("Incompleta Sem Metastase", "nCRT-NR Not-Metastatic" , Clinical$Response2)
Clinical$Response2 <- gsub("Incompleta Com Metastase", "nCRT-NR Metastatic" , Clinical$Response2)
Clinical$Response3 <- gsub("nao metastatico", "Not-Metastatic" , Clinical$Resposta3)
Clinical$Response3 <- gsub("metastatico", "Metastatic" , Clinical$Response3)
#==============================================================================#
## Selecionando amostas do Exoma 
Clinical <- Clinical[Clinical$WXS== "Y",]
## Listras de amostras excluidas pelo corpo clinico 
samples.excl <- c("ROP-116","ROP-118","ROP-120","ROP-121","ROP-122","ROP-123","ROP-125","ROP-126",
                  "ROP-127","ROP-130","ROP-131","ROP-132","ROP-23","ROP-84", "ROP-129","ROP-81","ROP-91")
## Excluindo amostras excluidas pelo corpo clinico 
Clinical <- subset(Clinical, !(ID_Exoma %in% samples.excl) )


somatic.69_Clinical <- merge(somatic_Samples, Clinical, by.x = "SAMPLE", by.y = "Sample")
somatic.69_Clinical$dMMR <-"sem_dMMR"
mutados<- unique(subset(somatic.69, class == "dMMR", select = ID_Exoma))
somatic.69_Clinical[somatic.69_Clinical$ID_Exoma %in% mutados$ID_Exoma]$dMMR <-"com_dMMR"


# CARREGANDO INFORMAÇÕES DE SVs
#===================================#
sv.events<- fread(file=paste0(output,"/tabela_all.SVs.tsv"), header=T, sep="\t", fill=TRUE, check.names = FALSE)

sv<- fread(file=paste0(output,"/tabela_count.SVs.tsv"), header=T, sep="\t", fill=TRUE, check.names = FALSE)
sv<- as.data.frame(sv)
sv<-sv %>% 
  mutate_all(replace_na, 0)

somatic.69_Clinical.sv <- somatic.69_Clinical

somatic.69_Clinical <- merge( somatic.69_Clinical, sv, by.x = "SAMPLE", by.y = "Var1", all.x = T)
somatic.69_Clinical$ALL[is.na(somatic.69_Clinical$ALL)] <- 0
somatic.69_Clinical$BND[is.na(somatic.69_Clinical$BND)] <- 0
somatic.69_Clinical$DEL[is.na(somatic.69_Clinical$DEL)] <- 0
somatic.69_Clinical$DUP[is.na(somatic.69_Clinical$DUP)] <- 0
somatic.69_Clinical$INV[is.na(somatic.69_Clinical$INV)] <- 0
somatic.69_Clinical$status_SV <- "sem SV"
somatic.69_Clinical[somatic.69_Clinical$ALL > 0, "status_SV"] <- "com SV"


#==============================================================================
# TESTE QUI-QUADRADO E FISHER
#==============================================================================
aux <-somatic.69_Clinical

aux$status_INV <- "sem_SV"
aux$status_DUP <- "sem_SV"
aux$status_DEL <- "sem_SV"
aux$status_BND <- "sem_SV"

aux[aux$INV >0, ]$status_INV <- "com_SV"
aux[aux$DUP >0, ]$status_DUP <- "com_SV"
aux[aux$DEL >0, ]$status_DEL <- "com_SV"
aux[aux$BND >0, ]$status_BND <- "com_SV"

aux <-aux
eventos_SV <- c("status_SV","status_INV","status_DUP","status_DEL","status_BND")
respostas <- c("Response1", "Response2", "Response3")

for (sv in eventos_SV) {
  for (resp in respostas) {
    # sv="status_SV"
    # resp= "Response1"
    
    dt3<- matrix(table(subset(aux, select = c(sv, resp))), nr=2)
    dt3
    if(ncol(dt3)>1){
      resultado<-fisher.test(dt3)
      resultado1<-chisq.test(dt3)
      valor_p <- resultado$p.value
      valor_p1 <- resultado1$p.value
      if((valor_p <= 0.05)||(valor_p1 <= 0.05)){
        print("=============================================================")
        print(sv)
        print(resp)
        print(as.matrix(table(subset(aux, select = c(sv, resp))), nr=2))
        #print("=============================================================")
        print(resultado)
        #print("=============================================================")
        print(chisq.test(dt3))
        print("=============================================================")
      }
    }
  }
}




# ATE AQUI ESTA OK
#==============================================================================
# PACIENTES COM SV
#==============================================================================
somatic.69<- somatic.69 %>% 
  mutate(Gene.refGene = strsplit(as.character(Gene.refGene), ";")) %>% 
  unnest(Gene.refGene)
merge(somatic.69,subset(somatic.69_Clinical.sv, ))

sv.genes=unique(subset(somatic.69, select=c("Gene.refGene","variable", "type")))
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

my_sample_col <- data.frame(subset(somatic.69_Clinical, 
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

aux <-somatic.69_Clinical

print(table(aux$cosmic_sig))
print(as.matrix(table(aux[,c("dMMR", "cosmic_sig")]), nr=2))
dt3<- matrix(table(aux[,c("dMMR", "cosmic_sig")]), nr=2)
print(fisher.test(dt3))
print(chisq.test(dt3))

aux<-as.data.frame(table(subset(final_Ass,  select=c("sample","ass_id","dMMR", "Resposta3"))), stringsAsFactors = F)
aux <-somatic.69_Clinical
List_Ass <-unique(final_Ass$ass_id)

for (id in List_Ass) {
  # id="6"
  print(id)
  
  aux <-somatic.69_Clinical
  aux$dMMR <-"sem_dMMR"
  mutados<- unique(subset(somatic.69, class == "dMMR", select = ID_Exoma))
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




#==============================================================================#
#         TESTANDO - StructuralVariantAnnotation
#==============================================================================
#suppressPackageStartupMessages(library(StructuralVariantAnnotation))
library(VariantAnnotation)

vcf <- readVcf("Result_Delly_SV.2023-10-28/ROP-Filter.vcf",
               genome = "hg38")
gr <- c(breakpointRanges(vcf), breakendRanges(vcf))
gr

library(VariantAnnotation)

vcf <- readVcf("D:/Github-Projects/Pipeline_Delly-SV/Result_Delly_SV.2023-10-28/ROP-Filter.vcf",
               genome = "hg38")

header(vcf)



library("intansv")

delly <- readDelly("Result_Delly_SV.2023-10-28/ROP-Filter.vcf")

sv_all_methods <- methodsMerge(delly,delly)
str(sv_all_methods)

anno.file.path <- system.file("extdata/chr05_chr10.anno.txt", package="intansv")
anno.file.path

msu_gff_v7 <- read.table(anno.file.path, head=TRUE, as.is=TRUE)

head(msu_gff_v7)

sv_all_methods.anno <- llply(delly,svAnnotation,genomeAnnotation=msu_gff_v7)
names(sv_all_methods.anno)

head(sv_all_methods.anno$del)

genome.file.path <- system.file("extdata/chr05_chr10.genome.txt", package="intansv")
genome <- read.table(genome.file.path, head=TRUE, as.is=TRUE)
plotChromosome(genome, delly,10000000)

head(msu_gff_v7, n=3)
plotRegion(delly,msu_gff_v7, "chr10", 1, 20000000)
