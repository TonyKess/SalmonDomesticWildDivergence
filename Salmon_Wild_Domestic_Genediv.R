setwd("~/Desktop/Projects/SalmonAQIntrogression2023/")
library(tidyverse)
library(data.table)
library(marmap)
library(ggman)
library(marmap)



#read in complete metadata, match to genotyped individuals
Meta <- fread("CIGENE_220K_Metadata_All_Aug_15_2023.txt") 

Fam <- fread("Ssal.31.220K.mapped.June2023merge.fam") %>%  
  mutate(FID = V1, IID = V2) %>% 
  select(FID, IID)

Joined <- inner_join(Fam, Meta)

Joined %>%  select(Dataset) %>%  distinct()

Chromconv <- fread("Ssal.31.chromconversion") %>%  mutate(Chromnumb = str_replace(Chromosome, "ssa0", "")) %>% 
  mutate(Chromnumb = str_replace(Chromnumb, "ssa", ""))

IcelandMeta <- fread("IceNorMeta.tsv")
Iceland_fam <- fread("IceAQ.PolyHigh.fam") %>% mutate(FID = V1, IID = V2) %>%  select(FID, IID)
Iceland_Joined <- inner_join(IcelandMeta, Iceland_fam)
IceNorMetaLatLon <- fread("IceNorMetawithLatLon.tsv")

Iceland_Joined <- inner_join(Iceland_Joined, IceNorMetaLatLon)


IcelandRiversTable60K <- Iceland_Joined %>% 
  filter(Group %in% "Iceland") %>% 
  mutate(Rivername = Population_River_breedingkernel) %>% 
  group_by(Rivername, Lat, Long) %>% 
  summarise(Count = n(), .groups = "drop")

write.table(IcelandRiversTable60K, "IcelandRiversTable60K.tsv", col.names = T, row.names = F, sep = "\t", quote = F)


#First we want to do some spatial assignments - we know from past work these are pretty accurate in terms of genetic assignment.
#We will PCA them too to confirm.
Salmon_NorthAM <- Joined %>%  filter(Region %in% "NorthAm")
Salmon_Norway <- Joined %>%  filter(Region %in% "Norway")

Salmon_NorthAM$Lat <- as.numeric(Salmon_NorthAM$Lat)
Salmon_NorthAM$Lon <- as.numeric(Salmon_NorthAM$Lon)
Salmon_Norway$Lat <- as.numeric(Salmon_Norway$Lat)
Salmon_Norway$Lon <- as.numeric(Salmon_Norway$Lon)


Salmon_NorthAM <- Salmon_NorthAM %>%  drop_na()
Salmon_Norway <- Salmon_Norway %>%  drop_na()

Clusters <- fread("SalmonRiverCluster.tsv") %>%  mutate(FID = location )


#We also want to split NA here into subgroups
Salmon_NorthAM_refset <- Salmon_NorthAM %>%  filter(Region %in%"NorthAm",
                                                Dataset %in% c("May 2021 5700 Recluster", "June 2021 Baseline"))

NothAmericaRiversTable220K <- Salmon_NorthAM_refset %>% 
  group_by(Rivername, Lat, Lon) %>% 
  summarise(Count = n(), .groups = "drop")

write.table(NothAmericaRiversTable220K, "NorthAmericanRiversTable.tsv", col.names = T, row.names = F, sep = "\t", quote = F)

#We also want to split NA here into subgroups
Salmon_Norway_refset <- Salmon_Norway %>%  filter(Region %in%"Norway",
                                                    Dataset %in% c("May 2021 5700 Recluster", "June 2021 Baseline"))

NorwayRiversTable220K<- Salmon_Norway_refset %>% 
  group_by(Rivername, Lat, Lon) %>% 
  summarise(Count = n(), .groups = "drop")

write.table(NorwayRiversTable220K, "NorwayRiversTable.tsv", col.names = T, row.names = F, sep = "\t", quote = F)



#map
Coords_NA_Nor <- bind_rows(Salmon_NorthAM_refset, Salmon_Norway_refset) %>%  select(FID, Region, Lat, Lon) %>%  distinct()
range(Coords_NA_Nor$Lat)
range(Coords_NA_Nor$Lon)
#map
bathydata <- getNOAA.bathy(-71,32,40,72, res=3,keep=T)
plot(bathydata)
map=autoplot.bathy(bathydata, geom=c("r", "c")) +
  scale_fill_etopo(guide = FALSE) +
  geom_contour(aes(z=z), breaks=c(-100, -200, -500, -1000, -2000, -4000), colour="grey90", size=0.01) +
  xlab("Degrees Longitude") +
  ylab("Degrees Latitude") +
  theme(legend.position = "right")

map + geom_point(data = Coords_NA_Nor, aes (y = Lat, x = Lon, colour = Region))


#Aquacultures
Aquaculture_NA <- Meta  %>%  filter(Region %in% "Aquaculture",
                          Dataset %in% c("May 2021 5700 Recluster", "June 2021 Baseline"))

write.table(Aquaculture_NA)

Aquaculture_Norway <- Meta  %>%  filter(Region %in% "Aquaculture_Norway")



AQ_NA_clusterfile <- bind_rows(Salmon_NorthAM_refset %>%  select(FID, IID, Region), Aquaculture_NA %>%  select(FID, IID, Region))
AQ_Norway_clusterfile <- bind_rows(Salmon_Norway_refset %>%  select(FID, IID, Region), Aquaculture_Norway %>%  select(FID, IID, Region))

write.table(AQ_NA_clusterfile, "AQ_NA_clusterfile.tsv", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(AQ_Norway_clusterfile, "AQ_Norway_clusterfile.tsv", col.names = F, row.names = F, quote = F, sep = "\t")

system("source ~/.zshrc
       conda activate plink
       plink --bfile Ssal.31.220K.mapped.June2023merge --maf 0.01 --fst  --keep AQ_NA_clusterfile.tsv --within AQ_NA_clusterfile.tsv --allow-extra-chr --make-bed --out Ssal.31.220K.mapped.June2023merge.NA.AQ")

system("source ~/.zshrc
       conda activate plink
       plink --bfile Ssal.31.220K.mapped.June2023merge --maf 0.01 --fst  --keep AQ_Norway_clusterfile.tsv --within AQ_Norway_clusterfile.tsv --allow-extra-chr --make-bed --out Ssal.31.220K.mapped.June2023merge.Norway.AQ")

system("source ~/.zshrc
       conda activate plink
       plink --bfile IceAQ --maf 0.01 --fst --within IceandAQ.groups.tsv  --extract ICESNPPoly.tsv --allow-extra-chr --chr-set 30  --out IceAQ.PolyHigh")



NA_AQ.fst <- fread("Ssal.31.220K.mapped.June2023merge.NA.AQ.fst")
NA.AQ.fst.outliers_top1K <-  NA_AQ.fst %>% slice_max(FST, n = 1000)
write.table(NA.AQ.fst.outliers_top1K  %>%  select(SNP),"NA.AQ.fst.outliers_top1KSNPs.tsv", col.names = F, row.names = F, sep = "\t", quote = F)


NA.AQ.fst.outliers_top1K.bed <- NA.AQ.fst.outliers_top1K %>% mutate(POS_left = POS-5000, POS_right = POS + 5000) %>%  select(CHR, POS_left, POS_right, SNP)
write.table(NA.AQ.fst.outliers_top1K.bed, "NA.AQ.fst.outliers_top1K.bed", col.names = F, row.names = F, sep = "\t", quote = F)
system("source ~/.zshrc
       conda activate align
       bedtools intersect -a Ssal_31.genes.bed -b NA.AQ.fst.outliers_top1K.bed -wb > NA.AQ.fstoutliers.top1K.geneverlap.tsv")

NA.AQ.fstoutliers.top1K.geneverlap <- fread("NA.AQ.fstoutliers.top1K.geneverlap.tsv") %>%  
  mutate(gene.id = V4, CHROM = V5, POS = V6+5000, SNP = V8) %>%  
  select(CHROM, POS, SNP, gene.id) %>% 
  inner_join(NA_AQ.fst)

NA.AQ.fstoutliers.top1K.geneverlap.genelist <- NA.AQ.fstoutliers.top1K.geneverlap  %>%  select(gene.id) %>% distinct()
write.table(NA.AQ.fstoutliers.top1K.geneverlap.genelist, "NA.AQ.fstoutliers.top1K.geneverlap.genelist", col.names = F, row.names = F, sep = "\t", quote = F)

system("grep -fNA.AQ.fstoutliers.top1K.geneverlap.genelist Ssal_31.genes.bed  >NorthAmerica_outliergenes.tsv")



Norway_AQ.fst <- fread("Ssal.31.220K.mapped.June2023merge.Norway.AQ.fst")
Norway.AQ.fst.outliers_top1K <-  Norway_AQ.fst%>% slice_max(FST, n = 1000)
write.table(Norway.AQ.fst.outliers_top1K %>%  select(SNP),"Norway.AQ.fst.outliers_top1KSNPs.tsv", col.names = F, row.names = F, sep = "\t", quote = F)

Norway.AQ.fst.outliers_top1K.bed <- Norway.AQ.fst.outliers_top1K%>% mutate(POS_left = POS-5000, POS_right = POS + 5000) %>%  select(CHR, POS_left, POS_right, SNP)
write.table(Norway.AQ.fst.outliers_top1K.bed, "Norway.AQ.fst.outliers_top1K.bed", col.names = F, row.names = F, sep = "\t", quote = F)
system("source ~/.zshrc
       conda activate align
       bedtools intersect -a Ssal_31.genes.bed -b Norway.AQ.fst.outliers_top1K.bed -wb > Norway.AQ.fstoutliers.top1K.geneverlap.tsv")



Norway.AQ.fstoutliers.top1K.geneverlap <- fread("Norway.AQ.fstoutliers.top1K.geneverlap.tsv") %>%  
  mutate(gene.id = V4, CHROM = V5, POS = V6 + 5000, SNP = V8) %>%  
  select(CHROM, POS, SNP, gene.id) %>% 
  inner_join(Norway_AQ.fst)

Norway.AQ.fstoutliers.top1K.geneverlap.genelist <- Norway.AQ.fstoutliers.top1K.geneverlap %>%  select(gene.id) %>% distinct()
write.table(Norway.AQ.fstoutliers.top1K.geneverlap.genelist , "Norway.AQ.fstoutliers.top1K.geneverlap.genelist", col.names = F, row.names = F, sep = "\t", quote = F)
system("grep -f Norway.AQ.fstoutliers.top1K.geneverlap.genelist Ssal_31.genes.bed  >Norway_outliergenes.tsv")


Ice_AQ.fst <- fread("IceAQ.PolyHigh.fst")
Ice.AQ.fst.outliers_top1K <-  Ice_AQ.fst %>% slice_max(FST, n = 1000)
write.table(Ice.AQ.fst.outliers_top1K  %>%  select(SNP),"Ice.AQ.fst.outliers_top1KSNPs.tsv", col.names = F, row.names = F, sep = "\t", quote = F)

SSA_genes <- fread("SSA_genes.bed") %>% mutate(Chromosome = V1) %>%  inner_join(Chromconv) %>%  select(Chromnumb, V2, V3, V4)
write.table(SSA_genes, "SSA_genes_IL_temp.bed", col.names = F, row.names = F, sep = "\t", quote = F)

Ice.AQ.fst.outliers_top1K.bed <- Ice.AQ.fst.outliers_top1K %>% mutate(POS_left = POS-5000, POS_right = POS + 5000, Chromnumb = CHR) %>% 
  select(Chromnumb, POS_left, POS_right, SNP) 
write.table(Ice.AQ.fst.outliers_top1K.bed, "Ice.AQ.fst.outliers_top1K.bed", col.names = F, row.names = F, sep = "\t", quote = F)
system("source ~/.zshrc
       conda activate align
       bedtools intersect -a SSA_genes_IL_temp.bed -b Ice.AQ.fst.outliers_top1K.bed -wb > Ice.AQ.fstoutliers.top1K.geneverlap.tsv")

Ice.AQ.fstoutliers.top1K.geneverlap <- fread("Ice.AQ.fstoutliers.top1K.geneverlap.tsv") %>%  
  mutate(gene.id = V4, CHROM = V5, POS = V6+5000, SNP = V8) %>%  
  select(CHROM, POS, SNP, gene.id) %>% 
  inner_join(Ice_AQ.fst)

Ice.AQ.fstoutliers.top1K.geneverlap.genelist <- Ice.AQ.fstoutliers.top1K.geneverlap  %>%  select(gene.id) %>% distinct()
write.table(Ice.AQ.fstoutliers.top1K.geneverlap.genelist, "Ice.AQ.fstoutliers.top1K.geneverlap.genelist", col.names = F, row.names = F, sep = "\t", quote = F)

system("grep -f Ice.AQ.fstoutliers.top1K.geneverlap.genelist SSA_genes_IL_temp.bed > Iceland_outliergenes.tsv")


inner_join(NA.AQ.fstoutliers.top1K.geneverlap.genelist , Norway.AQ.fstoutliers.top1K.geneverlap.genelist)
inner_join(Ice.AQ.fstoutliers.top1K.geneverlap.genelist , Norway.AQ.fstoutliers.top1K.geneverlap.genelist)
inner_join(NA.AQ.fstoutliers.top1K.geneverlap.genelist ,Ice.AQ.fstoutliers.top1K.geneverlap.genelist)


Norway.AQ.fstoutliers.top1K.geneverlap.noloc <- Norway.AQ.fstoutliers.top1K.geneverlap  %>%  filter(!str_starts(gene.id, "LOC"))
NA.AQ.fstoutliers.top1K.geneverlap.noloc <- NA.AQ.fstoutliers.top1K.geneverlap  %>%  filter(!str_starts(gene.id, "LOC"))
Ice.AQ.fstoutliers.top1K.geneverlap.noloc <- Ice.AQ.fstoutliers.top1K.geneverlap  %>%  filter(!str_starts(gene.id, "LOC"))


#for plotting

NA_AQ.fst_10K_downsamp <- NA_AQ.fst %>%  filter(!SNP %in% NA.AQ.fst.outliers_top1K$SNP) %>%  sample_n(10000) %>% bind_rows(NA.AQ.fst.outliers_top1K) %>% inner_join(Chromconv) 
Norway_AQ.fst_10K_downsamp <- Norway_AQ.fst %>%  filter(!SNP %in% Norway.AQ.fst.outliers_top1K$SNP) %>%  sample_n(10000) %>% bind_rows(Norway.AQ.fst.outliers_top1K)  %>% inner_join(Chromconv) 
Iceland_AQ.fst_10K_downsamp <- Ice_AQ.fst %>%  filter(!SNP %in% Ice.AQ.fst.outliers_top1K$SNP) %>%  sample_n(10000) %>% bind_rows(Ice.AQ.fst.outliers_top1K) %>%  mutate(Chromnumb = as.character(CHR)) %>%  select(-CHR) %>% inner_join(Chromconv) 

NAdsplot <- ggman(gwas = NA_AQ.fst_10K_downsamp , 
                      snp = "SNP", 
                      chrom = "Chromosome", 
                      bp = "POS", 
                      pvalue = "FST", 
                      logTransform = F, 
                      sigLine = 0, 
                      pointSize = 2.5, 
                      ymin = 0, 
                      ymax = 1, 
                      ylabel = "FST", 
                      xlabel = "Chromosome", 
                      title = "North America wild vs Aquaculture"
) + theme_classic()
NAdsplot<- ggmanHighlight(NAdsplot, highlight = NA.AQ.fst.outliers_top1K$SNP)

ggmanLabel(NAdsplot, labelDfm = NA.AQ.fstoutliers.top1K.geneverlap.noloc, snp = "SNP", label = "gene.id")


Norwaydsplot <- ggman(gwas = Norway_AQ.fst_10K_downsamp , 
                  snp = "SNP", 
                  chrom = "Chromosome", 
                  bp = "POS", 
                  pvalue = "FST", 
                  logTransform = F, 
                  sigLine = 0, 
                  pointSize = 2.5, 
                  ymin = 0, 
                  ymax = 1, 
                  ylabel = "FST", 
                  xlabel = "Chromosome", 
                  title = "Norway wild vs Aquaculture"
) + theme_classic()

Norwaydsplot<- ggmanHighlight(Norwaydsplot, highlight = Norway.AQ.fst.outliers_top1K$SNP)

ggmanLabel(Norwaydsplot, labelDfm = Norway.AQ.fstoutliers.top1K.geneverlap.noloc, snp = "SNP", label = "gene.id")




Icedsplot <- ggman(gwas = Iceland_AQ.fst_10K_downsamp, 
                      snp = "SNP", 
                      chrom = "Chromosome", 
                      bp = "POS", 
                      pvalue = "FST", 
                      logTransform = F, 
                      sigLine = 0, 
                      pointSize = 2.5, 
                      ymin = 0, 
                      ymax = 1, 
                      ylabel = "FST", 
                      xlabel = "Chromosome", 
                      title = "Norway wild vs Aquaculture"
) + theme_classic()

Icedsplot<- ggmanHighlight(Icedsplot, highlight = Ice.AQ.fstoutliers.top1K.geneverlap.noloc$SNP)

ggmanLabel(Icedsplot, labelDfm = Ice.AQ.fstoutliers.top1K.geneverlap.noloc, snp = "SNP", label = "gene.id")









NA.AQ.fstoutliers.top1K.geneverlap.noloc %>%  slice_max(FST, n = 5)
Norway.AQ.fstoutliers.top1K.geneverlap.noloc %>%  slice_max(FST, n = 5)

#PCA them 

system("source ~/.zshrc
       conda activate plink
       plink --bfile Ssal.31.220K.mapped.June2023merge --maf 0.01   --keep AQ_NA_clusterfile.tsv --extract NA.AQ.fst.outliers_top1KSNPs.tsv --allow-extra-chr --make-bed --out Ssal.31.220K.mapped.June2023merge.NA.AQ.top1K")

system("source ~/.zshrc
       conda activate plink
       plink --bfile Ssal.31.220K.mapped.June2023merge --maf 0.01  --keep AQ_Norway_clusterfile.tsv --extract Norway.AQ.fst.outliers_top1KSNPs.tsv --allow-extra-chr --make-bed --out Ssal.31.220K.mapped.June2023merge.Norway.AQ.top1K")

NAAQ.pcadapt <- read.pcadapt(input = "Ssal.31.220K.mapped.June2023merge.NA.AQ.bed" , type = "bed")
NAAQ.PCA <- pcadapt(NAAQ.pcadapt, K = 15)

NAAQPCA_meta <- fread("Ssal.31.220K.mapped.June2023merge.NA.AQ.fam") %>%
  mutate(FID = V1, IID = V2) %>% 
  select(FID, IID) %>%  
  inner_join(., AQ_NA_clusterfile) %>% 
  bind_cols(., NAAQ.PCA$scores) %>% 
  rename(PC1 = 4, PC2 = 5)

ggplot() + geom_point(data = NAAQPCA_meta, aes(x = PC1, y = PC2, colour = Region)) + theme_classic()



NorwayAQ.pcadapt <- read.pcadapt(input = "Ssal.31.220K.mapped.June2023merge.Norway.AQ.bed" , type = "bed")
NorwayAQ.PCA <- pcadapt(NorwayAQ.pcadapt, K = 15)

NorwayAQPCA_meta <- fread("Ssal.31.220K.mapped.June2023merge.Norway.AQ.fam") %>%
  mutate(FID = V1, IID = V2) %>% 
  select(FID, IID) %>%  
  inner_join(., AQ_Norway_clusterfile) %>% 
  bind_cols(., NorwayAQ.PCA$scores) %>% 
  rename(PC1 = 4, PC2 = 5)

ggplot() + geom_point(data = NorwayAQPCA_meta, aes(x = PC1, y = PC2, colour = Region)) + theme_classic()

(NorwayAQ.PCA$singular.values ^2) * 100

IceAQ.pcadapt <- read.pcadapt(input = "IceAQ.PolyHigh.bed" , type = "bed")
IceAQ.PCA <- pcadapt(IceAQ.pcadapt, K = 15)

IceAQPCA_meta <- fread("IceAQ.fam") %>%
  mutate(FID = V1, IID = V2) %>% 
  select(FID, IID) %>%  
  inner_join(., IcelandMeta) %>% 
  bind_cols(.,IceAQ.PCA$scores) %>% 
  rename(PC1 = 5, PC2 = 6)

ggplot() + geom_point(data = IceAQPCA_meta, aes(x = PC1, y = PC2, colour = FID_Meaning)) + theme_classic()




