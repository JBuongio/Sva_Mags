## Goal: Find distribution of Svalbard MAGs across transcriptomic samples

#load libraries
library(vegan)
library(factoextra)
library(ggplot2)
library(mclust)
library(dplyr)
library(randomForest)


mean_cov<-read.csv("Mean_coverage_all.csv", stringsAsFactors = FALSE) #In Excel, change cells with 0's be NA if you're going to transform data. 
str(mean_cov)
head(mean_cov)
nrow(mean_cov) #49 MAGs
ncol(mean_cov) #20 MTs

##going to use metaMDS for NMDS plotting of community abundances of different MAGs
#First, prepare the data for analysis
rownames(mean_cov)<-mean_cov[,1] #make MAG name rowname
mean_cov$Organism.Name<-NULL #remove MAG ID column
sapply(mean_cov, class) # numeric, good
which(rowSums(mean_cov[,1:20])==0) #No rows with all zeros

## create distance matrix of MAGs (not sites, since we are clustring MAGs not sites)
set.seed(2987)
library(scales)
mean_cov_scale<-scale(mean_cov, center=FALSE) #scale the values

hist(mean_cov_scale)
dist<-vegdist(mean_cov_scale, method = "bray", binary=FALSE) #weighted analysis on scaled data
NMDS<- metaMDS(dist, k = 2, trymax = 100, autotransform=TRUE, trace = FALSE)  
stressplot(NMDS) 
plot(NMDS, type="t", display="sites") #not actually sites, but MAGs
NMDS$stress #  0.07656367

#Build dataframe that has NMDS values for x and y plotting, along with sample attributes
NMDS1<-NMDS$points[,1]
NMDS2<-NMDS$points[,2]
MDS<-data.frame(MDS1=NMDS1, MDS2=NMDS2)


head(mean_cov_scale)

BIC<- mclustBIC(mean_cov_scale)
summary(BIC)
plot(BIC)


#Actually run the model based clustering

mclust_clusters <- Mclust(mean_cov_scale, x=BIC) #2 clusters of MAGs
summary(mclust_clusters, parameters = TRUE)
mean_cov$mclust<-as.factor(mclust_clusters$classification)

NMDS_w_clusters<-data.frame(MDS1=NMDS1, MDS2=NMDS2, Cluster = mean_cov$mclust)

# BIC values used for choosing the number of clusters
fviz_mclust(mclust_clusters, 
            "BIC", 
            palette = "npg",
            legend = "bottom")

#add column of MAG ID names to import taxnonomy
tax<-read.csv("MAG_GTDBtax.csv")
NMDS_w_clusters$NCBI_name<-rownames(NMDS_w_clusters)
NMDS_w_clusters<-left_join(NMDS_w_clusters, tax, by = "NCBI_name")
rownames(NMDS_w_clusters)<-NMDS_w_clusters[,4]


NMDS_plot<-ggplot(NMDS_w_clusters, aes(x=MDS1, y=MDS2, label=NCBI_name, col=Cluster)) +
  geom_point(aes(fill=Phylum), size=3, pch=21) + 
  #scale_shape_manual(values=c(21,22,23,24)) +
  #stat_ellipse() +
  geom_text()+
  theme_bw() +
  labs(title="MAG NMDS Bray diss, Mclust cluster")
NMDS_plot


#What samples are most influential in shaping MAG clusters?

randomforest_mod <- randomForest(mean_cov$mclust~.,data = mean_cov)
randomforest_mod$importance
varImpPlot(randomforest_mod) #Abundance of Micrbin229, Iainarch (x2) and Nitrosobin73 MAGs are the most important in determing the site clusters.
randomforest_mod


MAGs_clusters<-data.frame(MAG=rownames(mean_cov), Phylum=NMDS_w_clusters$Phylum, Class=NMDS_w_clusters$Class, Cluster = NMDS_w_clusters$Cluster)
write.csv(MAGs_clusters, "MAGs_clusters.csv")
freq<-as.data.frame(table(MAGs_clusters$Phylum, MAGs_clusters$Cluster))

ggplot(freq) +
  geom_bar(aes(x = Var1, y = Freq, fill = Var2), stat = 'identity') +
  labs(x = "Phylum",
       y = "Count") +
  guides(fill=guide_legend(title="Cluster")) +
  theme_bw() +
  theme(aspect.ratio = 0.9,
        axis.text = element_text(size=8, angle = 45, hjust=1),
        axis.title = element_text(size=10),
        plot.margin = ggplot2::margin(l=0.5, t=0.5, b=0.3, r=0.5, unit = 'cm')
  )
colnames(mean_cov)

#mean coverage values across depths
cov_vals<-aggregate(.~mclust, mean_cov[1:21], sum)
cov_vals_avg<-aggregate(.~mclust, mean_cov[1:21], mean)
write.csv(cov_vals_avg, "cov_vals_clusters_mean.csv")
write.csv(cov_vals, "cov_vals_clusters.csv")

#how about what is transcribed? Let's look at cazymes transcribed across fjords. Grouped together at fjord level and clumped all CAZyme categories together so data was not too sparse.
KF_caz<-read.csv("KF_num_cazymes_trans.csv")
VK_caz<-read.csv("VK_num_cazymes_trans.csv")
KF_caz_tax<-left_join(KF_caz, tax, by="Bin_name")
KF_caz_tax$X<-NULL
colnames(KF_caz_tax)[2]="KF_cazymes_transcribed"

VK_caz_tax<-left_join(VK_caz, tax, by="Bin_name")
VK_caz_tax$X<-NULL
colnames(VK_caz_tax)[2]="VK_cazymes_transcribed"

#make new table that has this cazyme data plus clusters (not worried about the mean coverage info now that we know AB and AC samples are the biggest drivers)
#need to left join?
new<-left_join(VK_caz_tax,KF_caz_tax[2:3], by="NCBI_name")
View(new3)
mean_cov$NCBI_name<-rownames(mean_cov) # add column to use for left_join that is name of MAG (rowname)
new2<-left_join(new, mean_cov[21:22], by="NCBI_name")
new3<-data.frame(new2$NCBI_name, new2$VK_cazymes_transcribed, new2$KF_cazymes_transcribed, new2$mclust)
