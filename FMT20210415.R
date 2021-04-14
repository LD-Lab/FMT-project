######  all figure in a paper order
setwd("F:/R file/DaiLab/wangying/FMT_project")
####### change the tree file , and use rbiom::unifrac to replace the phyloseq
## change color from red to orange.

file_name="figure_20210415"

####### figure names
## f for figure  + age + type , like f1 20m qpcr + metric
## sf for supplementary figure

### legend should be capital of the first alpha


## standardized the variable names

### remove some of the data like sample 1788

##### figure 1 ,2 pcoa baseline use O2,OO,YO
### figure 3 baseline use OO,YO
###  all baseline change name to baseline, in the significant ploting, different from 0910

#--------------
##### import data and preprocessing for 16S (no submit ) -------
meta <- read.csv("metadata_all.csv",header = T,row.names=1,na.strings = "",fill = TRUE,stringsAsFactors = FALSE)
meta$DSS <- as.factor(meta$DSS)
meta$treatment <- as.factor(meta$treatment)

# otu0 <-read.csv("feature_table_all.csv",header = T, row.names = 1)
# otu <- otu_filter(otu0,1/1000,0) 

# write.csv(otu,"feature_table_all_filter.csv")
otu= read.csv("feature_table_all_filter.csv",header = T, row.names = 1)

otu <- otu[,intersect(rownames(meta),colnames(otu))]



tax0 <- read.csv("taxonomy_all.csv",header = T,row.names = 1,na.strings = "",
                 fill = TRUE,stringsAsFactors = FALSE)
tax <- tax0[rownames(otu),]
tax= tax_detail(tax)




otu <- otu[rownames(tax),]
meta[is.na(meta)] <- -100
trefile <- read_tree("tree_all.nwk")
meta$donor2acceptor <- factor(meta$donor2acceptor,
                              levels = c("OD","YD","Y0","O0","Y2","O2","YY","OY","OO","YO","DSS2","DSS18",
                                         "Y3","O3","M0","M2","M3"))
meta= meta[-which(meta$id %in% c(55,56,1788,1784)),]



### tree refine
trefile=ape::drop.tip(trefile,setdiff(trefile$tip.label,rownames(otu)))
otu=otu[intersect(rownames(otu),trefile$tip.label),]

setdiff(trefile$tip.label,rownames(otu))


phylo_all <- phyloseq(otu_table(t(otu), taxa_are_rows = F),
                      sample_data(meta),tax_table(as.matrix(tax)),
                      phy_tree(trefile))

refined_data=list(meta=meta,otu=otu,tax=tax)
#######color data---------
# rgb2hsv(col2rgb("#F8766D"))
# rgb2hsv(col2rgb("orange"))

color_db <- na.omit(data.frame(donor2acceptor=unique(meta$donor2acceptor),value="0",test=1,stringsAsFactors = FALSE))
color_db$donor2acceptor = as.character(color_db$donor2acceptor) 
color_db[color_db$donor2acceptor=="O2","value"] <- hsv(0,1,1)
color_db[color_db$donor2acceptor=="O3","value"] <- hsv(0.9,0.6,1)
color_db[color_db$donor2acceptor=="Y2","value"] <- "#F8766D"
color_db[color_db$donor2acceptor=="Y3","value"] <- hsv(0.01,0.3,0.6)
color_db[color_db$donor2acceptor=="M2","value"] <-  "#00BA38"# green
color_db[color_db$donor2acceptor=="M3","value"] <- hsv(0.60,0.5,0.5)

color_db[color_db$donor2acceptor=="OO","value"] <- "purple" #purple
color_db["OO_DSS",]=c("OO_DSS",hsv(0.1,0.5,0.5),1)

color_db[color_db$donor2acceptor=="YO","value"] <-hsv(0.6,1,1)
color_db["YO_DSS",]=c("YO_DSS",hsv(0.10,1,1),1)  ## "orange"
# "#ADADAD"  gray
color_db[color_db$donor2acceptor=="O0"|color_db$donor2acceptor=="Y0"|color_db$donor2acceptor=="M0"
         ,"value"] <- "black"
color_db[color_db$donor2acceptor=="OD"|color_db$donor2acceptor=="YD","value"] <- "black"
color_db$donor2acceptor <- factor(color_db$donor2acceptor,
                                  levels = c("YD","OD","Y0","O0","Y2","O2","OO","YO",
                                             "Y3","O3","M0","M2","M3","OO_DSS","YO_DSS"))
color_db <- na.omit(color_db)
rownames(color_db) <- color_db$donor2acceptor
color_db <- color_db[order(color_db$donor2acceptor),]
ggplot(color_db,aes(donor2acceptor,test,color=donor2acceptor))+
  geom_point(size=10)+ scale_color_manual(values = color_db$value)


#### f1 qpcr -------------------
##### Old group
qpcr <- read.csv("QPCR/qPCR_data.csv")

qpcr$donor2acceptor <- factor(qpcr$donor2acceptor,
                              levels = c("O0","O2","O3"))
# qpcr=add_mix(qpcr,group = "donor2acceptor",time="timeFMT")

p<-errorbar_line_plot(df=qpcr,var1 = "donor2acceptor",var1_sub = c("O0","O2","O3"),
                      group = "donor2acceptor",time="timeFMT",signifi = TRUE,line="age",
                      sig_com = list(c("O0","O2"),c("O0","O3"))) +
  scale_color_manual(values =color_db[c("O0","O2","O3"),"value"],name="Treatment",
                     breaks = c("O0","O2","O3"),labels=c("Ctrl","ABH","ABL"))+
  labs(title="20 months",y="Bacterial load",x="Day")+
  guides(linetype="none")+
  # scale_linetype_manual(values=c( "solid"))+
  scale_x_continuous(breaks=c(0,1,4,7,14,21,28,35,42))

# h=p$data

# write.csv(h,"qpcr_pvalue.csv")

p


# ggsave(file=paste(file_name,"f1 20months QPCR.png",sep="/"), p,height = 8, width = 8)
# 
# pdf(file=paste(file_name,"f1 20months QPCR.pdf",sep="/"), height = 8, width = 8) #for pdf
# print(p)
# dev.off()







#### figure1 C  alpha diversity -----------
########* O group -------
phylo1 = phylo1_get(phylo_all,meta=refined_data$meta,var1 = "donor2acceptor",
                    var1_sub = c("Y2","O2","O0","Y0","Y3","O3"),var2 = "timeFMT",
                    var2_sub = c(0,1,4,7,14,21,28,35,42,49,56))





p1 <-  alpha_plot2(phylo_ls=phylo1,meta=meta,xlab="timeFMT",color_col = "donor2acceptor",line = "age",
                   var1 = "donor2acceptor",var1_sub = c("O0","O2","O3"),time = "timeFMT",
                   group="donor2acceptor",average = TRUE,sig_com = list(c("O0","O2"),c("O0","O3")),
                   title_name = "",measure = "Observed",scale=0.1,label_size = 5)+
  scale_color_manual(values =color_db[c("O0","O2","O3"),"value"],
                     name="treatment",
                     breaks = c("O0","O2","O3"),
                     labels=c("Ctrl","ABH","ABL"))+
  guides(linetype="none")+
  labs(title = "20 months ",x="Day",y="Observed Species")+
  scale_linetype_manual(values=c( "solid"))+
  scale_x_continuous(breaks=c(0,1,4,7,14,21,28,35,42,49,56))


# h=p1$data
p1

p2 <- alpha_plot2(phylo1,meta=refined_data$meta,xlab="timeFMT",color_col = "donor2acceptor",line = "age",
                  var1 = "donor2acceptor",var1_sub = c("O0","O2","O3"),time = "timeFMT",
                  group="donor2acceptor",average = TRUE,sig_com = list(c("O0","O2"),c("O0","O3")),
                  title_name = "",measure = "Shannon",scale=0.1,label_size = 5)+
  scale_color_manual(values =color_db[c("O0","O2","O3"),"value"],
                     name="treatment",
                     breaks = c("O0","O2","O3"),
                     labels=c("Ctrl","ABH","ABL"))+
  guides(linetype="none")+
  labs(title = "20 months ",x="Day",y="Shannon Index")+
  scale_linetype_manual(values=c( "solid"))+
  scale_x_continuous(breaks=c(0,1,4,7,14,21,28,35,42,49,56))


p1/p2 +plot_layout(guides = "collect")

p = list(p1,p2)


# 
# 
# pdf(file=paste(file_name,"f1 20months alpha.pdf",sep = "/"), height = 8, width = 8) #for pdf
# print(p)
# dev.off()





######* FMT old  group --------
#### figure3 A  alpha diversity, old
phylo1 = phylo1_get(phylo_all,meta,var1 = "donor2acceptor",
                    var1_sub = c("O0","O2","YO","OO"),var2 = "timeFMT",
                    var2_sub = c(0,1,4,7,14,21))

p1 <-  alpha_plot3(phylo_ls=phylo1,meta=refined_data$meta,xlab="timeFMT",color_col = "donor2acceptor",line = "age",
                   var1 = "donor2acceptor",var1_sub = c("O0","O2","YO","OO"),time = "timeFMT",
                   group="donor2acceptor",average = TRUE,
                   sig_com = list(c("O2","YO"),c("O2","OO")),sig_com2=list(c("O0","OO"),c("O0","YO")),
                   title_name = "",measure = "Observed",scale = 0.1)+
  labs(title = "FMT",x="Day",y="Observed Species")+ 
  scale_color_manual(values = color_db[c("O0","O2","OO","YO"),"value"],
                     name="treatment",
                     breaks = c("O0","O2","OO","YO"),
                     labels=c("ctrl","spon","aFMT","hFMT"))+
  scale_x_continuous(breaks=c(0,1,4,7,14,21))

p1
# h1=p1$data

p2 <- alpha_plot3(phylo_ls=phylo1,meta=refined_data$meta,xlab="timeFMT",color_col = "donor2acceptor",line = "age",
                  var1 = "donor2acceptor",var1_sub = c("O0","O2","YO","OO"),time = "timeFMT",
                  group="donor2acceptor",average = TRUE,
                  sig_com = list(c("O2","YO"),c("O2","OO")),sig_com2=list(c("O0","OO"),c("O0","YO")),
                  title_name = "",measure = "Shannon",scale = 0.1)+
  labs(x="Day",title="FMT",y="Shannon Index")+ 
  scale_color_manual(values = color_db[c("O0","O2","OO","YO"),"value"],
                     name="treatment",
                     breaks = c("O0","O2","OO","YO"),
                     labels=c("ctrl","spon","aFMT","hFMT"))+
  scale_x_continuous(breaks=c(0,1,4,7,14,21))

p1/p2 +plot_layout(guides = "collect")

p <- list(p1,p2)

# 
# pdf(file=paste(file_name,"f2 FMT alpha.pdf",sep = "/"), height = 8, width = 8) #for pdf
# 
# print(p)
# dev.off()















########## figure3 D humann2 ##########################################
###########import data
# setwd("F:/R file/DaiLab/wangying/metagenome/humann2")
library(dplyr)
pathway = read.csv("humann2/pathway_refine.csv",row.names = 1)

# pathway1  = read.csv("F:/R file/DaiLab/wangying/metagenome/humann2/pathway_refine.csv",row.names = 1)
# pathway1 = data.frame(prim=rownames(pathway1),split= sapply(str_split(rownames(pathway1),": "),"[",2))

rownames(pathway) <- sapply(str_split(rownames(pathway),": "),"[",2)

colSums(pathway)

pathway_norm = sweep(pathway,2,colSums(pathway),"/")

colSums(pathway_norm)

meta <- read.csv("humann2/meta_genome.csv",
                 row.names = 1,stringsAsFactors = FALSE)

meta[is.na(meta)] <- -100
meta<- subset(meta,rownames(meta)%in% colnames(pathway_norm))


meta$mix=paste(meta$donor2acceptor,meta$timeFMT,sep = "_")

humann2_list=list(meta=meta,pathway_norm=pathway_norm)

####*pcoa bray-curtis -----
meta$mix = paste(meta$donor2acceptor,meta$timeFMT,sep = "_")
res= add_group_mean(counts=pathway_norm,meta,time = "timeFMT",group = "donor2acceptor",
                    mean_second = list(c("OO","YO","O2"),c("YY","OY")),remove_baseline1 = TRUE)

pathway_norm2=res$counts
colSums(pathway_norm)
meta2 = res$meta

metric="bray"
all = pcoa_product(pathway_norm2,meta2,metric=metric)


data_ord=all$data_ord
pcoa=all$pcoa
unique(data_ord$donor2acceptor)


data_ord$donor2acceptor<-factor(data_ord$donor2acceptor,levels=c("YO","OO","O2","OD","YD"))

data_ord = data_ord[-which(rownames(data_ord) %in% c("YY_0~mean","OY_0~mean")),]

data_ord[c("O2_0~mean","OO_0~mean","YO_0~mean"),c("pc1","pc2","pc3")]=data_ord["OD",c("pc1","pc2","pc3")]
### add mean baseline 
p=pcoa_plot_size(data_ord,pcoa,metric=metric,title_name = "PCoA--humann2--pathway",size="timeFMT",
                 color = "donor2acceptor",shape = "donor2acceptor",shape_sub = c(21,22),
                 text="timeFMT")+ scale_fill_manual(values = color_db[c("YO","OO","O2","OD","YD"),"value"])+
  scale_color_manual(values = color_db[c("YO","OO","O2","OD","YD"),"value"])+
  scale_size_continuous(range=c(5,12),breaks = c(0,56),labels = c(0,56))
p


data1 <- data_ord[data_ord$donor2acceptor=="YD",]
data2 <- data_ord[data_ord$donor2acceptor=="OD",]

data_mix <-rbind(data1,data2) 
data_mix$donor2acceptor <-  factor(data_mix$donor2acceptor,
                                   levels = c("YO","OO","O2","OD","YD"))

p_data_mean=p$data 
p_data_mean[,c("pc1","pc2")]=p_data_mean[,c("pc1_mean","pc2_mean")]

new_data=rbind(p$data,p_data_mean)

new_data$donor2acceptor


p1=p+  geom_point(data=p$data,aes(pc1,pc2,color=p$data$donor2acceptor),size=3,alpha=0.5) +
  geom_line(data=new_data,aes(pc1,pc2,group=new_data$Row.names,color=new_data$donor2acceptor),
            size=0.5,alpha=0.5)+
  geom_point(data=data_mix,aes(x=pc1,y=pc2,shape=data_mix[,"donor2acceptor"]),
                 stroke=5,size=10)+ scale_shape_manual(values = c(4,3))+
  labs(fill="Treatment")+
  guides(color="none")
 
  
# stat_ellipse(data=p$data,show.legend = FALSE,aes(x=pc1,y=pc2,fill=donor2acceptor),
#              geom = "polygon",
#              alpha=0.2,level = 0.95) 


p1


# pdf(file=paste(file_name,paste("f3 humann2 PCoA",metric,".pdf"),sep = "/"), height = 8, width = 8) #for pdf
# print(p1)
# dev.off()




####* pathway heatmap ---------
#########metacyc  calculate the level7 and find the difference with wilcox test

metric="bray"
pathway3 = pathway
pathway3$path=rownames(pathway)

p_melt= reshape2::melt(data = pathway3,id.vars="path")
p_melt = merge(p_melt,meta[,c("mix","age")],by.x="variable",by.y="row.names")
p_melt[p_melt$mix %in% c("YY_0","OY_0"),"mix"]="YD"
p_melt[p_melt$mix %in% c("OO_0","YO_0"),"mix"]="OD"
unique(p_melt$mix)
p_melt$mix= factor(p_melt$mix,levels = c("YD","OD","O2_56","OO_56","YO_56"))
p_melt=na.omit(p_melt)

p_melt1=p_melt
# group_list=list(c("OD","YO_56"),c("OD","OO_56"),c("OD","O2_56"))

group_list=list(c("OD","O2_56"))  ## pathway difference based on OD and day 56

for (i in 1:length(group_list)) {
  
  p_melt1=sig_table_get(df=p_melt1,time="path",group = "mix",group_base = group_list[[i]][1],
                        group_change =  group_list[[i]][2])
  p_melt1=padjust_get(p_melt1,time="path")
}


### p_adjust 0.05

pmelt2= na.omit(p_melt1) %>%
  filter(p_adjust <0.05)

###get the DEP pathway for heatmap
path_diff= unique(pmelt2[,c("path","age")])
rownames(path_diff)= path_diff$path

res=add_group_mean(counts = pathway,meta=meta,time = "timeFMT",
                   mean_second = list(c("YO","OO"),c("YY","OY")),
                   remove_baseline1=TRUE,remove_baseline2=TRUE,
                   group = "mix")

pathway1= res$counts
meta1=res$meta
meta1$mix=factor(meta1$mix,levels = c("YD","OD","O2_56","OO_56","YO_56"))

rowSums(pathway1)

meta1 = meta1[order(meta1$mix),]



p=pheatmap_plot(counts = pathway1,meta = meta1 ,DEG = path_diff,nchar=100,cluster_cols = FALSE,
                title_name = "pathway difference based on OD and O2_56 ",Annotation = TRUE,
                scale="row",
                row_colors = color_db[c("Y2","O0","O2","OO","YO"),"value"])

p

# write.csv(pathway1,paste(file_name,"metagenome_pathway.csv",sep = "/"))
# 
# 
# 
# pdf(file=paste(file_name,paste("f3 Humann2 pheatmap level7.pdf"),sep = "/"), height = 8, width = 18) #for pdf
# print(p)
# dev.off()













##########   sf figure  #################################################
##########* lefse O0 vs Yo ---------
lefse_data=read.csv("lefse/lefse_data_ctrl_O0_Y0/R_result.csv",
                    row.names = 1)
lefse_data$group=as.factor(lefse_data$group)
lefse_data$taxonomy=factor(lefse_data$taxonomy,level=lefse_data$taxonomy)

p <- ggplot(lefse_data,aes(x=taxonomy,y=lefse_score,fill=group,color=group))+
  #facet_grid(cols = vars(res3[,"group"]),scale="free")+
  geom_bar(stat = "identity") + 
  scale_color_manual(values = c("black","black")) +
  scale_fill_manual(values=c("black","white"))+
  coord_flip()+ labs(title="lefse O0 vs Yo")+
  theme( panel.grid.minor = element_blank(),
         panel.background = element_blank(),panel.border = element_blank(),
         axis.line = element_line(colour = "black",linetype="solid",size = 1),
         panel.grid.major=element_line(colour=NA),axis.text=element_text(size = 15, face= "bold",family = ""),
         legend.text = element_text(size = 10),legend.title = element_text(size=15),axis.title=element_text(size=15))
p


# pdf(file=paste(file_name,paste("Sf1 lefse_O0-Y0",metric,".pdf"),sep = "/"), height = 8, width = 8) #for pdf
# print(p)
# dev.off()



#####* young group----------
### qpcr 
qpcr <- read.csv("QPCR/qPCR_data.csv")

qpcr$donor2acceptor <- factor(qpcr$donor2acceptor,
                              levels = c("Y0","Y2","Y3"))


p<-errorbar_line_plot(df=qpcr,var1 = "donor2acceptor",var1_sub = c("Y0","Y2","Y3"),
                      group = "donor2acceptor",time="timeFMT",signifi = TRUE,line="age",
                      sig_com = list(c("Y0","Y2"),c("Y0","Y3"))) +
  scale_color_manual(values =color_db[c("Y0","Y2","Y3"),"value"],name="Treatment",
                     breaks = c("Y0","Y2","Y3"),labels=c("Ctrl","ABH","ABL"))+
  labs(title="2 months",y="Bacterial load",x="Day")+
  guides(linetype="none")+
  # scale_linetype_manual(values=c( "solid"))+
  scale_x_continuous(breaks=c(0,1,4,7,14,21,28,35,42))

p


pdf(file=paste(file_name,"sf6 2months QPCR.pdf",sep="/"), height = 8, width = 8) #for pdf
print(p)
dev.off()
########* alpha plot young group-------
phylo1 = phylo1_get(phylo_all,meta=refined_data$meta,var1 = "donor2acceptor",
                    var1_sub = c("Y2","O2","O0","Y0","Y3","O3"),var2 = "timeFMT",
                    var2_sub = c(0,1,4,7,14,21,28,35,42,49,56))


p1 <-  alpha_plot2(phylo1,meta=meta,xlab="timeFMT",color_col = "donor2acceptor",line = "age",
                   var1 = "donor2acceptor",var1_sub = c("Y0","Y2","Y3"),time = "timeFMT",
                   group="donor2acceptor",average = TRUE,sig_com = list(c("Y0","Y2"),c("Y0","Y3")),
                   title_name = "",measure = "Observed",scale=0.1)+
  scale_color_manual(values =color_db[c("Y0","Y2","Y3"),"value"],
                     name="treatment",
                     breaks = c("Y0","Y2","Y3"),
                     labels=c("Ctrl","ABH","ABL"))+
  guides(linetype="none")+
  labs(title = "2 months ",x="Day",y="Observed Species")+
  scale_linetype_manual(values=c( "solid"))+
  scale_x_continuous(breaks=c(0,1,4,7,14,21,28,35,42,49,56))



p1

p2 <- alpha_plot2(phylo1,meta=meta,xlab="timeFMT",color_col = "donor2acceptor",line = "age",
                  var1 = "donor2acceptor",var1_sub = c("Y0","Y2","Y3"),time = "timeFMT",
                  group="donor2acceptor",average = TRUE,sig_com = list(c("Y0","Y2"),c("Y0","Y3")),
                  title_name = "",measure = "Shannon",scale=0.1)+
  scale_color_manual(values =color_db[c("Y0","Y2","Y3"),"value"],
                     name="treatment",
                     breaks = c("Y0","Y2","Y3"),
                     labels=c("Ctrl","ABH","ABL"))+
  guides(linetype="none")+
  labs(title = "2 months ",x="Day",y="Shannon Index")+
  scale_linetype_manual(values=c( "solid"))+
  scale_x_continuous(breaks=c(0,1,4,7,14,21,28,35,42,49,56))


p1/p2 +plot_layout(guides = "collect")

p = list(p1,p2)

# pdf(file=paste(file_name,"sf6 2 mons alpha.pdf",sep = "/"), height = 8, width = 8) #for pdf
# print(p)
# dev.off()





#####* figure weight   with integral areas ---------

weight <- read.csv("weight_DAI/FMT_DSS_weight.csv")
weight$value<- round(as.numeric(weight$value),2)*100
weight <- na.omit(weight)

unique(weight$timeFMT)

weight1 <- subset(weight,timeFMT %in% c(1,4,7,10,13,16,19,24,28,31,35,41,49))

p <- errorbar_line_plot(df=weight1,var1 = "donor2acceptor",var1_sub = c("YO","OO"),
                        group = "donor2acceptor",time="timeFMT",signifi = TRUE,
                        sig_com = list(c("YO","OO")),dis_index = 10)+
  scale_color_manual(values =color_db[c("YO","OO"),"value"] )+labs(y="% weight change")+
  scale_y_continuous(breaks = c(0,-5,-10,-15,-20,-25,-30))+
  scale_x_continuous(position = "top")


p1=p+geom_point(aes(x=p$data$timeFMT-p$data$Dodge),size=3,shape=21,fill="white",stroke=2)
p1

# pdf(file=paste(file_name,"weight.pdf",sep = "/"), height = 8, width = 10) #for pdf
# # 
# print(p1)
# dev.off()





##### figure5 D DAI -------
dai <- read.csv("weight_DAI/FMT_DSS_DAI.csv")
dai$value<- as.numeric(dai$value)

dai <- na.omit(dai)

p <- errorbar_line_plot(df=dai,var1 = "donor2acceptor",var1_sub = c("YO","OO"),
                        group = "donor2acceptor",time="timeFMT",signifi = TRUE,sig_com = list(c("YO","OO")))+
  scale_color_manual(values =color_db[c("YO","OO"),"value"] )+labs(y="DAI")+
  labs(x="Day")

# p$data[p$data$timeFMT!=14,"pval"]<-NA
# p$data$pval <- round(p$data$pval,3)

p1 <- p +
  # geom_text(aes(y=position+2 ,label=pval),size=5) +
     geom_point(aes(x=p$data$timeFMT-p$data$Dodge),size=3,shape=21,fill="white",stroke=2)

p1

# pdf(file=paste(file_name,"DAI.pdf",sep = "/"), height = 8, width = 8) #for pdf
# 
# print(p1)
# dev.off()











