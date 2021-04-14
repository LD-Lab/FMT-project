######  all figure in a paper order
setwd("F:/R file/DaiLab/wangying/FMT_project")
####### change the tree file , and use rbiom::unifrac to replace the phyloseq
## change color from red to orange.

file_name="figure20210415"





######## relative abundance preprocessing -------------
phylo1 = phylo1_get(phylo_all,meta,var1 = "donor2acceptor",
                    var1_sub = c("O2","Y2","O3","Y3","O0","Y0"),var2 = "timeFMT",
                    var2_sub = c(0,1,4,7,14,21,28,35,42,49,56))

phylo2 = phylo1_get(phylo_all,meta,var1 = "donor2acceptor",
                    var1_sub = c("OO","YO","OD","YD"),var2 = "timeFMT",
                    var2_sub = c(-1,0,1,4,7,14,21,28,35,42,49,56))
###merge
phylo3 <- merge_phyloseq(phylo1,phylo2)

otu_tax1 <- otu_tax1_get(phylo_all=phylo3,tax_Name="Genus",meta=meta)



# h=merge(tax,data.frame(Genus=colnames(otu_tax1)),by="Genus")


meta1 <- meta[rownames(otu_tax1),]

tax_melt<- relative_abun_get(otu_tax1,meta1,time="timeFMT",time_type="num",
                             var1="donor2acceptor",var2 = "label_id")

### top based on all data 底层级包含在高层级
tax_melt1=relative_abun_second_get(tax_melt,tax0=tax,Time="timeFMT")


#######remove some of the taxonomy and shorten the names 
tax_level=data.frame(old=levels(tax_melt1$variable))
tax_level$new=sapply(str_split(tax_level$old,"[-gut]*-group"),"[",1)
tax_level[tax_level$new=="G_gut-metagenome","new"]="others"
##### relevels again , base on the specified taxonomy of each taxonomy 每个进化层次只保留一个
tax_melt1$donor2acceptor <- factor(tax_melt1$donor2acceptor,levels = c("O2","Y2","OO","YO","YD",
                                                                       "OD","O3","Y3","O0","Y0"))
tax_melt1$label_id = str_sub(tax_melt1$label_id,2,2)

tax_melt1$variable=factor(tax_melt1$variable,levels = tax_level$old,labels = tax_level$new)

tax_melt1 = merge(tax_melt1,meta[,c("id","cage")],by.x="sampleId",by.y="row.names")

tax_melt1$timeFMT=as.numeric(as.matrix(tax_melt1$timeFMT))



color2taxonomy <- data.frame(taxonomy=levels(tax_melt1$variable),
                             genus_color=c("grey",distinctive_colors(length(unique(tax_melt1$variable))-8),"purple","red","brown","cyan","blue","green","orange"))

tax_melt2 <- subset(tax_melt1, donor2acceptor %in% c("O2","Y2","O3","Y3","YO","OO","Y0","O0"))

#### add the same  taxonomy of each day each mouse before get the mean value , some taxonomy may has two value in a day of a mouse
tax_melt2=tax_melt2 %>% group_by(donor2acceptor,label_id,sampleId,variable,timeFMT,cage) %>%
  summarise(value=sum(value)) %>% ungroup() %>%as.data.frame()


# write.csv(tax_melt2,file =paste(file_name,"relative abundance.csv",sep = "/") )




#########* sf day0 of control -----
p1 = relative_abun_plot(tax_melt1,var1="donor2acceptor",var1_sub = c("O0","Y0"),
                        var2 = "timeFMT",var2_sub = 0,xvalue_order = FALSE,
                        xvalue="label_id",facet="donor2acceptor",tax="Genus",title_name = "day0 control young vs aging mice",
                        color_g = length(unique(tax_melt1$variable)))+ 
  labs(x="mouse id", y="Relative Abundance")

p1




# ggsave(file=paste(file_name,"sf control day0 relative abundance.png",sep = "/"), p1,height = 8, width = 18)
# 
# 
# pdf(file=paste(file_name,"sf control day0 relative abundance.pdf",sep = "/"), height = 8, width = 18) #for pdf
# print(p1)
# dev.off()


######*stack ploting figure1 H relative abundance and figure3 relative abundance ---------
################# for genus relative abundance , old or young
########* old group ABX ---------


p1 = relative_abun_plot(tax_melt1,var1="donor2acceptor",var1_sub = "O3",xvalue_order = FALSE,
                        xvalue="label_id",facet="timeFMT",tax="Genus",title_name = "20 months ABL",
                        color_g = length(unique(tax_melt1$variable)))+ guides(fill="none")+
  labs(x="mouse id", y="Relative Abundance")

p2 = relative_abun_plot(tax_melt1,var1="donor2acceptor",var1_sub = "O2",xvalue_order = FALSE,
                        xvalue="label_id",facet="timeFMT",tax="Genus",title_name = "20 months ABH",
                        color_g = length(unique(tax_melt1$variable)))+ guides(fill="none")+
  labs(x="mouse id", y="Relative Abundance")

p3 = relative_abun_plot(tax_melt1,var1="donor2acceptor",var1_sub = "O2",xvalue_order = FALSE,
                        xvalue="label_id",facet="timeFMT",tax="Genus",title_name = "O2--Genus",
                        color_g = length(unique(tax_melt1$variable))) 

p3
p1/p2

p <- list(p1,p2,p3)


# pdf(file=paste(file_name,"f1  20 mons relative abundance.pdf",sep = "/"), height = 8, width = 18) #for pdf
# print(p)
# dev.off()

########* young group ---------
p1 = relative_abun_plot(tax_melt1,var1="donor2acceptor",var1_sub = "Y3",xvalue_order = FALSE,
                        xvalue="label_id",facet="timeFMT",tax="Genus",title_name = "2 months ABL",
                        color_g = length(unique(tax_melt1$variable)))+ guides(fill="none")+
  labs(x="mouse id", y="Relative Abundance")

p2 = relative_abun_plot(tax_melt1,var1="donor2acceptor",var1_sub = "Y2",xvalue_order = FALSE,
                        xvalue="label_id",facet="timeFMT",tax="Genus",title_name = "2 months ABH",
                        color_g = length(unique(tax_melt1$variable)))+ guides(fill="none")+
  labs(x="mouse id", y="Relative Abundance")

p3 = relative_abun_plot(tax_melt1,var1="donor2acceptor",var1_sub = "Y2",xvalue_order = FALSE,
                        xvalue="label_id",facet="timeFMT",tax="Genus",title_name = "Y2--Genus",
                        color_g = length(unique(tax_melt1$variable))) 

p1/p2
p3
p <- list(p1,p2,p3)


# pdf(file=paste(file_name,"sf6 2 mons relative.pdf",sep = "/"), height = 8, width = 18) #for pdf
# print(p)
# dev.off()

#######* FMT ---------
####### OO, YO 



p1 = relative_abun_plot(tax_melt1,var1="donor2acceptor",var1_sub = "OO",xvalue_order = FALSE,
                        var2 = "cage",var2_sub = "J",xvalue="label_id",facet="timeFMT",tax="Genus",title_name = "20 months aFMT",
                        color_g = length(unique(tax_melt1$variable)))+ guides(fill="none")+
  labs(x="mouse id", y="Relative Abundance")

p2 = relative_abun_plot(tax_melt1,var1="donor2acceptor",var1_sub = "YO",xvalue_order = FALSE,
                        var2 = "cage",var2_sub = "G", xvalue="label_id",facet="timeFMT",tax="Genus",title_name = "20 months hFMT",
                        color_g = length(unique(tax_melt1$variable)))+ guides(fill="none")+
  labs(x="mouse id", y="Relative Abundance")

p3 = relative_abun_plot(tax_melt1,var1="donor2acceptor",var1_sub = "OO",xvalue_order = FALSE,
                        var2 = "cage",var2_sub = "J", xvalue="label_id",facet="timeFMT",tax="Genus",title_name = "20 months aFMT",
                        color_g = length(unique(tax_melt1$variable))) 

p1/p2
p3
p <- list(p1,p2,p3)


# pdf(file=paste(file_name,"sf4  20 mons FMT OO YO relative.pdf",sep = "/"), height = 8, width = 18) #for pdf
# print(p)
# dev.off()



#######* Control -------
####### O0,Y0



p1 = relative_abun_plot(tax_melt1,var1="donor2acceptor",var1_sub = "Y0",xvalue_order = FALSE,
                        xvalue="label_id",facet="timeFMT",tax="Genus",title_name = "2 months ctrl",
                        color_g = length(unique(tax_melt1$variable)))+ guides(fill="none")+
  labs(x="mouse id", y="Relative Abundance")

p2 = relative_abun_plot(tax_melt1,var1="donor2acceptor",var1_sub = "O0",xvalue_order = FALSE,
                        xvalue="label_id",facet="timeFMT",tax="Genus",title_name = "20 months ctrl",
                        color_g = length(unique(tax_melt1$variable)))+ guides(fill="none")+
  labs(x="mouse id", y="Relative Abundance")

p3 = relative_abun_plot(tax_melt1,var1="donor2acceptor",var1_sub = "Y0",xvalue_order = FALSE,
                        xvalue="label_id",facet="timeFMT",tax="Genus",title_name = "ctrl--Genus",
                        color_g = length(unique(tax_melt1$variable))) 

p1/p2

p3
p <- list(p1,p2,p3)


# pdf(file=paste(file_name,"ctrl O0 and Y0 relative abundance.pdf",sep = "/"), height = 8, width = 18) #for pdf
# print(p)
# dev.off()
















###### sf  absolute abundace ------------------
qpcr <- read.csv("QPCR/qpcr_refine.csv")
qpcr$raw_reads=qpcr$res_log10.^10
qpcr1=merge(qpcr,meta[,c("label_id","id")],by.x="sampleId",by.y="row.names")
# qpcr1$timeFMT=paste("D",qpcr1$timeFMT,sep = "_")


# h=demelt(data=qpcr1[,c("label_id","raw_reads","timeFMT")],value="raw_reads",y = "label_id",x = "timeFMT")



# h=dcast(timeFMT~sampleId,qpcr1[,c("sampleId","raw_reads","timeFMT")])
## calculate the absolute value with qpcr
tax_melt2$timeFMT= as.numeric(as.matrix(tax_melt2$timeFMT))

tax_melt3= full_join(tax_melt2,qpcr[,c("sampleId","donor2acceptor","timeFMT","res_log10.")]) %>% 
  na.omit()
 
tax_melt3$ab_value= tax_melt3$value * (tax_melt3$res_log10.)^10

###### old group  ABX
p1 <- line_plot2(tax_melt = tax_melt3,var1 = "donor2acceptor",var1_sub = "O2",color_level=color2taxonomy,
                 Value = "ab_value",
                 group = "variable",time = "timeFMT",tax_num = length(unique(tax_melt3$variable)),mad_all = TRUE )+
  scale_x_continuous(breaks=c(0,1,4,7,14,21,28,35,42,49,56)) +guides(color="none")+
  labs(x="Days",y="absolute Abundance *relative abundance",title = "20 months ABH")+
  # scale_y_log10(limits=c(1,1e+12)) +
  ylim(NA,7e+10)

p1
p2<-   line_plot2(tax_melt = tax_melt3,var1 = "donor2acceptor",var1_sub = "O3",color_level=color2taxonomy,
                   Value = "ab_value",
                   group = "variable",time = "timeFMT",tax_num = length(unique(tax_melt3$variable)),mad_all = TRUE )+
  scale_x_continuous(breaks=c(0,1,4,7,14,21,28,35,42,49,56)) +guides(color="none")+
  labs(x="Days", y="absolute Abundance *relative abundance",title = "20 months ABL")+
  ylim(NA,7e+10)
p2


p3<-   line_plot2(tax_melt = tax_melt3,var1 = "donor2acceptor",var1_sub = "O3",color_level=color2taxonomy,
                  Value = "ab_value",
                  group = "variable",time = "timeFMT",tax_num = length(unique(tax_melt3$variable)),mad_all = TRUE )+
  scale_x_continuous(breaks=c(0,1,4,7,14,21,28,35,42,49,56)) +
  labs(x="Days", y="absolute Abundance *relative abundance",title = "20 months ABL")+
  ylim(NA,7e+10)


p <- list(p1,p2,p3)

p1+p2

pdf(file=paste(file_name,"sf2 absolute abundance single_genus .pdf",sep = "/"), height = 8, width = 8) #for pdf
print(p)
dev.off()




######## Sf figure day21 and donor OO and YO ---------
################# for genus relative abundance


phylo1 = phylo1_get(phylo_all,meta,var1 = "donor2acceptor",
                    var1_sub = c("O2","Y2"),var2 = "timeFMT",
                    var2_sub = c(0,1,4,7,14,21,28,35,42))

phylo2 = phylo1_get(phylo_all,meta,var1 = "donor2acceptor",
                    var1_sub = c("OO","YO","OD","YD"),var2 = "timeFMT",
                    var2_sub = c(-1,21))
###merge
phylo3 <- merge_phyloseq(phylo1,phylo2)




otu_tax1 <- otu_tax1_get(phylo3,tax_Name="Genus",meta)
meta1 <- meta[rownames(otu_tax1),]

tax_melt<- relative_abun_get(otu_tax1,meta1,time="timeFMT",time_type="num",
                             var1="donor2acceptor",var2 = "label_id")

### top based on all data
tax_melt1=relative_abun_second_get(tax_melt,tax0=tax,Time="timeFMT")


levels(tax_melt1$variable)

##### relevels again 
tax_melt1$donor2acceptor <- factor(tax_melt1$donor2acceptor,levels = c("O2","Y2","OO","YO","YD",
                                                                       "OD","O3","Y3","M2","M3"))


tax_melt2 <- subset(tax_melt1, donor2acceptor %in% c("OO","YO","YD","OD"))

# colnames(tax_melt2)
### sum the same taxonomy based on each time points of each mouse id 
tax_melt2<- aggregate(value~sampleId + variable + donor2acceptor + label_id + timeFMT,data = tax_melt2,sum)

meta2 <- subset(meta1,donor2acceptor %in% c("OO","YO","YD","OD"))
### get distance , bray-curtis, to donor
### 
pcoa_otu <- demelt(df= tax_melt2,var1 = "variable",var_group = "sampleId")

all <- pcoa_product(Data = pcoa_otu ,meta=meta1,metric = "bray")

pcoa_distance <- all$Distance  %>% as.matrix()
pcoa_distance <- pcoa_distance[rownames(meta2),rownames(meta2)]
meta2$mix <- paste(meta2$donor2acceptor,meta2$label_id,sep = "_")
rownames(pcoa_distance)<- meta2$mix
colnames(pcoa_distance ) <-meta2$mix 

OO_OD <- pcoa_distance[!is.na(str_match(colnames(pcoa_distance),"OO")),!is.na(str_match(colnames(pcoa_distance),"OD"))]
OO_OD <- OO_OD[order(OO_OD,decreasing = FALSE)]%>% as.matrix()
rownames(OO_OD) <- sapply(str_split(rownames(OO_OD),"_"),"[",2)

YO_YD <- pcoa_distance[!is.na(str_match(colnames(pcoa_distance),"YO")),!is.na(str_match(colnames(pcoa_distance),"YD"))]
YO_YD<- YO_YD[order(YO_YD,decreasing = FALSE)]%>% as.matrix()
rownames(YO_YD) <- sapply(str_split(rownames(YO_YD),"_"),"[",2)


tax_melt2$label_id <- factor(tax_melt2$label_id,levels = c("OD1","YD1",rownames(OO_OD),rownames(YO_YD) ))



p1 <-relative_abun_plot(tax_melt2,var1="donor2acceptor",var1_sub = c("OO","OD"),xvalue_order = FALSE,
                        xvalue="label_id",facet="timeFMT",tax="Genus",
                        title_name = "OO--Genus",color_g = length(unique(tax_melt1$variable)))


p2 <-relative_abun_plot(tax_melt2,var1="donor2acceptor",var1_sub = c("YO","YD"),xvalue_order = FALSE,
                        xvalue="label_id",facet="timeFMT",tax="Genus",
                        title_name = "YO--Genus",color_g = length(unique(tax_melt1$variable)))


p1/p2

p <- list(p1,p2)
# pdf(file=paste(file_name,"sf3 FMT day21 vsdonor relative.pdf",sep = "/"), height = 8, width = 14) #for pdf
# 
# print(p)
# dev.off()
