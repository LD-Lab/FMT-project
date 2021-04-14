######  all figure in a paper order
setwd("F:/R file/DaiLab/wangying/FMT_project")
####### change the tree file , and use rbiom::unifrac to replace the phyloseq


file_name="figure_20210415"

##### all the PCoA used the same metric, only for wunifrac and bray 
metric="bray"


####
#### FMT baseline 使用30 个数据，figure1 D 使用control的5个数据




###### figure1. D PCoA and signifi-----------
##### O2 ,with different type of point ploting , FMT all group and distance

phylo1 = phylo1_get(phylo_all,meta,var1 = "donor2acceptor",
                    var1_sub = c("O2","O0"))




all = pcoa_product(Data=phylo1,meta,metric=metric,trefile=trefile)
data_ord=all[[1]]
pcoa=all[[2]]
unique(data_ord$donor2acceptor)

data_ord$donor2acceptor <- as.character(data_ord$donor2acceptor)


data_ord$timeFMT<- as.numeric(data_ord$timeFMT)

######## make the control time to day0, and day0 of O2 and O3 to control
for (i in c("O0")) {
  data_ord[data_ord$donor2acceptor==i,"timeFMT"]=0
}




for (i in c("O2","O3")) {
  ii=paste(substr(i,1,1),"0",sep="")
  data_ord[data_ord$donor2acceptor==i& data_ord$timeFMT==0,"donor2acceptor"]=ii
}



data_ord1 = dd_rowmean2(data1=subset(data_ord,donor2acceptor %in% c("O0")),
                        var1 = "pc1",var2 = "pc2",var3 = "pc3",
                        group1 = "donor2acceptor")

data_ord1=data_ord1[c("O0"),]
data_ord1$donor2acceptor = "O0"
data_ord2=data_ord1[c("O0"),]
data_ord2$donor2acceptor="O2"

data_ord1 = rbind(data_ord1,data_ord2,data_ord[data_ord$timeFMT!=0,])

data_ord1$donor2acceptor = factor(data_ord1$donor2acceptor,levels =c("O0","O2","O3") )
data_ord1$timeFMT <- as.numeric(data_ord1$timeFMT)

### all group 
p <- pcoa_plot_size(data_ord1,pcoa,metric=metric,title_name = "PCoA--all",size="timeFMT",
                    color = "donor2acceptor",shape = "donor2acceptor",shape_sub = c(21,22),
                    var1 = "donor2acceptor",var1_sub = c("O0","O2"),
                    # var1 = "timeFMT",var1_sub = c(0,4,21,42),
                    text="timeFMT")+ 
  # geom_point(data=data_ord1[data_ord1$timeFMT==0,],aes(x=pc1,y=pc2,shape=donor2acceptor))+
  # scale_shape_manual(values=c(21))+
  scale_fill_manual(values = color_db[c("O0","O2"),"value"],
                    name="treatment",breaks = c("O0","O2"),labels=c("Ctrl","ABH")
                    )+
  scale_color_manual(values = color_db[c("O0","O2"),"value"])+
  scale_size_continuous(range=c(1,12),breaks = seq(1,11),
                        labels = c(0,1,4,7,14,21,28,35,42,49,56))

p
p_data_mean=p$data 
p_data_mean[,c("pc1","pc2")]=p_data_mean[,c("pc1_mean","pc2_mean")]

new_data=rbind(p$data,p_data_mean)



p1=p+ geom_point(data=subset(data_ord1,donor2acceptor %in% c("O0")),
                aes(x=pc1,y=pc2,color=donor2acceptor,shape=donor2acceptor),
                stroke=3,size=10)+ 
  scale_shape_manual(values = c(4))+
  geom_point(data=p$data,aes(pc1,pc2),color=color_db[c("O2"),"value"],size=3,alpha=0.5) +
  guides(color="none")+
  geom_line(data=new_data,aes(pc1,pc2,group=new_data$Row.names),color=color_db[c("O2"),"value"],
            size=0.5,alpha=0.5) + 
  labs(size="Day",shape="Treatment")


p1

# pdf(file=paste(file_name,paste("f1 O2 pcoa",metric,".pdf"),sep = "/"), height = 8, width = 8) #for pdf
# print(p1)
# dev.off()


##### O3 ,with different type of point ploting , FMT all group and distance

phylo1 = phylo1_get(phylo_all,meta,var1 = "donor2acceptor",
                    var1_sub = c("O3","O0"))



# metric="wunifrac"
all = pcoa_product(phylo1,meta,metric=metric,trefile=trefile)
data_ord=all[[1]]
pcoa=all[[2]]
unique(data_ord$donor2acceptor)

data_ord$donor2acceptor <- as.character(data_ord$donor2acceptor)


data_ord$timeFMT<- as.numeric(data_ord$timeFMT)

######## make the control time to day0
for (i in c("O0")) {
  data_ord[data_ord$donor2acceptor==i,"timeFMT"]=0
}

for (i in c("O2","O3")) {
  ii=paste(substr(i,1,1),"0",sep="")
  data_ord[data_ord$donor2acceptor==i& data_ord$timeFMT==0,"donor2acceptor"]=ii
}



data_ord1 = dd_rowmean2(data1=subset(data_ord,donor2acceptor %in% c("O0")),
                        var1 = "pc1",var2 = "pc2",var3 = "pc3",
                        group1 = "donor2acceptor")

data_ord1=data_ord1[c("O0"),]
data_ord1$donor2acceptor = rownames(data_ord1)
data_ord2=data_ord1[c("O0"),]
data_ord2$donor2acceptor="O3"


data_ord1 = rbind(data_ord1,data_ord2,data_ord[data_ord$timeFMT!=0,])

data_ord1$donor2acceptor = factor(data_ord1$donor2acceptor,levels =c("O0","O2","O3") )
data_ord1$timeFMT <- as.numeric(data_ord1$timeFMT)

### all group 
p <- pcoa_plot_size(data_ord1,pcoa,metric=metric,title_name = "PCoA--all",size="timeFMT",
                    color = "donor2acceptor",shape = "donor2acceptor",shape_sub = c(21,22),
                    var1 = "donor2acceptor",var1_sub = c("O0","O3"),
                    # var1 = "timeFMT",var1_sub = c(0,4,21,42),
                    text="timeFMT")+ 
  scale_fill_manual(values = color_db[c("O0","O3"),"value"],
                    name="treatment",breaks = c("O0","O3"),labels=c("Ctrl","ABL"))+
  scale_color_manual(values = color_db[c("O0","O3"),"value"])+
  scale_size_continuous(range=c(1,12),breaks = seq(1,11),
                        labels = c(0,1,4,7,14,21,28,35,42,49,56))

p
p_data_mean=p$data 
p_data_mean[,c("pc1","pc2")]=p_data_mean[,c("pc1_mean","pc2_mean")]

new_data=rbind(p$data,p_data_mean)

# data3 <- unique(p$data[p$data$donor2acceptor=="OO",c("pc1_mean","pc2_mean","timeFMT")])
# data4 <- unique(p$data[p$data$donor2acceptor=="Y  O",c("pc1_mean","pc2_mean","timeFMT")])
# data5 <- unique(p$data[p$data$donor2acceptor=="O2",c("pc1_mean","pc2_mean","timeFMT")])

p1=p+ geom_point(data=subset(data_ord1,donor2acceptor %in% c("O0")),
                aes(x=pc1,y=pc2,color=donor2acceptor,shape=donor2acceptor),
                stroke=3,size=10)+ scale_shape_manual(values = c(4,3,2))+
  geom_point(data=p$data,aes(pc1,pc2),color=color_db[c("O3"),"value"],size=3,alpha=0.5) +
  geom_line(data=new_data,aes(pc1,pc2,group=new_data$Row.names),color=color_db[c("O3"),"value"],
            size=0.5,alpha=0.5) + 
  guides(color="none")+
  # scale_fill_manual(values = color_db[c("O0",rep("O3",12)),"value"])+
  labs(size="Day",shape="Treatment")


p1

# pdf(file=paste(file_name,paste("f1 O3 pcoa",metric,".pdf"),sep = "/"), height = 8, width = 8) #for pdf
# print(p1)
# dev.off()



#######* distance and significance  -------
####figure1 D.2  significant  box plot  , baseline: day56 of control
phylo1 = phylo1_get(phylo_all,meta,var1 = "donor2acceptor",
                    var1_sub = c("O3","O2","O0"),var2 = "timeFMT",var2_sub = c(0,56))



# metric="wunifrac"
all = pcoa_product(phylo1,meta,metric=metric,trefile=trefile)
data_ord=all[[1]]
pcoa=all[[2]]
unique(data_ord$donor2acceptor)


pcoa_distance <- all[[3]] %>% as.matrix()

dis_melt = distance_refine(pcoa_distance,meta,mouse_id="label_id",
                           group="donor2acceptor",time="timeFMT",time_pos=3,time_base = "0")

# ##### change the name of O3 day0 and O2 day0 to the same name baseline before run baseline get 使用5个点
dis_melt[dis_melt$variable==0 & dis_melt$donor2acceptor %in% c("O0"),"donor2acceptor"]="O_base"


#### chagne the name of O3 and O2 and O0 day 0 to the same name 
# meta1=meta

meta1=data.frame(phylo_all@sam_data )


meta1$donor2acceptor=as.character(meta1$donor2acceptor)
meta1[meta1$timeFMT==0 & meta1$donor2acceptor %in% c("O0","O2","O3"),"donor2acceptor"]="O_base"


base_distance<- distance_base_product(data=phylo1, meta=meta1,var1="timeFMT",var1_sub="0",
                                      group="donor2acceptor",group_sub=c("O_base"),
                                      metric=metric, trefile=trefile)

dis_melt$mix=paste(str_sub(dis_melt$id,1,2),
                   dis_melt$donor2acceptor,"0",sep = "_")

#### change the value of the day0 
for (i in rownames(base_distance)) {
  dis_melt[dis_melt$mix==i& dis_melt$variable==0,"value"] <- base_distance[i,"value"]
}

dis_melt$donor2acceptor=factor(dis_melt$donor2acceptor,levels = c("O_base","O2","O3"),
                               labels = c("baseline","ABH","ABL"))

###remove the other value of each group
dis_melt=dis_melt[dis_melt$value!=0,]
### with p adjust
p=sig_box_plot(df=dis_melt,var1 = "donor2acceptor",var1_sub = c("baseline","ABH","ABL"),
               # com=list(c("O2","O3"),c("O_base","O2"),c("O_base","O3")),
               var2 = "variable",var2_sub = c(0,56),
               time="donor2acceptor") +
  scale_color_manual(breaks = c("baseline","ABH","ABL"),labels=c("Ctrl","ABH","ABL"),
                     values = color_db[c("O0","O2","O3"),"value"])+
  labs(title = paste(metric,"day56 O2 vs O3"),y="Bray-Curtis distance",color="Treatment")+
  ylim(0.001,1.5) + labs(x="treatment")

p 



# h=p$data
# h_sig=p$layers[[3]][["data"]]
# 
# write.csv(h,paste(file_name,"distance_O2 vs O3.csv",sep = "/"))
# 
# write.csv(h_sig,paste(file_name,"distance_O2 vs O3_pval.csv",sep = "/"))
# 
# 
# ggsave(file=paste(file_name,paste("f1 significance O2 O3",metric,".png"),sep = "/"), p,height = 8, width = 8)
# ####
# pdf(file=paste(file_name,paste("f1 significance O2 O3",metric,".pdf"),sep = "/"), height = 8, width = 8) #for pdf
# print(p)
# dev.off()



##### figure2 PCoA  -----------
##########*  O acceptor comparision, baseline use one donor sample---------

# phylo1 = phylo1_get(phylo_all,meta,var1 = "donor2acceptor",
#                     var1_sub = c("O2","OO","YO","YD","OD"),var2 = "timeFMT",var2_sub = c(-1,0,1,4,7,14,21))

phylo1 = phylo1_get(phylo_all,meta,var1 = "donor2acceptor",
                    var1_sub = c("O2","OO","YO","OY","YY"),var2 = "timeFMT",var2_sub = c(0,21))


# metric="wunifrac"
all = pcoa_distance_product2(data=phylo1,meta,metric=metric,time="timeFMT",remove_baseline = FALSE,
                             trefile=trefile,mean_second = list(c("OO","YO"),c("YY","OY")),
                             group = "donor2acceptor",group_sub=c("OO","YO","YD","OD","O2"))  ## distance to the mean value
data_ord=all$data_ord
pcoa=all$pcoa
unique(data_ord$donor2acceptor)
levels(data_ord$donor2acceptor)


data_ord$donor2acceptor <- factor(data_ord$donor2acceptor,
                                  levels = c("OD","YD","O2","OO","YO"))


data_ord1=subset(data_ord,mix %in% c("O2_21","YO_21","OO_21"))


p <- pcoa_plot_size(data_ord1,pcoa,metric=metric,title_name = "day 0",size="timeFMT",
                    var1 = "timeFMT",var1_sub = c(0,21),
                    color = "donor2acceptor",shape = "donor2acceptor",text="timeFMT")+
  scale_fill_manual(values = color_db[c("O2","OO","YO"),"value"])+
  # # scale_color_gradient(high=gray(0.5),low = gray(0.8))+
  scale_color_manual(values = color_db[c("O2","OO","YO"),"value"])+
  scale_size_continuous(range=c(6,12),breaks = c(4,21))+
  # coord_cartesian(xlim=range(-0.5,0.45),ylim=range(-0.5,0.4)) +
  # coord_cartesian(xlim=range(-0.5,1),ylim=range(-0.5,0.4)) +
  labs(title=paste(metric,"day 21 FMT"))


# p2=p


p



data_ord1=data_ord[c("OD","YD"),]
data_ord1$donor2acceptor = rownames(data_ord1)


data_ord1 = rbind(data_ord1,data_ord[data_ord$timeFMT!=0,])


# data_ord1$donor2acceptor = factor(data_ord1$donor2acceptor,levels =c("O0","O2","O3") )
data_ord1$timeFMT <- as.numeric(data_ord1$timeFMT)

p_data_mean=p$data 
p_data_mean[,c("pc1","pc2")]=p_data_mean[,c("pc1_mean","pc2_mean")]

new_data=rbind(p$data,p_data_mean)

p1=p+ geom_point(data=subset(data_ord1,donor2acceptor %in% c("OD","YD")),
                aes(x=pc1,y=pc2,shape=donor2acceptor),
                stroke=3,size=10)+ scale_shape_manual(values = c(4,3,2) ) +
  geom_point(data=p$data,aes(pc1,pc2,color=donor2acceptor),size=3,alpha=0.5) +
  # annotate("text",x=p$data$pc1,y=p$data$pc2,label=p$data$label_id,size=5)+
  geom_line(data=new_data,aes(pc1,pc2,group=new_data$Row.names,color=donor2acceptor),
            size=0.5,alpha=0.5) +
  guides(color="none")+
  labs(fill="Treatment")
# 



p1 



# pdf(file=paste(file_name,paste("f2.FMT day21 pcoa",metric,".pdf"),sep = "/"), height = 8, width = 8) #for pdf
# print(p1)
# dev.off()

######* distance boxplot ----------
###### day21 to baseline, keep the baseline of each time point, and use OO YO as baseline

phylo1 = phylo1_get(phylo_all,meta,var1 = "donor2acceptor",
                    var1_sub = c("O2","OO","YO","OY","YY"),var2 = "timeFMT",var2_sub = c(0,21))
# metric="wunifrac"
all = pcoa_distance_product2(data=phylo1,meta,metric=metric,time="timeFMT", remove_baseline = FALSE,
                             trefile=trefile,mean_second = list(c("OO","YO"),c("YY","OY")),
                             group = "donor2acceptor",group_sub=c("OO","YO","YD","OD","O2"))


data_ord=all$data_ord
pcoa=all$pcoa
unique(data_ord$donor2acceptor)
levels(data_ord$donor2acceptor)

pcoa_distance <- all$Distance %>% as.matrix()
pcoa_distance =  pcoa_distance[-which(!is.na(str_match(rownames(pcoa_distance),"mean"))),
                               -which(!is.na(str_match(rownames(pcoa_distance),"mean")))] #%>% as.data.frame()




dis_melt = distance_refine3(pcoa_distance,meta=data_ord,mouse_id="label_id",
                           group="donor2acceptor",time="timeFMT",time_pos=3,time_base = "0")

# ##### change the name of O3 day0 and O2 day0 to the same name baseline before run baseline get 
dis_melt[dis_melt$variable==0 & dis_melt$donor2acceptor %in% c("OO","YO"),"donor2acceptor"]="O_base"

#### chagne the name of OO and YO day 0 to the same name 
meta1=refined_data$meta
meta1$donor2acceptor=as.character(meta1$donor2acceptor)
meta1[meta1$timeFMT==0 & meta1$donor2acceptor %in% c("OO","YO"),"donor2acceptor"]="O_base"
meta1$mix=paste(meta1$donor2acceptor,meta1$timeFMT,sep = "_")

base_distance<- distance_base_product(data=phylo1, meta=meta1,var1="timeFMT",var1_sub="0",
                                      metric, trefile=trefile,group="donor2acceptor",group_sub=c("O_base"))

dis_melt$mix=paste(str_sub(dis_melt$id,1,2),
                   dis_melt$donor2acceptor,"0",sep = "_")
#### change the value of the day0 
for (i in rownames(base_distance)) {
  dis_melt[dis_melt$mix==i& dis_melt$variable==0,"value"] <- base_distance[i,"value"]
}


dis_melt[dis_melt$timeFMT==0,"donor2acceptor"]="O_base"
dis_melt=dis_melt[dis_melt$value>0,]

p=sig_box_plot(dis_melt,var1 = "donor2acceptor",var1_sub = c("O_base","YO","O2","OO"),
               # com=list(c("O_base","OO"),c("O_base","YO"),c("OO","YO")),
               var2 = "timeFMT",var2_sub = c("0","21"),time = "donor2acceptor") + 
  labs(title = "day21 vs day0",color="treatment",x="treatment")+
  scale_color_manual(values = color_db[c("O0","O2","OO","YO"),"value"])+ylim(NA,1.5)

p                      


h=p$data
h_sig=p$layers[[3]][['data']]

pdf(file=paste(file_name,paste(metric,"figure2E.signif 21_vs 0.pdf"),sep = "/"), height = 8, width = 8)

print(p)
dev.off()














####* distance lineplot -----------
##########  all day distance to baseline of OO and YO as old donor baseline
######## day21 to baseline, keep the baseline of each time point, and use OO YO as baseline

phylo1 = phylo1_get(phylo_all,refined_data$meta,var1 = "donor2acceptor",
                    var1_sub = c("O2","OO","YO","OY","YY"),var2 = "timeFMT",var2_sub = c(0,1,4,7,14,21,28,35,42,49,56))
metric="bray"
all = pcoa_distance_product2(data=phylo1,meta,metric=metric,time="timeFMT", remove_baseline = FALSE,
                             trefile=trefile,mean_second = list(c("OO","YO"),c("YY","OY")),
                             group = "donor2acceptor",group_sub=c("OO","YO","YD","OD","O2"))


data_ord=all$data_ord
pcoa=all$pcoa
unique(data_ord$donor2acceptor)
levels(data_ord$donor2acceptor)

pcoa_distance <- all$Distance %>% as.matrix()
pcoa_distance =  pcoa_distance[-which(!is.na(str_match(rownames(pcoa_distance),"mean"))),
                               -which(!is.na(str_match(rownames(pcoa_distance),"mean")))] #%>% as.data.frame()



######### distance to old mean donor or young mean donor
dis_melt = rbind(distance_refine3(pcoa_distance,meta=data_ord,mouse_id="label_id",
                            group="donor2acceptor",time="timeFMT",time_pos=3,time_base = "0"),
                 distance_refine3_2(pcoa_distance,meta=data_ord,mouse_id="label_id",
                                  group="donor2acceptor",time="timeFMT",time_pos=3,time_base = "0"))
                 

dis_melt$id=sapply(str_split(dis_melt$id,"_"),"[",1)
dis_melt$cage = str_sub(dis_melt$id,1,1)
dis_melt$mix = paste(dis_melt$donor2acceptor,dis_melt$distanceTo,sep = "_")

dis_melt1=mean_line_get(dis_melt,time="timeFMT",group="mix",value="value")
dis_melt1$timeFMT=as.numeric(dis_melt1$timeFMT)

# dis_melt1$timeFMT=as.factor(dis_melt1$timeFMT)
unique(dis_melt1$timeFMT)
# dis_melt2=dis_melt1[dis_melt1$cage!=k ]

###  old group to old mean donor
p=ggplot(dis_melt1[dis_melt1$distanceTo=="OD",],aes(x=timeFMT,y=value,color=donor2acceptor))+
  geom_line(aes(y=mean,group=donor2acceptor),size=2)+
  geom_point(size=2,position = position_dodge2(width = 0.2))+
  theme( panel.grid.minor = element_blank(),
         panel.background = element_blank(),panel.border = element_blank(),
         axis.line = element_line(colour = "black",linetype="solid",size = 1),
         panel.grid.major=element_line(colour=NA),axis.text=element_text(size = 15, face= "bold",family = ""),
         legend.text = element_text(size = 10),legend.title = element_text(size=15),axis.title=element_text(size=15))+
  scale_color_manual(values = color_db[c("O2","OO","YO"),"value"])+
  labs(title=paste(metric,"FMT distance to OD baseline"),x="Day",
                   y="Bray-Curtis distance",color="Treatment")+
  scale_x_continuous(breaks=c(0,1,4,7,14,21,28,35,42,49,56))

p

#

# pdf(file=paste(file_name,paste("sf4 O2 OO YO distance to old baseline",metric,".pdf"),sep = "/"), height = 8, width = 8)
# 
# print(p)
# dev.off()


###  old group to young mean donor
p=ggplot(dis_melt1[dis_melt1$distanceTo=="YD",],aes(x=timeFMT,y=value,color=donor2acceptor))+
  geom_line(aes(y=mean,group=donor2acceptor),size=2)+
  geom_point(size=2,position = position_dodge2(width = 0.2))+
  theme( panel.grid.minor = element_blank(),
         panel.background = element_blank(),panel.border = element_blank(),
         axis.line = element_line(colour = "black",linetype="solid",size = 1),
         panel.grid.major=element_line(colour=NA),axis.text=element_text(size = 15, face= "bold",family = ""),
         legend.text = element_text(size = 10),legend.title = element_text(size=15),axis.title=element_text(size=15))+
  scale_color_manual(values = color_db[c("O2","OO","YO"),"value"])+
  labs(title=paste(metric,"FMT distance to YD baseline"),x="Day",
       y="Bray Curtis dissimilarity")

p

#

# pdf(file=paste(file_name,paste(metric,"figure O2 OO YO distance to young baseline.pdf"),sep = "/"), height = 8, width = 8)
# 
# print(p)
# dev.off()



##### no minus, figure2  bray dissimilarity------
dis_21=dis_melt1[dis_melt1$timeFMT==21,]
dis_21$mix = paste(dis_21$donor2acceptor,dis_21$distanceTo,sep = "_")

p=base_boxplot(dis_21,x="mix",y="value",fill = "distanceTo",
               com=list(c("O2_OD","O2_YD"),c("OO_OD","OO_YD"),c("YO_OD","YO_YD")))+
  labs(title=paste(metric,"day 21 distance to OD minus distance to YD"),
       x="treatment",y="delta Bray-Curties distance")+
  ylim(NA,1.3)+
  scale_color_manual(values = c("black","grey"),name="Treatment")
  # scale_color_manual(values = color_db[c("O2","OO","YO"),"value"],name="Treatment")


p

h_sig=p$layers[[3]][['data']]

# pdf(file=paste(file_name,paste("f2 day21 distance to OD and YD",metric,".pdf"),sep = "/"), height = 8, width = 8)
# 
# print(p)                                                                                                
# dev.off()
# 
# 
# 
# 



############# FMT PCoA  -------------------
##########*  O receiver comparision  day56——---------

phylo1 = phylo1_get(phylo_all,meta,var1 = "donor2acceptor",
                    var1_sub = c("O2","OO","YO","OY","YY"),var2 = "timeFMT",var2_sub = c(0,56))


# metric="wunifrac"
all = pcoa_distance_product2(data=phylo1,meta,metric=metric,time="timeFMT",remove_baseline = FALSE,
                             trefile=trefile,mean_second = list(c("OO","YO"),c("YY","OY")),
                             group = "donor2acceptor",group_sub=c("OO","YO","YD","OD","O2"))
data_ord=all$data_ord
pcoa=all$pcoa
unique(data_ord$donor2acceptor)
levels(data_ord$donor2acceptor)

# 
# data_ord[data_ord$donor2acceptor=="OO"& data_ord$timeFMT==0,"donor2acceptor"] <-"OD"
# 
# data_ord[data_ord$donor2acceptor=="YO"& data_ord$timeFMT==0,"donor2acceptor"] <-"YD"

data_ord$donor2acceptor <- factor(data_ord$donor2acceptor,
                                  levels = c("OD","YD","O2","OO","YO"))




p <- pcoa_plot_size(data_ord,pcoa,metric=metric,title_name = "day 0",size="timeFMT",
                    var1 = "timeFMT",var1_sub = c(0,4,21,56),
                    # var1 = "timeFMT",var1_sub = c(0,1,4,7,14,21,28,35,42,49,56),
                    color = "donor2acceptor",shape = "donor2acceptor",text="timeFMT")+
  scale_fill_manual(values = color_db[c("OD","YD","O2","OO","YO"),"value"])+
  # # scale_color_gradient(high=gray(0.5),low = gray(0.8))+
  scale_color_manual(values = color_db[c("OD","YD","O2","OO","YO"),"value"])+
  scale_size_continuous(range=c(6,12),breaks = c(4,21))+
  # coord_cartesian(xlim=range(-0.5,0.45),ylim=range(-0.5,0.4)) +
  # coord_cartesian(xlim=range(-0.5,1),ylim=range(-0.5,0.4)) +
  labs(title=paste(metric,"FMT"))

p


data_ord1=data_ord[c("OD","YD"),]
data_ord1$donor2acceptor = rownames(data_ord1)


data_ord1 = rbind(data_ord1,data_ord[data_ord$timeFMT!=0,])

# data_ord1$donor2acceptor = factor(data_ord1$donor2acceptor,levels =c("O0","O2","O3") )
data_ord1$timeFMT <- as.numeric(data_ord1$timeFMT)

p_data_mean=p$data 
p_data_mean[,c("pc1","pc2")]=p_data_mean[,c("pc1_mean","pc2_mean")]

new_data=rbind(p$data,p_data_mean)


p1=p+ geom_point(data=subset(data_ord1,donor2acceptor %in% c("OD","YD")),
                aes(x=pc1,y=pc2,shape=donor2acceptor),
                stroke=3,size=10)+ scale_shape_manual(values = c(4,3,2) ) +
  geom_point(data=p$data,aes(pc1,pc2,color=p$data$donor2acceptor),size=3,alpha=0.5) +
  geom_line(data=new_data,aes(pc1,pc2,group=new_data$Row.names,color=donor2acceptor),
            size=0.5,alpha=0.5)  #+ #for confident area

p1



# pdf(file=paste(file_name,paste(metric,"figure3.FMT56_pcoa.pdf"),sep = "/"), height = 8, width = 8)
# print(p)
# dev.off()


#######*  day 56 distance to baseline---------------

phylo1 = phylo1_getlo_a(phyll,meta,var1 = "donor2acceptor",
                    var1_sub = c("O2","OO","YO","OY","YY"),var2 = "timeFMT",var2_sub = c(0,56))
# metric="wunifrac"
all = pcoa_distance_product2(data=phylo1,meta,metric=metric,time="timeFMT", remove_baseline = FALSE,
                             trefile=trefile,mean_second = list(c("OO","YO"),c("YY","OY")),
                             group = "donor2acceptor",group_sub=c("OO","YO","YD","OD","O2"))


data_ord=all$data_ord
pcoa=all$pcoa
unique(data_ord$donor2acceptor)
levels(data_ord$donor2acceptor)

pcoa_distance <- all$Distance %>% as.matrix()
pcoa_distance =  pcoa_distance[-which(!is.na(str_match(rownames(pcoa_distance),"mean"))),
                               -which(!is.na(str_match(rownames(pcoa_distance),"mean")))]




dis_melt = distance_refine3(pcoa_distance,meta=data_ord,mouse_id="label_id",
                           group="donor2acceptor",time="timeFMT",time_pos=3,time_base = "0")

##### change the name of O3 day0 and O2 day0 to the same name baseline before run baseline get 
dis_melt[dis_melt$variable==0 & dis_melt$donor2acceptor %in% c("OO","YO"),"donor2acceptor"]="O_base"

#### chagne the name of OO and YO day 0 to the same name 
meta1=meta
meta1$donor2acceptor=as.character(meta1$donor2acceptor)
meta1[meta1$timeFMT==0 & meta1$donor2acceptor %in% c("OO","YO"),"donor2acceptor"]="O_base"


base_distance<- distance_base_product(data=phylo1, meta=meta1,var1="timeFMT",var1_sub="0",metric = metric,
                                       trefile=trefile,group="donor2acceptor",group_sub=c("O_base"))

dis_melt$mix=paste(str_sub(dis_melt$id,1,2),
                   dis_melt$donor2acceptor,"0",sep = "_")
#### change the value of the day0 
for (i in rownames(base_distance)) {
  dis_melt[dis_melt$mix==i& dis_melt$variable==0,"value"] <- base_distance[i,"value"]
}


dis_melt$donor2acceptor
dis_melt[dis_melt$timeFMT==0,"donor2acceptor"]="O_base"
dis_melt=dis_melt[dis_melt$value>0,]

p=sig_box_plot(dis_melt,var1 = "donor2acceptor",var1_sub = c("O_base","OO","YO","O2"),
               # com=list(c("0","56")),
               var2 = "timeFMT",var2_sub = c("0","56"),time="donor2acceptor") + 
  labs(title = paste(metric,"day56 vs day0"),color="Day",x="Day")+
  scale_color_manual(values = color_db[c("O0","O2","OO","YO"),"value"])+ylim(NA,1)

p


# pdf(file=paste(file_name,paste(metric,"figure3.signif 56_vs 0.pdf"),sep = "/"), height = 8, width = 8)
# 
# print(p)
# dev.off()


##### figure4  YO FMT with or not DSS--------------
metric="bray"

meta1= sub_select(meta,var1 = "cage",
                  var1_sub = c("H","G","YD","OD") )

meta1[meta1$cage=="G","mix_time"]=meta1[meta1$cage=="G","timeFMT"]
meta1[meta1$cage=="H","mix_time"]=meta1[meta1$cage=="H","raw_time"]
meta1[meta1$cage=="OD"|meta1$cage=="YD","mix_time"]=0

### remove some of the time points,## remove day0 exclude the donor
sort(unique(as.numeric(meta1$mix_time)))
meta1 =rbind(subset(meta1,mix_time %in% c(21,28,29,35,36,39,42,46,49,53,56,60,67,74)),
             meta1[meta1$cage=="OD"|meta1$cage=="YD",])




phylo1 = phylo1_get(phylo_all,meta1,var1 = "cage",
                    var1_sub = c("H","G","YD","OD"),)



all = pcoa_product(Data=phylo1,meta1,metric=metric,trefile=trefile)

data_ord=all[[1]]
pcoa=all[[2]]
unique(data_ord$donor2acceptor)




# data_ord[data_ord$cage=="G","mix_time"]=data_ord[data_ord$cage=="G","timeFMT"]
# data_ord[data_ord$cage=="H","mix_time"]=data_ord[data_ord$cage=="H","raw_time"]
# data_ord[data_ord$cage=="OD"|data_ord$cage=="YD","mix_time"]=0




# unique(data_ord$timeDSS)


unique(data_ord$mix_time)

unique(data_ord$cage)
data_ord$cage = factor(data_ord$cage,levels = c("G","H","OD","YD"),labels = c("YO_G_FMT","YO_H_FMT_DSS","OD","YD"))
data_ord=data_ord[data_ord$mix_time!=25,]

data_ord$donor2acceptor
length(unique(data_ord$timeDSS))
p<- pcoa_plot_size(data_ord,pcoa,metric=metric,title_name = "PCoA-FMT-with or not DSS",size="mix_time" ,
                   # var1 = "donor2acceptor",var1_sub = c("OO","YO"),
                   color = "cage",shape = "donor2acceptor",shape_sub = c(21,22),
                   text="mix_time")+ scale_fill_manual(values = color_db[c("YO","YO_DSS","OD","YD"),"value"])+
  scale_color_manual(values = color_db[c("YO","YO_DSS","OD","YD"),"value"])+
  scale_size_continuous(range=c(1,12),breaks = c(seq(1,6)),labels = c(0,1,4,7,11,14))


p_data_mean=p$data 
p_data_mean[,c("pc1","pc2")]=p_data_mean[,c("pc1_mean","pc2_mean")]

new_data=rbind(p$data,p_data_mean)

data_donor = subset(p$data,cage %in% c("OD","YD"))

p
p1=p+ 
  geom_point(data=data_donor,aes(shape=donor2acceptor),
              stroke=3,size=10)+ scale_shape_manual(values = c(4,3))# +
  # geom_point(data=p$data,aes(pc1,pc2,color=p$data$cage),size=3,alpha=0.5) +
  # geom_line(data=new_data,aes(pc1,pc2,group=new_data$Row.names,color=cage),
  #           size=0.5,alpha=0.5)

      

p1


# pdf(file=paste(file_name,paste("f4.YO_FMT+DSS",metric,".pdf"),sep = "/"), height = 8, width = 8)
# 
# print(p1)
# dev.off()


#######*    OO FMT with or not DSS--------
###### DSS after FMT and FMT only to plot PCOA ,OO group


phylo1 = phylo1_get(phylo_all,meta,var1 = "cage",
                    var1_sub = c("J","K","YD","OD"))

# metric="wunifrac"
all = pcoa_product(Data=phylo1,meta,metric=metric,trefile=trefile)

data_ord=all$data_ord
pcoa=all$pcoa
unique(data_ord$donor2acceptor)


data_ord$mix_time = 0

data_ord[data_ord$cage=="J","mix_time"]=data_ord[data_ord$cage=="J","timeFMT"]
data_ord[data_ord$cage=="K","mix_time"]=data_ord[data_ord$cage=="K","raw_time"]
data_ord[data_ord$cage=="OD"|data_ord$cage=="YD","mix_time"]=0

unique(data_ord$mix_time)

data_ord$cage = factor(data_ord$cage,levels = c("J","K","OD","YD"),labels = c("OO_J_FMT","OO_K_FMT_DSS","OD","YD"))

unique(data_ord$cage)

data_ord=data_ord[data_ord$mix_time!=25,]


p<- pcoa_plot_size(data_ord,pcoa,metric=metric,title_name = "PCoA-FMT with or not DSS",size="mix_time" ,
                   # var1 = "donor2acceptor",var1_sub = c("OO","YO"),
                   color = "cage",shape = "donor2acceptor",shape_sub = c(21,22),
                   text="mix_time")+ scale_fill_manual(values = color_db[c("OO","OO_DSS","OD","YD"),"value"])+
  scale_color_manual(values = color_db[c("OO","OO_DSS","OD","YD"),"value"])+
  scale_size_continuous(range=c(1,12),breaks = c(seq(1,6)),labels = c(0,1,4,7,11,14))






p_data_mean=p$data 
p_data_mean[,c("pc1","pc2")]=p_data_mean[,c("pc1_mean","pc2_mean")]

new_data=rbind(p$data,p_data_mean)

data_donor = subset(p$data,cage %in% c("OD","YD"))


p1=p+ geom_point(data=data_donor,aes(shape=donor2acceptor),
                stroke=3,size=10)+ scale_shape_manual(values = c(4,3)) #+
  # geom_point(data=p$data,aes(pc1,pc2,color=p$data$cage),size=3,alpha=0.3) +
  # geom_line(data=new_data,aes(pc1,pc2,group=new_data$Row.names,color=cage),
  #           size=0.5,alpha=0.5)
p1



# pdf(file=paste(file_name,paste("f4 OO_FMT+DSS", metric,".pdf"),sep = "/"), height = 8, width = 8)
# 
# print(p1)
# dev.off()



#######* distance to donor boxplot -------------- 
# phylo1=phylo1_get(phylo_all,meta,var1 = "cage",
#                   var1_sub = c("J","K","H","G"))
metric="bray"

phylo1=phylo1_get(phylo_all,meta,var1 = "donor2acceptor",var1_sub = c("YY","OY","YO","OO"))

otu_d1= otu_table(phylo1) %>% as.data.frame()
meta_d1 = meta[rownames(otu_d1),]
meta_d1$mix_time = 0


####  make  a new time columns with timeFMT and timeDSS
meta_d1[meta_d1$FMT==1,"mix_time"]=meta_d1[meta_d1$FMT==1,"timeFMT"]
meta_d1[meta_d1$DSS==1,"mix_time"]=meta_d1[meta_d1$DSS==1,"timeDSS"]


##### select day 49
meta_d2= subset(meta_d1,mix_time %in% c(0,-25,49))
meta_d2[meta_d2$mix_time==-25,"mix_time"]=0
unique(meta_d2$mix_time)
otu_d2=otu_d1[rownames(meta_d2),]


phylo1=phyloseq(otu_table(otu_d2, taxa_are_rows = F),
                sample_data(meta_d2),tax_table(as.matrix(tax0)),
                phy_tree(trefile))

### distance calculation
# metric="unifrac"
all = pcoa_distance_product2(data=phylo1,meta_d2,metric=metric,time="mix_time",
                             trefile=trefile,mean_second = list(c("OO","YO"),c("YY","OY")),
                             group = "donor2acceptor",group_sub=c("OO","YO","YD","OD"))
data_ord=all$data_ord
pcoa=all$pcoa
unique(data_ord$donor2acceptor)



data1 <- data_ord[data_ord$donor2acceptor=="YD",]
data2 <- data_ord[data_ord$donor2acceptor=="OD",]

data_mix <-rbind(data1,data2)
data_mix$donor2acceptor <-  factor(data_mix$donor2acceptor,
                                   levels = c("OD","YD","O2","OO","YO"))



distance=all$Distance %>% as.data.frame()
distance=distance[rownames(data_ord),rownames(data_ord)]


distance1 <- distance[rownames(data_mix),rownames(data_ord)]

# distance1 <- distance1[,rownames(data_ord)]

colnames(distance1) == rownames(data_ord)
colnames(distance1) <- paste(data_ord$donor2acceptor,data_ord$label_id,data_ord$mix_time,sep = "_")

distance1$group <- rownames(distance1)


dis_melt <- reshape2::melt(distance1,id.var="group")
dis_melt$donor2acceptor <- sapply(str_split(dis_melt$variable,"_"),"[",1)
dis_melt = subset(dis_melt , donor2acceptor %in% c("YO","OO"))

dis_melt$cage=  str_sub(dis_melt$variable,4,4) %>% 
  factor(levels=c("K","J","G","H"),
         labels = c("OO_K_FMT","OO_J_FMT_DSS","YO_G_FMT","YO_H_FMT_DSS"))

# dis_melt$cage = factor(dis_melt$cage,levels=c("K","J","G","H"),
#                        labels = c("OO_K_FMT","OO_J_FMT_DSS","YO_G_FMT","YO_H_FMT_DSS"))


#### compare FMT and FMT+DSS
base_boxplot(df=dis_melt,x="cage",y="value",fill = "donor2acceptor",facet = "group",

             com = list(c("OO_K_FMT","OO_J_FMT_DSS"),c("YO_G_FMT","YO_H_FMT_DSS")))+
  scale_color_manual(values=color_db[c("OO","YO"),"value"])+ylim(NA,0.85)+
  labs(y="Distance to YD or OD",title=paste(metric,"distance from day 49 of FMT or FMT+DSS to average donor baseline"))




# h_sig=p$layers[[3]][['data']]
# 
# pdf(file=paste(file_name,paste("f4.FMT+DSS signif",metric,".pdf"),sep = "/"), height = 8, width = 8)
# 
# print(p)
# dev.off()




#################### supplementary ---------------------
#####* O0 and Y0 ---------

phylo1 = phylo1_get(phylo_all,meta,var1 = "donor2acceptor",
                    var1_sub = c("Y0","O0"))




all = pcoa_product(Data=phylo1,meta,metric=metric,trefile=trefile)
data_ord=all[[1]]
pcoa=all[[2]]
unique(data_ord$donor2acceptor)

data_ord$donor2acceptor <- as.character(data_ord$donor2acceptor)


data_ord$timeFMT<- as.numeric(data_ord$timeFMT)


### all group 
p <- pcoa_plot_size(data_ord,pcoa,metric=metric,title_name = "PCoA--all",size="timeFMT",
                    color = "donor2acceptor",shape = "donor2acceptor",shape_sub = c(10,22),
                    var1 = "donor2acceptor",var1_sub = c("O0","Y0"),
                    # var1 = "timeFMT",var1_sub = c(0,4,21,42),
                    text="timeFMT")+ 
  # geom_point(data=data_ord1[data_ord1$timeFMT==0,],aes(x=pc1,y=pc2,shape=donor2acceptor))+
  # scale_shape_manual(values=c(21))+
  scale_fill_manual(values = c("white","black"),
                    name="treatment",breaks = c("O0","Y0"),labels=c("O0","Y0"))+
  scale_color_manual(values = c("black","black"))+
  scale_size_continuous(range=c(1,12),breaks = seq(1,11),
                        labels = c(0,1,4,7,14,21,28,35,42,49,56))+xlim(-0.4,0.35)+ylim(-0.4,0.35)

p




# 
# pdf(file=paste(file_name,paste("sf1.O0 vs Y0 -pcoa-",metric,".pdf"),sep = "/"), height = 8, width = 8) #for pdf
# print(p)
# dev.off()


##################### supplementary figure ###############
###### sf6
#####* Y2 ,with different type of point ploting , FMT all group and distance-------
phylo1 = phylo1_get(phylo_all,meta,var1 = "donor2acceptor",
                    var1_sub = c("Y2","Y0"))




all = pcoa_product(Data=phylo1,meta,metric=metric,trefile=trefile)
data_ord=all[[1]]
pcoa=all[[2]]
unique(data_ord$donor2acceptor)

data_ord$donor2acceptor <- as.character(data_ord$donor2acceptor)


data_ord$timeFMT<- as.numeric(data_ord$timeFMT)

######## make the control time to day0, and day0 of Y2  to control
for (i in c("Y0")) {
  data_ord[data_ord$donor2acceptor==i,"timeFMT"]=0
}




for (i in c("Y2","Y3")) {
  ii=paste(substr(i,1,1),"0",sep="")
  data_ord[data_ord$donor2acceptor==i& data_ord$timeFMT==0,"donor2acceptor"]=ii
}



data_ord1 = dd_rowmean2(data1=subset(data_ord,donor2acceptor %in% c("Y0")),
                        var1 = "pc1",var2 = "pc2",var3 = "pc3",
                        group1 = "donor2acceptor")

data_ord1=data_ord1[c("Y0"),]
data_ord1$donor2acceptor = "Y0"
data_ord2=data_ord1[c("Y0"),]
data_ord2$donor2acceptor="Y2"

data_ord1 = rbind(data_ord1,data_ord2,data_ord[data_ord$timeFMT!=0,])

data_ord1$donor2acceptor = factor(data_ord1$donor2acceptor,levels =c("Y0","Y2","Y3") )
data_ord1$timeFMT <- as.numeric(data_ord1$timeFMT)

### all group 
p <- pcoa_plot_size(data_ord1,pcoa,metric=metric,title_name = "PCoA--all",size="timeFMT",
                    color = "donor2acceptor",shape = "donor2acceptor",shape_sub = c(21,22),
                    var1 = "donor2acceptor",var1_sub = c("Y0","Y2"),
                    # var1 = "timeFMT",var1_sub = c(0,4,21,42),
                    text="timeFMT")+ 
  # geom_point(data=data_ord1[data_ord1$timeFMT==0,],aes(x=pc1,y=pc2,shape=donor2acceptor))+
  # scale_shape_manual(values=c(21))+
  scale_fill_manual(values = color_db[c("Y0","Y2"),"value"],
                    name="treatment",breaks = c("Y0","Y2"),labels=c("Ctrl","ABH")
  )+
  scale_color_manual(values = color_db[c("Y0","Y2"),"value"])+
  scale_size_continuous(range=c(1,12),breaks = seq(1,11),
                        labels = c(0,1,4,7,14,21,28,35,42,49,56))

p
p_data_mean=p$data 
p_data_mean[,c("pc1","pc2")]=p_data_mean[,c("pc1_mean","pc2_mean")]

new_data=rbind(p$data,p_data_mean)



p1=p+ geom_point(data=subset(data_ord1,donor2acceptor %in% c("Y0")),
                 aes(x=pc1,y=pc2,color=donor2acceptor,shape=donor2acceptor),
                 stroke=3,size=10)+ 
  scale_shape_manual(values = c(4))+
  geom_point(data=p$data,aes(pc1,pc2),color=color_db[c("Y2"),"value"],size=3,alpha=0.5) +
  guides(color="none")+
  geom_line(data=new_data,aes(pc1,pc2,group=new_data$Row.names),color=color_db[c("Y2"),"value"],
            size=0.5,alpha=0.5) + 
  labs(size="Day",shape="Treatment")


p1

# pdf(file=paste(file_name,paste("sf6 Y2 pcoa",metric,".pdf"),sep = "/"), height = 8, width = 8) #for pdf
# print(p1)
# dev.off()



#####* Y3 ,with different type of point ploting , FMT all group and distance-------


phylo1 = phylo1_get(phylo_all,meta,var1 = "donor2acceptor",
                    var1_sub = c("Y3","Y0"))




all = pcoa_product(Data=phylo1,meta,metric=metric,trefile=trefile)
data_ord=all[[1]]
pcoa=all[[2]]
unique(data_ord$donor2acceptor)

data_ord$donor2acceptor <- as.character(data_ord$donor2acceptor)


data_ord$timeFMT<- as.numeric(data_ord$timeFMT)

######## make the control time to day0, and day0 of Y2  to control
for (i in c("Y0")) {
  data_ord[data_ord$donor2acceptor==i,"timeFMT"]=0
}




for (i in c("Y2","Y3")) {
  ii=paste(substr(i,1,1),"0",sep="")
  data_ord[data_ord$donor2acceptor==i& data_ord$timeFMT==0,"donor2acceptor"]=ii
}



data_ord1 = dd_rowmean2(data1=subset(data_ord,donor2acceptor %in% c("Y0")),
                        var1 = "pc1",var2 = "pc2",var3 = "pc3",
                        group1 = "donor2acceptor")

data_ord1=data_ord1[c("Y0"),]
data_ord1$donor2acceptor = "Y0"
data_ord2=data_ord1[c("Y0"),]
data_ord2$donor2acceptor="Y3"

data_ord1 = rbind(data_ord1,data_ord2,data_ord[data_ord$timeFMT!=0,])

data_ord1$donor2acceptor = factor(data_ord1$donor2acceptor,levels =c("Y0","Y2","Y3") )
data_ord1$timeFMT <- as.numeric(data_ord1$timeFMT)

### all group 
p <- pcoa_plot_size(data_ord1,pcoa,metric=metric,title_name = "PCoA--all",size="timeFMT",
                    color = "donor2acceptor",shape = "donor2acceptor",shape_sub = c(21,22),
                    var1 = "donor2acceptor",var1_sub = c("Y0","Y3"),
                    # var1 = "timeFMT",var1_sub = c(0,4,21,42),
                    text="timeFMT")+ 
  # geom_point(data=data_ord1[data_ord1$timeFMT==0,],aes(x=pc1,y=pc2,shape=donor2acceptor))+
  # scale_shape_manual(values=c(21))+
  scale_fill_manual(values = color_db[c("Y0","Y3"),"value"],
                    name="treatment",breaks = c("Y0","Y3"),labels=c("Ctrl","ABL")
  )+
  scale_color_manual(values = color_db[c("Y0","Y3"),"value"])+
  scale_size_continuous(range=c(1,12),breaks = seq(1,11),
                        labels = c(0,1,4,7,14,21,28,35,42,49,56))

p
p_data_mean=p$data 
p_data_mean[,c("pc1","pc2")]=p_data_mean[,c("pc1_mean","pc2_mean")]

new_data=rbind(p$data,p_data_mean)



p1=p+ geom_point(data=subset(data_ord1,donor2acceptor %in% c("Y0")),
                 aes(x=pc1,y=pc2,color=donor2acceptor,shape=donor2acceptor),
                 stroke=3,size=10)+ 
  scale_shape_manual(values = c(4))+
  geom_point(data=p$data,aes(pc1,pc2),color=color_db[c("Y3"),"value"],size=3,alpha=0.5) +
  guides(color="none")+
  geom_line(data=new_data,aes(pc1,pc2,group=new_data$Row.names),color=color_db[c("Y3"),"value"],
            size=0.5,alpha=0.5) + 
  labs(size="Day",shape="Treatment")


p1

# pdf(file=paste(file_name,paste("sf6 Y3 pcoa",metric,".pdf"),sep = "/"), height = 8, width = 8) #for pdf
# print(p1)
# dev.off()



#######* distance and significance ---------
#baseline: day42 of control
phylo1 = phylo1_get(phylo_all,meta,var1 = "donor2acceptor",
                    var1_sub = c("Y3","Y2","Y0"),var2 = "timeFMT",var2_sub = c(0,42))

# metric="wunifrac"
all = pcoa_product(phylo1,meta,metric=metric,trefile=trefile)
data_ord=all[[1]]
pcoa=all[[2]]
unique(data_ord$donor2acceptor)


pcoa_distance <- all[[3]] %>% as.matrix()

dis_melt = distance_refine(pcoa_distance,meta,mouse_id="label_id",
                           group="donor2acceptor",time="timeFMT",time_pos=3,time_base = "0")

# ##### change the name of O3 day0 and O2 day0 to the same name baseline before run baseline get 
dis_melt[dis_melt$variable==0 & dis_melt$donor2acceptor %in% c("Y0"),"donor2acceptor"]="Y_base"


#### chagne the name of O3 and O2 and O0 day 0 to the same name 
# meta1=meta

meta1=data.frame(phylo_all@sam_data )


meta1$donor2acceptor=as.character(meta1$donor2acceptor)
meta1[meta1$timeFMT==0 & meta1$donor2acceptor %in% c("Y0","Y2","Y3"),"donor2acceptor"]="Y_base"


base_distance<- distance_base_product(data=phylo1, meta=meta1,var1="timeFMT",var1_sub="0",
                                      group="donor2acceptor",group_sub=c("Y_base"),
                                      metric=metric, trefile=trefile)

dis_melt$mix=paste(str_sub(dis_melt$id,1,2),
                   dis_melt$donor2acceptor,"0",sep = "_")

#### change the value of the day0 
for (i in rownames(base_distance)) {
  dis_melt[dis_melt$mix==i& dis_melt$variable==0,"value"] <- base_distance[i,"value"]
}

dis_melt$donor2acceptor=factor(dis_melt$donor2acceptor,levels = c("Y_base","Y2","Y3"),
                               labels = c("baseline","ABH","ABL"))

p=sig_box_plot(dis_melt,var1 = "donor2acceptor",var1_sub = c("baseline","ABH","ABL"),
               # com=list(c("O2","O3"),c("O_base","O2"),c("O_base","O3")),
               var2 = "variable",var2_sub = c(0,42),
               time="donor2acceptor") +
  scale_color_manual(breaks = c("baseline","ABH","ABL"),labels=c("Ctrl","ABH","ABL"),
                     values = color_db[c("Y0","Y2","Y3"),"value"])+
  labs(title = paste(metric,"day42 Y2 vs Y3"),y="Bray-Curtis distance",color="Treatment")+
  ylim(0.001,1.2) + labs(x="treatment")





p 

# pdf(file=paste(file_name,paste("sf6 Y2 Y3 signi",metric,".pdf"),sep = "/"), height = 8, width = 8) #for pdf
# print(p)
# dev.off()




##########################################  other
#####* OY FMT with or not DSS ---------




phylo1 = phylo1_get(phylo_all,meta,var1 = "cage",
                    var1_sub = c("W","V","YD","OD"))


# metric="bray"
all = pcoa_product(Data=phylo1,meta,metric=metric,trefile=trefile)

data_ord=all[[1]]
pcoa=all[[2]]
unique(data_ord$donor2acceptor)

data_ord[data_ord$cage=="W","mix_time"]=data_ord[data_ord$cage=="W","timeFMT"]
data_ord[data_ord$cage=="V","mix_time"]=data_ord[data_ord$cage=="V","raw_time"]


data_ord[data_ord$cage=="OD"|data_ord$cage=="YD","mix_time"]=0

# unique(data_ord$timeDSS)


unique(data_ord$mix_time)

unique(data_ord$cage)
data_ord$cage = factor(data_ord$cage,levels = c("W","V","OD","YD"),labels = c("OY_W_FMT","OY_V_FMT_DSS","OD","YD"))
data_ord=data_ord[data_ord$mix_time!=25,]

data_ord$donor2acceptor
length(unique(data_ord$timeDSS))
p<- pcoa_plot_size(data_ord,pcoa,metric=metric,title_name = "OY  PCoA-FMT-with or not DSS",size="mix_time" ,
                   # var1 = "donor2acceptor",var1_sub = c("OO","YO"),
                   color = "cage",shape = "donor2acceptor",shape_sub = c(21,22),
                   text="mix_time")+ scale_fill_manual(values = color_db[c("YO","YO_DSS","OD","YD"),"value"])+
  scale_color_manual(values = color_db[c("YO","YO_DSS","OD","YD"),"value"])+
  scale_size_continuous(range=c(1,12),breaks = c(seq(1,6)),labels = c(0,1,4,7,11,14))


p_data_mean=p$data 
p_data_mean[,c("pc1","pc2")]=p_data_mean[,c("pc1_mean","pc2_mean")]

new_data=rbind(p$data,p_data_mean)

data_donor = subset(p$data,cage %in% c("OD","YD"))

p
p1=p+ 
  geom_point(data=data_donor,aes(shape=donor2acceptor),
             stroke=3,size=10)+ scale_shape_manual(values = c(4,3))
# geom_point(data=p$data,aes(pc1,pc2,color=p$data$cage),size=3,alpha=0.5) +
# geom_line(data=new_data,aes(pc1,pc2,group=new_data$Row.names,color=cage),
#           size=0.5,alpha=0.5)



p1
# pdf(file=paste(file_name,paste("other OY FMT-DSS pcoa",metric,".pdf"),sep = "/"), height = 8, width = 8) #for pdf
# print(p1)
# dev.off()



#####* YY FMT with or not DSS -------

phylo1 = phylo1_get(phylo_all,meta,var1 = "cage",
                    var1_sub = c("U","S","YD","OD"))


# metric="bray"
all = pcoa_product(Data=phylo1,meta,metric=metric,trefile=trefile)

data_ord=all[[1]]
pcoa=all[[2]]
unique(data_ord$donor2acceptor)

data_ord[data_ord$cage=="S","mix_time"]=data_ord[data_ord$cage=="S","timeFMT"]
data_ord[data_ord$cage=="U","mix_time"]=data_ord[data_ord$cage=="U","raw_time"]


data_ord[data_ord$cage=="OD"|data_ord$cage=="YD","mix_time"]=0

# unique(data_ord$timeDSS)


unique(data_ord$mix_time)

unique(data_ord$cage)
data_ord$cage = factor(data_ord$cage,levels = c("S","U","OD","YD"),labels = c("OY_S_FMT","OY_U_FMT_DSS","OD","YD"))
data_ord=data_ord[data_ord$mix_time!=25,]

data_ord$donor2acceptor
length(unique(data_ord$timeDSS))
p<- pcoa_plot_size(data_ord,pcoa,metric=metric,title_name = "YY  PCoA-FMT-with or not DSS",size="mix_time" ,
                   # var1 = "donor2acceptor",var1_sub = c("OO","YO"),
                   color = "cage",shape = "donor2acceptor",shape_sub = c(21,22),
                   text="mix_time")+ scale_fill_manual(values = color_db[c("YO","YO_DSS","OD","YD"),"value"])+
  scale_color_manual(values = color_db[c("YO","YO_DSS","OD","YD"),"value"])+
  scale_size_continuous(range=c(1,12),breaks = c(seq(1,6)),labels = c(0,1,4,7,11,14))


p_data_mean=p$data 
p_data_mean[,c("pc1","pc2")]=p_data_mean[,c("pc1_mean","pc2_mean")]

new_data=rbind(p$data,p_data_mean)

data_donor = subset(p$data,cage %in% c("OD","YD"))

p
p1=p+ 
  geom_point(data=data_donor,aes(shape=donor2acceptor),
             stroke=3,size=10)+ scale_shape_manual(values = c(4,3))
# geom_point(data=p$data,aes(pc1,pc2,color=p$data$cage),size=3,alpha=0.5) +
# geom_line(data=new_data,aes(pc1,pc2,group=new_data$Row.names,color=cage),
#           size=0.5,alpha=0.5)



p1
# pdf(file=paste(file_name,paste("other YY FMT-DSS pcoa",metric,".pdf"),sep = "/"), height = 8, width = 8) #for pdf
# print(p1)
# dev.off()
