library(dplyr)
library(phyloseq)
library(ggplot2)
library(reshape2)
library(stringr)
library(patchwork)
sessionInfo()
#######
#### remove the data and keep the function and test time comsuming
options(warn=-1)
rm(list = setdiff(ls(), lsf.str()))


###
# otu：row for ID, col for sample ID
# tax(taxonomy)：row for sample ID, col for taxonomy
# metadata：row for sample ID, col for features(age, group, treatment)
# df：dataframe, row for row number, col for features, might need 'mix' for some variables
# counts：row for gene names, transcript ID or pathway, col for sample ID
# mkegg ： row for gene names，col for level1, level2/ enrich , row for pathway id, col for Description,geneID
# 
# list: many types of data (phyloseq data labels as phylo_ls( otu: row for sample ID, tax is matrix type)
# Data：many tpyes of data, phyloseq or df data

####

#####data prepocessing------------------
otu_filter=function(otu,col_lim,row_lim=0){
  ### col_lim : filter richness of asv in each sample， row_lim : filter richness of asv between samples
  ## col filter
  otu1 <- sweep(otu,2,colSums(otu),'/')
  otu1<- ifelse(otu1>col_lim,1,0)
  otu2 <- otu*otu1
  h <- rowSums(otu2)%>%as.data.frame()
  colnames(h) <- "Row_sum"
  #row filter
  h$rate <- sweep(h,1,colSums(h),"/")
  select_name <- rownames(h[h$rate>row_lim,])
  otu3 <- otu2[select_name,] 
  
  otu4<- round(sweep(otu3,2,colSums(otu3)/min(colSums(otu3)),"/"))
  
  return(otu4)
  
  
}

#### taxonomy, from qiime2 output dataset
tax_make = function(file,target_name="taxonomy.csv"){
  tax <- read.table(file,sep = "\t",header = TRUE)
  tax <- tax[,1:2]
  tax$Kindom <- sapply(str_split(tax$Taxon,";"),"[",1)
  tax$Phylum <- sapply(str_split(tax$Taxon,";"),"[",2)
  tax$Class <- sapply(str_split(tax$Taxon,";"),"[",3)
  tax$Order <- sapply(str_split(tax$Taxon,";"),"[",4)
  tax$Family <- sapply(str_split(tax$Taxon,";"),"[",5)
  tax$Genus <- sapply(str_split(tax$Taxon,";"),"[",6)
  tax$Species <- sapply(str_split(tax$Taxon,";"),"[",7)
  tax <- tax[,-2]
  colnames(tax)[1]<- "Feature_ID"
  write.csv(tax,target_name,row.names = FALSE)
  
}

tax_detail=function(tax,lefse=FALSE,prefix=TRUE){
  #### remember to add na.strings = "",fill = TRUE ,no factor for the string
  ### add prefix and fill the na with upper taxonomy
  tax<- tax[-which(is.na(tax$Phylum)),]
  # tax <- tax[!is.na(tax$Phylum),]
  if(is.factor(tax[,1])){
    stop("it's factor")
  }else if(isFALSE(lefse)){
    for(i in 1:nrow(tax)){
      for (j in 2:ncol(tax)) { # match the first not with un
        if (is.na(str_match(tax[i,j],"^[uU]n"))==FALSE | tax[i,j]=="NA" | is.na(tax[i,j])){
          tax[i,j]=tax[i,j-1]
        }else if (is.na(tax[i,j])){
          tax[i,j]=tax[i,j-1]
        }else{
          if(isTRUE(prefix)){
            tax[i,j]=paste(substr(colnames(tax[j]),1,1),tax[i,j],sep = "_")
          }else{
            tax[i,j]=tax[i,j]
          }
          
        }
        
      }
    }
    
  }else { # for lefse , don't replace the na or unculture, but keep them
    for(i in 1:nrow(tax)){
      for (j in 2:ncol(tax)) {
        if ( tax[i,j]=="NA" | is.na(tax[i,j])){
          next
        }else if(is.na(str_match(tax[i,j],"[uU]n"))==FALSE ){
          tax[i,j] = "NA"
        }
        else{
          tax[i,j]=paste(substr(colnames(tax[j]),1,1),tax[i,j],sep = "_")
        }
      }
    }
  }
  
  return(tax)
}


####### mean or other simple  functions
SEM=function(x){ ##x refers to array, or a col of a dataframe
  x1=sd(x)/(length(x)^0.5)
  return(x1)
}

### data type change
factor2character=function(data){
  for (i in 1:ncol(data)) {
    data[,i] <- as.matrix(data[,i]) 
  }
  return(data)
}

#### only return the mean value  and sem value
get_mean=function(df,var1,var2,var3={},var4={},value="value"){
  ## mean and sem based on group and time for dataframe
  df$mix <- paste(df[,var1],df[,var2],df[,var3],df[,var4],sep = "_")
  
  df1 <- data.frame(mean=tapply(df[,value],df$mix,mean),
                    se=tapply(df[,value],df$mix,SEM))
  
  df1[,var1]<- sapply(str_split(rownames(df1),"_"),"[",1)
  df1[,var2] <- sapply(str_split(rownames(df1),"_"),"[",2)
  if(!is.null(var3)){
    df1[,var3] <- sapply(str_split(rownames(df1),"_"),"[",3)
    if(!is.null(var4)){
      df1[,var4] <- sapply(str_split(rownames(df1),"_"),"[",4)
    }else{
      df1$mix <- paste(df1[,var1],df1[,var2],df1[,var3],sep = "_")
    }
    
  }else{
    df1$mix <- paste(df1[,var1],df1[,var2],sep = "_")
  }
  
  
  return(df1)
}


#### add mean row, for FMT mean data based on age
add_group_mean=function(counts,meta,time,group,mean_second={},mean_time=".*_0",remove_baseline1=FALSE,
                        remove_baseline2=FALSE){
  ### chang row and col for mean adding
  counts=as.data.frame(t(counts))
  if(c("mix") %in% colnames(meta)){}else{
    meta$mix=paste(meta[,group],meta[,time],sep = "_")
  }
  meta1=meta
  ##### mean each mix first time
  for (i in unique(na.omit(str_match(meta1$mix,mean_time)))) {
    counts_mean <- apply( counts[rownames(meta1[meta1$mix==i,]),],2,mean) 
    counts=as.data.frame(counts)
    
    counts[paste(i,"mean",sep = "~"),]<- counts_mean
    meta1[paste(i,"mean",sep = "~"),]<- paste(i,"mean",sep = "~")
    
    meta1[paste(i,"mean",sep = "~"),time] <- str_split(i,"_")[[1]][2]
    
    meta1[paste(i,"mean",sep = "~"),group] <- str_split(i,"_")[[1]][1]
    if(isTRUE(remove_baseline1)){ ### whether to remove the raw data used for mean
      counts <- counts[-which(rownames(counts) %in% rownames(meta1[meta1$mix==i,])),]
      meta1 <- meta1[-which(rownames(meta1) %in% rownames(meta1[meta1$mix==i,])),]
    }
  }
  ####second mean, group some of the first mean to second mean, only used fo donor data
  if(!is.null(mean_second)){
    for (j in mean_second) {
      name= paste(j,"_0~mean",sep = "")
      meta2=meta1[rownames(subset(meta1,rownames(meta1) %in% name),j),]
      counts_mean2 <- t(data.frame(apply( counts[rownames(meta2),],2,mean)))
      if(str_split(name[1],"")[[1]][2]=="Y"){
        new_name <- "YD"
        
      }else if (str_split(name[1],"")[[1]][2]=="O"){
        new_name <- "OD"
      }
      rownames(counts_mean2)<- new_name
      counts <- rbind(counts,counts_mean2)
      meta1[new_name,group]<-  new_name
      meta1[new_name,time] <- 0
      
      
    }
    if(isTRUE(remove_baseline2)){ ### whether to remove the first mean data used for mean
      counts <- counts[-which(!is.na(str_match(rownames(counts),"mean"))),]
      meta1 <- meta1[-which(!is.na(str_match(rownames(meta1),"mean"))),]
    }
    
  }
  counts=as.data.frame(t(counts))
  return(list(counts=counts,meta=meta1) )
}

#### add mean based on group1, mean var1,var2,var3 individually
dd_rowmean = function(data1,var1,var2={},var3={},
                      group1,group2={}){
  ###move the zero factors 
  if(!is.null(group2)){
    data1$mix <- paste(data1[,group1],data1[,group2],sep = "_")
    h1 <- rowsum(data1[,c(var1,var2,var3)],group = data1$mix)/summary(as.factor(data1$mix))
    colnames(h1) <- paste(c(var1,var2,var3),"mean",sep = "_")
    res <- merge(data1,h1,by.x="mix",by.y="row.names")
  }else{
    divided <- summary(as.factor(data1[,group1]),maxsum = Inf) ### set to Inf
    divided <- divided[divided != 0]
    
    h1 <- rowsum(data1[,c(var1,var2,var3)],group = data1[,group1]) / divided 
    colnames(h1) <- paste(c(var1,var2,var3),"mean",sep = "_")
    res <- merge(data1,h1,by.x=group1,by.y="row.names")
  }
  
  return(res)
}


dd_rowmean2 = function(data1,var1,var2={},var3={},
                       group1){
  
  divided <- summary(as.factor(data1[,group1]),maxsum = Inf) ### set to Inf
  divided <- divided[divided != 0]
  
  h1 <- rowsum(data1[,c(var1,var2,var3)],group = data1[,group1]) / divided 
  colnames(h1) <- paste(c(var1,var2,var3))
  data1[rownames(h1),]<- 0
  data1[rownames(h1),colnames(h1)] <- h1
  
  return(data1)
}


#### make a dataframe from two var to make a matric like , to mean from row for dd_rowmean
## only use for mad calculation, var1 for row, var_group for col
demelt = function(df,value="value",var1,var_group={}){ ## var1 for col, var_group for row
  if(any(df[,var1]=="other")){  ###remove the other for it's so large
    df <- df[df[,var1]!="other",]
  }
  ### select the biggest cols length
  dmelt_tax <- as.data.frame(matrix(nrow = length((unique(df[,var1]))),
                                  ncol =max(summary(df$variable))))
  
  
  rownames(dmelt_tax) <- unique(df[,var1])
  
  df <- df[order(df[,var1]),]## ensure the variable in the same cols
  
  for (i in unique(df[,var1])) {
    d1 <-df[df[,var1]==i,value]
    dmelt_tax[i,1:length(d1)] <- df[df[,var1]==i,value]
  }
  if(!is.null(var_group)){
    colnames(dmelt_tax) <- unique(df[,var_group])
  }
  return(dmelt_tax)
} 


##### dodge in the ploting
x_dodge=function(df,group,scale=1){
  ## used for figure significant label
  h1=as.numeric(as.factor(df[,group]))
  df$Dodge=(h1-mean(unique(h1)))*scale
  return(df)
}

### selecting function

sub_select = function(df,var1={},var1_sub={},
                      var2={},var2_sub={},var3={},var3_sub={}){
  if (!is.null(var1)){
    df <- subset(df,df[,var1] %in% var1_sub)
    if(!is.null(var2)){
      df <- subset(df,df[,var2] %in% var2_sub)
      if(!is.null(var3)){
        df <-subset(df,df[,var3] %in% var3_sub)
      } 
    }
  }
  return(df)
}


phylo1_get = function(phylo_ls,meta,var1={},var1_sub={},
                      var2={},var2_sub={},var3={},var3_sub={}){
  if (!is.null(var1)){
    meta <- subset(meta,meta[,var1] %in% var1_sub)
    if(!is.null(var2)){
      meta <- subset(meta,meta[,var2] %in% var2_sub)
      if(!is.null(var3)){
        meta <-subset(meta,meta[,var3] %in% var3_sub)
      } 
    }
  }
  sample_name <- rownames(meta)
  
  phylo1 <- prune_samples(sample_name, phylo_ls)
  return(phylo1)
}



## statistics--------------
### t test or wilcox.test,remove NA, test groups could not be exactly the same , or all zero
####  p adjust  for RNA-seq, metagenome
sig_table_get = function(df,time,group,group_base="wt",group_change="str",
                         method="wilcox",sig_label=TRUE,sig_sign="*"){
  ### select base and change for comparison
  if(is.null(time)){
    df$df_time=1
    time="df_time"
  }
  for(i in unique(df[,time])) {
    sig1<-df[df[,time]==i,]
    base_value <-sig1[sig1[,group]==group_base,]
    change_value <- sig1[sig1[,group]==group_change,]
    
    if(sum(base_value$value)+sum(change_value$value)==0){
      next  ## if all the value  is zero, no need for significant test 
    }else{
      if(method=="t_test"){
        res <- t.test(x=  base_value$value,y=change_value$value)
      }else if(method=="wilcox"){
        res <- wilcox.test(x=  base_value$value,y=change_value$value)
      }
      df[df[,time]==i &df[,group]==group_change,"comparison"]=paste(group_base,group_change,sep = "~")
      df[df[,time]==i &df[,group]==group_change ,"pval"]=res$p.value
      df[df[,time]==i &df[,group]==group_change,"se_value"]=SEM(change_value$value)
      
      if(isTRUE(sig_label)){
        sig_sign <- sig_sign
        
        if(res$p.value < 0.0001){
          df[df[,group]==group_change & df[,time]==i,"sig_label"] =strrep(sig_sign,4)
          
        }else if(res$p.value < 0.001){
          df[df[,group]==group_change & df[,time]==i,"sig_label"] =strrep(sig_sign,3)
          
        }else if (res$p.value < 0.01){
          df[df[,group]==group_change & df[,time]==i,"sig_label"] =strrep(sig_sign,2)
          
        }else if (res$p.value < 0.05){
          df[df[,group]==group_change & df[,time]==i,"sig_label"] =strrep(sig_sign,1)
        }else {
          df[df[,group]==group_change & df[,time]==i,"sig_label"] <- ""
        }
      }
    }
    
    
    
  }
  return(df)
}

padjust_get=function(df,time,sig_sign="*"){
  for (i in 1:length(unique(df[,time]))) {
    df1=df[df[,time]==unique(df[,time])[i],]
    df1$p_adjust=p.adjust(df1$pval,method = "fdr")
    
    if(i==1){
      df_all=df1
    }else{
      df_all=rbind(df_all,df1)
    }
  }
  
  
  
  df_all$adjust_sig_label=""
  df_all=df_all[!is.na(df_all$pval),]
  #######
  df_all[df_all$p_adjust < 0.05,"adjust_sig_label"] =strrep(sig_sign,1)
  df_all[df_all$p_adjust < 0.01,"adjust_sig_label"] =strrep(sig_sign,2)
  df_all[df_all$p_adjust < 0.001,"adjust_sig_label"] =strrep(sig_sign,3)
  df_all[df_all$p_adjust < 0.0001,"adjust_sig_label"] =strrep(sig_sign,4)
  
  df=full_join(df,df_all)
  
}




################ ploting function --------------

errorbar_line_plot=function(df, var1,var1_sub,group,time,value="value",line={},
                            signifi=TRUE,padjust=TRUE,sig_com=list(),dis_index=1){ #sig_com based on group var
  ## selecting data
  data_sub <- sub_select(df,var1,var1_sub )
  ## based on two cols for mean and sd calculation
  data_sub$mix <- paste(data_sub[,group],data_sub[,time],sep = "_")
  
  
  
  
  stats=get_mean(df=data_sub,var1=group,var2=time,value=value)
  stats[,time] <- as.numeric(stats[,time])
  ### significant label
  
  if(isTRUE(signifi)){
    pdata <- data_sub
    
    ## for significant label position
    y_min <- min(pdata$value)
    y_max <- max(pdata$value)
    y_low <- y_min - (y_max-y_min)/10
    ## space for different comparison 
    
    for (i in 1:length(sig_com)) {
      
      tax2_sig <- sig_table_get(df = pdata,group = group,time = time,
                                group_base=sig_com[[i]][1],group_change=sig_com[[i]][2],
                                method="wilcox",sig_sign = "*")
      
      tax2_sig <- subset(tax2_sig,tax2_sig[,group] %in% sig_com[[i]])
      
      tax2_sig$position <-  y_low+(y_max-y_min)/50*dis_index
      
      
      dis_index=dis_index+5
      if(i==1){
        tax_sig_all = unique(tax2_sig[tax2_sig[,group]==sig_com[[i]][2],c("pval","sig_label","mix","position",time)])
      }else{
        tax_sig_all =  full_join(tax_sig_all,unique(tax2_sig[tax2_sig[,group]==sig_com[[i]][2],c("pval","sig_label","mix","position",time)]))
      }
    }
    
    
    if(isTRUE(padjust)){
      tax_sig_all=padjust_get(df=tax_sig_all,time=time)
      tax_sig_all=tax_sig_all[,c("pval","sig_label","mix","position","p_adjust","adjust_sig_label")]
    }
    stats <-  full_join(unique(tax_sig_all[,c("pval","sig_label","mix","position","p_adjust","adjust_sig_label")]),stats)
  }else{
    stats$sig_label=0
    stats$position=0
  }
  
  
  
  
  
  stats[,group] <- factor(stats[,group],level=levels(df[,group]))
  
  if(!is.null(line)){  ## line for linetype
    stats <- left_join(stats,data_sub[,c("mix",line)])
    
    stats=x_dodge(stats,group = group,scale=0.5)
    
    
    
    ggplot(stats,aes(stats[,time],mean,color=stats[,group],group=stats[,group])) +
      
      
      geom_line(aes(y=mean,linetype=as.factor(stats[,line])),size=1.5)+
      
      geom_point(aes(x=stats[,time]-stats$Dodge),
                 position=position_dodge2(width=0.2),
                 size=3,shape=21,fill="white",stroke=2)+
      
      geom_errorbar(aes(x=stats[,time]-stats$Dodge,ymin=mean-se,ymax=mean+se),
                    # position = position_dodge2(width = 0.2),
                    width=1,size=1.5)+
      
      # scale_color_manual(values =color_db[c("Y0","Y2"),"value"] )+
      geom_text(aes(y=position ,label=sig_label,color=stats[,group]),size=10) +  ## use raw p value 
      
      labs(x="Days",y="log10 16S rRNA gene copies\n (per gram of wet weight)",
           color=group,linetype=line)+
      
      theme( panel.grid.minor = element_blank(),
             panel.background = element_blank(),panel.border = element_blank(),
             axis.line = element_line(colour = "black",linetype="solid",size = 1),
             panel.grid.major=element_line(colour=NA),axis.text=element_text(size = 15, face= "bold",family = ""),
             legend.text = element_text(size = 15),legend.title = element_text(size=15),axis.title=element_text(size=15))
    
  }else{
    stats[,group] <- factor( stats[,group],levels=levels(data_sub[,group]))
    
    stats=x_dodge(stats,group = group,scale=0.5)
    
    ggplot(stats,aes(stats[,time],mean,color=stats[,group],group=stats[,group])) +
      geom_line(aes(y=mean),size=1.5)+
      
      
      geom_errorbar(aes(x=stats[,time]-stats$Dodge,ymin=mean-se,ymax=mean+se),
                    # position = position_dodge2(width = 0.2),
                    width=1,size=1.5)+
      # scale_color_manual(values =color_db[c("Y0","Y2"),"value"] )+
      geom_text(aes(y=position ,label=sig_label,color=stats[,group]),size=10) +
      
      labs(x="Days",y="log10 16S rRNA gene copies\n (per gram of wet weight)",
           color=group,linetype=line)+
      
      theme( panel.grid.minor = element_blank(),
             panel.background = element_blank(),panel.border = element_blank(),
             axis.line = element_line(colour = "black",linetype="solid",size = 1),
             panel.grid.major=element_line(colour=NA),axis.text=element_text(size = 15, face= "bold",family = ""),
             legend.text = element_text(size = 15),legend.title = element_text(size=15),axis.title=element_text(size=15))
    
    
  }
}


### boxplot, only box_plot,with significant,melt type 
## com should be list like com={c(0,1),c(0,4)}, the ylim should change if no signal
base_boxplot= function(df,x,y,com={},fill={},facet={}){
  
  ## remove some of the levels with no data 
  df=arrange(df,by_group=df[,x])
  df[,x]=factor(df[,x],levels = unique(df[,x]))
  ### significant comparison group
  if(is.null(com)){
    com=combn(unique(as.character(df[,x])),2,simplify = FALSE)
  }
  
  ## significant calculation
  for (i in 1:length(com)) {
    sig_df=sig_table_get(df,time={},group = x, group_base = com[[i]][1],group_change = com[[i]][2]) %>% na.omit()
    if(i==1){
      sig_all=sig_df
    }else{
      sig_all=rbind(sig_all,sig_df)
    }
  }
  sig_all1=sig_all[,c(x,"comparison","pval","sig_label")]%>% unique()
  sig_all1$y=seq(from=max(df$value),to=1.5*max(df$value),by=0.1*max(df$value))[1:nrow(sig_all1)]
  # sig_all1[,x]=factor(sig_all1[,x])
  # sig_all1$x=factor(sig_all1[,x],)
  sig_all1$xend=sapply(str_split(sig_all1$comparison,"~"),"[",1) %>% factor(.,levels = levels(sig_all1[,x]))
  sig_all1$sig_label_y= sig_all1$y+0.01*max(df$value)
  sig_all1$sig_label_x=(as.numeric(sig_all1$xend)+as.numeric(sig_all1[,x]))/2
  ## p adjust
  sig_all1=padjust_get(df=sig_all1,time = "xend") %>% filter(adjust_sig_label != "")  
  
  
  
  ### make fill the same as x if no input
  if(is.null(fill)){
    fill="v1"
    df[,fill]<-as.factor(df[,x])
  }
  
  # return(sig_all1)
  
  if(is.null(facet)){
    
    
    
    y_max= 1.5*max(df[,y])
    ggplot(df,aes(x=df[,x],y=df[,y],color=df[,fill]))+
      geom_boxplot(aes(color=df[,fill]),position = position_dodge(0.8),size=1.5)+
      # guides(fill="none")+
      # facet_grid(.~df[,facet],scale="free")+
      labs(color=x,x=x,y=y)+
      ylim(NA,y_max*1.5)+
      geom_segment(data=sig_all1,aes(x=sig_all1[,x],y=y,xend=xend,yend=y),color="black",size=1.5)+ ## signifi
      geom_text(data=sig_all1,aes(x=sig_label_x,y=sig_label_y,label=adjust_sig_label),
                color="black",size=10)+
      # ggpubr::stat_compare_means(size=10, label='p.signif',bracket.size = 1.5,vjust=-.3,
      #                            # ggpubr::stat_compare_means(size=10, label='p.format',bracket.size = 1.5,
      #                            # label.y=1.0,tip.length = 0.1,
      #                            # comparisons = combn(unique(df[,x]),2,simplify = FALSE),
      #                            # comparisons = list(c("0","1"),c("0","4"),c("0","42")),
      #                            comparisons=com,
      #                            symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
      #                            method = "wilcox.test",hide.ns=TRUE)+
      geom_point(aes(color=df[,fill]),size=5,position = position_dodge(0.8))+
      theme( panel.grid.minor = element_blank(), panel.grid.major=element_line(colour=NA),
             panel.background = element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"),
             axis.text.x=element_text(angle=45,hjust=1,size = 15),
             legend.text = element_text(size = 15),legend.title = element_text(size=20),# for legend
             axis.text.y=element_text(size=15),axis.title = element_text(size=15), # for axis and its title
             strip.text=element_text(size=15))
    
    
    
  }else{
    ggplot(df,aes(x=df[,x],y=df[,y],color=df[,fill]))+
      geom_boxplot(aes(color=df[,fill]),position = position_dodge(0.8),size=1.5)+
      guides(fill="none")+
      facet_grid(.~df[,facet],scale="free")+
      ggpubr::stat_compare_means(size=10, label='p.signif',bracket.size = 1.5,vjust=-.3,
                                 # ggpubr::stat_compare_means(size=5, label='p.format',bracket.size = 1.5,
                                 tip.length = 0.1,
                                 symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                                 comparisons = com,
                                 method = "wilcox.test")+
      labs(color=fill,x=x,y=y)+
      geom_point(aes(color=df[,fill]),size=5,position = position_dodge(0.8))+
      theme( panel.grid.minor = element_blank(), panel.grid.major=element_line(colour=NA),
             panel.background = element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"),
             axis.text.x=element_text(angle=45,hjust=1,size = 15),
             legend.text = element_text(size = 15),legend.title = element_text(size=20),# for legend
             axis.text.y=element_text(size=15),axis.title = element_text(size=15), # for axis and its title
             strip.text=element_text(size=15))
    
  }
  
  
  
}


###significat boxplot
sig_box_plot=function(df,var1,var1_sub={},var2={},var2_sub={},
                      title=var1_sub,com={},time="variable",Value="value",facet={}){
  df1 <- sub_select (df,var1 ,var1_sub,
                     var2,var2_sub)
  df1 <- na.omit(df1)
  base_boxplot(df = df1,x=time,y=Value,com=com,facet=facet)+ 
    labs(title=title,y=paste(metric,"dissimilarity",sep = " ")) #+
  # ylim(NA,1.2)
}


############     alpha ploting  ---------------
### for one side significant label at the bottom
alpha_plot2 = function(phylo_ls,meta,xlab,color_col,measure="Observed",group={},var1={},var1_sub={},
                       line=group,signifi=TRUE,padjust=TRUE,sig_com={},scale=1,
                       var2={},var2_sub={},time={},average=FALSE,title_name,x_name=xlab,
                       var3={},var3_sub={},label_size=10){
  
  ### select data
  phylo1 <- phylo1_get(phylo_ls,meta,var1,var1_sub ,var2,var2_sub ,var3,var3_sub)
  
  
  ### get the data from richness function
  p<-plot_richness(phylo1, x=xlab,measures = measure)
  
  all <- p$data 
  rownames(all) <- paste("sample",p$data$id,sep = "")
  meta1 = meta[rownames(all),]
  
  
  
  if(!is.null(time)){
    p$data[,time] <- as.numeric(as.matrix(p$data[,time]))
  }
  if(isTRUE(average)){
    p1 <- p
    ##get mean value from value based on mix 
    p1$data= get_mean(df=all,var1=xlab,var2=color_col,var3="variable",var4=line)
    p1$data[,time]  <- as.numeric(as.matrix(p1$data[,time]))
    
    
    
    ### significant label position
    y_min <- min(p1$data$mean)
    y_max <- max(p1$data$mean)
    y_low <- y_min - 2*(y_max-y_min)/5
    
    
    if(isTRUE(signifi)){ ##  use raw data, no mean
      pdata <- p$data
      dis_index=1
      for (i  in 1:length(sig_com)) { ## method default wilcox test
        tax2_sig <- sig_table_get(df = pdata,time = time,group = group,
                                  group_base=sig_com[[i]][1],group_change=sig_com[[i]][2],
                                  method="wilcox",sig_sign = "*")
        tax2_sig <- subset(tax2_sig,tax2_sig[,group] %in% sig_com[[i]])
        
        tax2_sig$position <-  y_low+(y_max-y_min)/5/10*dis_index
        
        tax2_sig$mix<-  paste(tax2_sig[,xlab],tax2_sig[,color_col],tax2_sig$variable,tax2_sig[,line],sep = "_")
        
        dis_index=dis_index+5
        
        if(i==1){
          tax_sig_all = unique(tax2_sig[tax2_sig[,group]==sig_com[[i]][2],c("pval","sig_label","mix","position",time)])
        }else{
          tax_sig_all =  full_join(tax_sig_all,unique(tax2_sig[tax2_sig[,group]==sig_com[[i]][2],c("pval","sig_label","mix","position",time)]))
        }
        
      }
      
      if(isTRUE(padjust)){
        tax_sig_all=padjust_get(df=tax_sig_all,time=time)
        tax_sig_all=tax_sig_all[,c("pval","sig_label","mix","position","p_adjust","adjust_sig_label")]
      }
      
      
    }
    ### merge to p1 data for ploting 
    p1$data <- merge(p1$data,tax_sig_all,by.x="row.names",by.y = "mix",all.x=TRUE)
    
    
    
    p1$data=x_dodge(p1$data,group = group,scale=scale)
    
    ggplot(p1$data,aes(x=p1$data[,time],y=mean))+
      labs(x=paste(var1_sub,xlab,sep = " "),size=3,color=color_col,fill=color_col,linetype=line,
           title = title_name,x=x_name,y=measure)+
      geom_line(aes(group=p1$data[,group],color=p1$data[,color_col],
                    linetype=p1$data[,line]),size=1.5)+
      
      
      
      geom_errorbar(aes(x=p1$data[,time]-p1$data$Dodge,ymin=mean-se,ymax=mean+se,color=p1$data[,color_col]),
                    # position = position_dodge2(width = 0.2),
                    width=1,size=1.5)+
      geom_point(aes(x=p1$data[,time]-p1$data$Dodge,color=p1$data[,color_col]),
                 position=position_dodge2(width=0.2),
                 size=3,shape=21,fill="white",stroke=2)+
      
      
      geom_text(aes(y=position ,label=sig_label,color=p1$data[,group]),
                position = position_dodge2(width = 0.2),size=label_size) +   ##raw pvalue
      ylim(y_low,NA)+
      
      theme( panel.grid.minor = element_blank(),
             panel.background = element_blank(),panel.border = element_blank(),
             axis.line = element_line(colour = "black",linetype="solid",size = 1),
             panel.grid.major=element_line(colour=NA),axis.text=element_text(size = 15, face= "bold",family = ""),
             legend.text = element_text(size = 15),legend.title = element_text(size=15),axis.title=element_text(size=15))
    
    
  }
  
  
}

###### for upper and bottom side of the figure
alpha_plot3 = function(phylo_ls,meta,xlab,color_col,measure="Observed",group={},var1={},var1_sub={},
                       ###with 2 significant symbol
                       line=group,signifi=TRUE,padjust=TRUE,sig_com={},sig_com2={}, scale=1,
                       var2={},var2_sub={},time={},average=FALSE,title_name,x_name=xlab,
                       var3={},var3_sub={}){
  
  
  ### select data
  phylo1 <- phylo1_get(phylo_ls,meta,var1,var1_sub ,var2,var2_sub ,var3,var3_sub)
  
  
  ### get the data from richness function
  p<-plot_richness(phylo1, x=xlab,measures = measure)
  
  all <- p$data 
  rownames(all) <- paste("sample",p$data$id,sep = "")
  meta1 = meta[rownames(all),]
  
  
  if(!is.null(time)){
    p$data[,time] <- as.numeric(as.matrix(p$data[,time]))
  }
  if(isTRUE(average)){
    p1 <- p
    ##get mean value from value based on mix 
    p1$data= get_mean(df=all,var1=xlab,var2=color_col,var3="variable",var4=line)
    p1$data[,time]  <- as.numeric(as.matrix(p1$data[,time]))
    
    
    
    ### significant label
    y_min <- min(p1$data$mean)
    y_max <- max(p1$data$mean)
    y_low <- y_min - 2*(y_max-y_min)/5
    
    y_high <- y_max +  2*(y_max-y_min)/5
    
    
    if(isTRUE(signifi)){
      pdata <- p$data
      dis_index=1
      for (i in 1:length(sig_com)) {
        tax2_sig <- sig_table_get(df = pdata,group = group,time = time,
                                  group_base=sig_com[[i]][1],group_change=sig_com[[i]][2],
                                  method="wilcox",sig_sign = "#")
        tax2_sig <- subset(tax2_sig,tax2_sig[,group] %in% sig_com[[i]])
        
        tax2_sig$position <-  y_low+(y_max-y_min)/5/10*dis_index
        
        tax2_sig$mix<-  paste(tax2_sig[,xlab],tax2_sig[,color_col],tax2_sig$variable,tax2_sig[,line],sep = "_")
        
        dis_index=dis_index+5
        
        if(i==1){
          tax_sig_all = unique(tax2_sig[tax2_sig[,group]==sig_com[[i]][2],c("pval","sig_label","mix","position",time)])
        }else{
          tax_sig_all =  full_join(tax_sig_all,unique(tax2_sig[tax2_sig[,group]==sig_com[[i]][2],c("pval","sig_label","mix","position",time)]))
        }
        
      }
      
      if(isTRUE(padjust)){
        tax_sig_all=padjust_get(df=tax_sig_all,time=time)
        tax_sig_all=tax_sig_all[,c("pval","sig_label","mix","position","p_adjust","adjust_sig_label")]
      }
      
      if(!is.null(sig_com2)){ ## if two more sig_com2 placing in a different place from the sig_com, at the top
        dis_index=1
        for (i in 1:length(sig_com2)) {
          tax2_sig <- sig_table_get(df = pdata,group = group,time = time,
                                    group_base=sig_com2[[i]][1],group_change=sig_com2[[i]][2],
                                    method="wilcox")
          tax2_sig <- subset(tax2_sig,tax2_sig[,group] %in% sig_com2[[i]])
          
          tax2_sig$position2 <-  y_high-(y_max-y_min)/5/10*dis_index
          
          tax2_sig$mix<-  paste(tax2_sig[,xlab],tax2_sig[,color_col],tax2_sig$variable,tax2_sig[,line],sep = "_")
          
          dis_index=dis_index+5
          
          if(i==1){
            
            tax_sig_all2 = unique(tax2_sig[tax2_sig[,group]==sig_com2[[i]][2],c("pval","sig_label","mix","position2",time)])
          }else{
            tax_sig_all2 =  full_join(tax_sig_all2,unique(tax2_sig[tax2_sig[,group]==sig_com2[[i]][2],c("pval","sig_label","mix","position2",time)]))
          }
          
        }
        if(isTRUE(padjust)){
          tax_sig_all2=padjust_get(df=tax_sig_all2,time=time)
          tax_sig_all2=tax_sig_all2[,c("pval","sig_label","mix","position2","p_adjust","adjust_sig_label")]
          colnames(tax_sig_all2)<- c("pval2","sig_label2","mix","position2","p_adjust2","adjust_sig_label2")
        }
        
        tax_sig_all <- full_join(tax_sig_all,tax_sig_all2)
      }
      
    }
    
    
    p1$data <- merge(p1$data,tax_sig_all,by.x="row.names",by.y = "mix",all.x=TRUE)
    
    p1$data=x_dodge(p1$data,group = group,scale=scale)
    
    
    #### 
    ggplot(p1$data,aes(x=p1$data[,time],y=mean))+
      labs(x=paste(var1_sub,xlab,sep = " "),size=3,color=color_col,fill=color_col,linetype=line,
           title = title_name,x=x_name,y=measure)+
      geom_line(aes(group=p1$data[,group],color=p1$data[,color_col],
                    linetype=p1$data[,line]),size=1.5)+
      
      geom_errorbar(aes(x=p1$data[,time]-p1$data$Dodge,ymin=mean-se,ymax=mean+se,color=p1$data[,color_col]),
                    # position = position_dodge2(width = 0.2),
                    width=0.5,size=1.5)+
      geom_point(aes(x=p1$data[,time]-p1$data$Dodge,color=p1$data[,color_col]),
                 position=position_dodge2(width=0.2),
                 size=3,shape=21,fill="white",stroke=2)+
      labs(title = title_name,x=x_name,y=measure)+
      
      geom_text(aes(y=position ,label=sig_label,color=p1$data[,group]),size=5) + ## raw pvalue 
      geom_text(aes(y=position2 ,label=sig_label2,color=p1$data[,group]),size=10) +
      ylim(y_low,y_high)+
      
      theme( panel.grid.minor = element_blank(),
             panel.background = element_blank(),panel.border = element_blank(),
             axis.line = element_line(colour = "black",linetype="solid",size = 1),
             panel.grid.major=element_line(colour=NA),axis.text=element_text(size = 15, face= "bold",family = ""),
             legend.text = element_text(size = 15),legend.title = element_text(size=15),axis.title=element_text(size=15))
    
    
  }
  
  
}



################################  relative abundance-------------- 
## main function,the genus should be add with numerber to avoid lost

otu_tax1_get = function(phylo_all,tax_Name,meta,var={},var_sub={}){
  #########select a section from the whole
  if(!is.null(var)){
    sample_name <- rownames(meta[meta[,var] ==var_sub,])
    phylo1 <- prune_samples(sample_name, phylo_all)
  }else{
    phylo1 <- phylo_all
  }
  phylo1.tax = tax_glom(phylo1, tax_Name)
  
  otu_tax1=otu_table(phylo1.tax)
  tax <- tax_table(phylo1.tax) 
  # make a table , the row for sample, and the col for taxa ranks
  otu_tax1 <- as.data.frame(otu_tax1)
  colnames(otu_tax1) = tax[colnames(otu_tax1),tax_Name]
  if(tax_Name=="Genus"){
    colnames(otu_tax1) <- paste(colnames(otu_tax1),seq(1,ncol(otu_tax1)),sep = "_")
  }
  
  
  
  return(otu_tax1)
}

# for related abundance, to make tax_melt,   melting the time and var1 2,3
# meta1 should countain at least 2 vars besides rownames 
relative_abun_get=function(otu_tax1,meta1,time,time_type="num",var1={},var2={},
                           var3={},var4={}){
  otu_tax_rel = sweep(otu_tax1,1,rowSums(otu_tax1),'/') 
  otu_tax_rel$sampleId = rownames(otu_tax_rel)
  
  if (!is.null(time)){
    otu_tax_rel$V0 <- meta1[,time]
    colnames(otu_tax_rel)[ncol(otu_tax_rel)] <- time
    tax_melt=melt(otu_tax_rel,id.vars = c("sampleId",time))
  }
  if (!is.null(var1)){
    otu_tax_rel$V1 <- meta1[,var1]
    colnames(otu_tax_rel)[ncol(otu_tax_rel)] <- var1
    tax_melt=melt(otu_tax_rel,id.vars = c("sampleId",var1,time))
    if (!is.null(var2)){
      otu_tax_rel$V2 <- meta1[,var2]
      colnames(otu_tax_rel)[ncol(otu_tax_rel)]<- var2
      tax_melt=melt(otu_tax_rel,id.vars = c("sampleId",var1,var2,time))
      if (!is.null(var3)){
        otu_tax_rel$V3 <- meta1[,var3]
        colnames(otu_tax_rel)[ncol(otu_tax_rel)]<- var3
        tax_melt=melt(otu_tax_rel,id.vars = c("sampleId",var1,var2,var3,time))
        if (!is.null(var4)){
          otu_tax_rel$V4 <- meta1[,var4]
          colnames(otu_tax_rel)[ncol(otu_tax_rel)]<- var4
          tax_melt=melt(otu_tax_rel,id.vars = c("sampleId",var1,var2,var3,var4,time))
        }
      }
    }
  }
  meta1[,time] <- factor(meta1[,time])
  time_level = levels(meta1[,time])
  if(time_type=="num"){
    tax_melt[,time]= factor(as.numeric(as.matrix(tax_melt[,time])),levels = time_level)
  }else{
    tax_melt[,time]= factor(tax_melt[,time],levels = time_level)
  }
  
  return(tax_melt)
  
}

### select after relative_abun_get, 
relative_abun_second_get=function(tax_melt,tax0=tax0,Time="timeFMT",time_type="num",
                                  plot_number=15){
  top <- top_tax_get(tax_melt,Time=Time)
  
  
  ###only need  to factor , because the other will add up in the barplot
  top_all <- top[[1]]
  
  #### provide larger than 7 colors in the ploting,based on all top  
  tax_melt1 <- top[[2]]
  
  ### top based on each time 
  tax_melt1 <- top[[2]]
  
  top2 <- Top_taxonomy_time(tax_melt,time = Time,select_num = plot_number)
  top2$tax2 <- paste(sapply(str_split(top2$variable,"_"),"[",1) ,
                     sapply(str_split(top2$variable,"_"),"[",2),sep = "_" )
  top_all <- top[[1]]
  top_all2 <- subset(top_all,name %in% unique(top2$tax2))
  top_all2$order <- seq(1,nrow(top_all2))
  
  plot_number = nrow(top_all2)
  
  top_all <- full_join(top_all,top_all2)
  
  top_all[is.na(top_all$order),"order"] <-plot_number+1
  
  top_all <- top_all[order(top_all$order),]
  
  
  top_all_sub <- top_all[1:plot_number,]
  
  ### add the value from the lower level taxonomy to the upper level taxonomy
  h <- subset(tax0[,c("Family","Genus")],Family %in% top_all_sub[!is.na(str_match(top_all_sub$name,"F_")),"name"]) %>% unique()
  
  colnames(h) <- c("new_name","name")
  
  top_all_sub <- left_join(top_all_sub,h)
  
  for (i in 1:nrow(top_all_sub)) {
    if(is.na(top_all_sub$new_name[i])){
      top_all_sub$new_name[i] <- top_all_sub$name[i]
    }
  }
  
  tax_melt1$variable <- factor(tax_melt1$variable,levels = c(top_all$name[(plot_number+1):nrow(top_all)],rev(top_all_sub$name[1:plot_number])),
                               labels = c(rep("others",length(top_all$name[(plot_number+1):nrow(top_all)])),rev(top_all_sub$new_name[1:plot_number])))
  
  ##### relevels again 
  tax_melt1$variable <- factor(tax_melt1$variable,levels = c("others",rev(unique(top_all_sub$new_name))))
  ## if time is number then change it 
  if(Time=="num"){
    tax_melt1[,Time] <- as.numeric(as.matrix(tax_melt1[,Time]))
  }
  
  
  return(tax_melt1)  
  
}





# to get the top taxonomy ,based on which group is avalible, time is always used and individual is used
Top_taxonomy = function(tax_melt,top_number=10,top_percent={},time,time_sub={},individual={}){
  if(!is.null(time_sub)){
    tax_melt <- tax_melt[tax_melt[,time]==time_sub,]
  }
  if(!is.null(individual)){
    tax_melt$fullname=paste(tax_melt[,time],tax_melt[,individual],tax_melt$variable,sep = "~")
    #mean
    tax_melt2<-tapply(tax_melt[,'value'],tax_melt$fullname,mean)%>% as.data.frame()
    split_name <- t(as.data.frame(str_split(rownames(tax_melt2),"~")))
    tax_melt2$V2 <- split_name[,1]
    tax_melt2$V3 <- split_name[,2]
    tax_melt2$V4 <- split_name[,3]
    colnames(tax_melt2)<- c("value",time,individual,"variable")
    
  }else{
    tax_melt$fullname=paste(tax_melt[,time],tax_melt$variable,sep = "~")
    #mean based on the fullname to mean 
    tax_melt2<-tapply(tax_melt[,'value'],tax_melt$fullname,mean)%>% as.data.frame()
    split_name <- t(as.data.frame(str_split(rownames(tax_melt2),"~")))
    tax_melt2$V2 <- split_name[,1]
    tax_melt2$V3 <- split_name[,2]
    colnames(tax_melt2)<- c("value",time,"variable")
  }
  
  
  ## based on the variable to mean 
  hh <- tapply(tax_melt2$value,tax_melt2$variable,mean) #%>% as.data.frame() 
  
  hh <- sort(hh,decreasing = TRUE)
  if(!is.null(top_percent)){
    top <- rownames(hh[hh > top_percent]) %>% as.character(.)
  }else{
    top <- hh[1:top_number]  %>% as.data.frame(.)
  }
  colnames(top) <- "percent"
  top$variable <- rownames(top)
  return(top)
}

###use for getting each day's top taxonomy ,top 10 each day,no more than 50
Top_taxonomy_time=function(tax_melt,time,select_num=10){
  count=1
  for (i in unique(tax_melt[,time])) {
    
    if(count==1){
      top_all =  Top_taxonomy(tax_melt,top_number = 10,time = time,time_sub = i)
      colnames(top_all) <- "percent_sum"
      genus_select <- rownames(top_all)[1:select_num] %>% as.data.frame()
      colnames(genus_select) <- "variable"
      genus_select$time <- i
      genus_select$percent_sum <- top_all[1:select_num,"percent_sum"]
    }else{
      top_all =  Top_taxonomy(tax_melt,top_number = 50,time = time,time_sub = i)
      colnames(top_all) <- "percent_sum"
      h1 <- rownames(top_all)[1:select_num] %>% as.data.frame()
      colnames(h1) <- "variable"
      h1$time <- i
      h1$percent_sum<- top_all[1:select_num,"percent_sum"]
      genus_select <- rbind(genus_select,h1)
    }
    count=count+1
    
  }
  return(genus_select)
} 


top_tax_get= function(tax_melt,Time,mad){
  length(unique(tax_melt$variable))
  genus_select = Top_taxonomy(tax_melt,top_number = length(unique(tax_melt$variable)),time=Time)
  
  
  tax_melt$variable <- paste(sapply(str_split(tax_melt$variable,"_"),"[",1),
                             sapply(str_split(tax_melt$variable,"_"),"[",2),sep = "_")
  
  genus_select$variable <- paste(sapply(str_split(genus_select$variable,"_"),"[",1),
                                 sapply(str_split(genus_select$variable,"_"),"[",2),sep = "_")
  
  # length(unique(genus_select1$variable))
  top_all = tapply(genus_select$percent,genus_select$variable,mean)
  top_all = top_all[order(top_all,decreasing = TRUE)] %>% as.data.frame()
  top_all$name <- rownames(top_all)
  return(list(top_all,tax_melt))
}



# to get no more than 74 colors
distinctive_colors <- function(number = 74) {
  library(RColorBrewer)
  qual_col_pals <-  brewer.pal.info[brewer.pal.info$category == 'qual',]
  distinctive_colors <-  unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  distinctive_colors <- distinctive_colors[1:number]
  return(distinctive_colors)
}




# plot the relative abundance figure
### when treated with phylum, remeber to remove only kingdom at the otu_tax step 
###tax only plot seven color and grey for other 
###sort the level

relative_abun_plot=function(tax_melt1,var1,xvalue,facet,facet2={},tax="Genus",
                            plot_sort={},mean=FALSE,var1_sub={},var2={},var2_sub={},var3={},var3_sub={},
                            color_label=TRUE,color_g=74,title_name={},xvalue_order=TRUE,
                            genus_color=c("grey",distinctive_colors(color_g-8),
                                          "purple","red","brown","yellow","blue","green","orange")){
  if(!is.null(var1_sub)){
    data <- subset(tax_melt1,tax_melt1[,var1] %in% var1_sub)
    if(!is.null(var2_sub)){
      data <- subset(data,data[,var2] %in% var2_sub)
      if(!is.null(var3_sub)){
        data <- subset(data,data[,var3] %in% var3_sub)
      }
    }
  }else{
    data <- tax_melt1
  }
  ###### sort the levels and order the data in order to plot in encreasing order in phylum
  if(tax=="Phylum"){
    tax_sum <- aggregate(value~variable,data,sum)
    tax_sum <- tax_sum[order(tax_sum$value),]
    data$variable<-factor(data$variable,levels = tax_sum$variable)
    data<- data[order(data$variable),]
    
    
    
  }
  
  ### to order the xvalue 
  if(isTRUE(xvalue_order)){
    data[,xvalue] <-factor(data[,xvalue],levels = unique(data[,xvalue])[order(as.numeric(str_extract(unique(data[,xvalue]),"[0-9]+")))])
  }
  
  
  if(isTRUE(mean)){
    data$mix <- paste(data[,facet],data[,"variable"],sep = "_")
    h1 <- tapply(data$value,data$mix,mean)
    data2 <-h1 %>% as.data.frame()
    data2$v1 <- sapply(str_split(rownames(h1),"_"),"[",1)
    data2$v2 <- paste(sapply(str_split(rownames(h1),"_"),"[",2),
                      sapply(str_split(rownames(h1),"_"),"[",3),sep = "_")
    colnames(data2) <- c("value",facet,"variable")
    if(tax=="Genus"){
      ggplot(data2, aes(x = data2[,facet], y = value, fill = variable)) + 
        geom_bar(stat = "identity") + 
        scale_fill_manual(values = distinctive_colors(74))+
        labs(x=facet,y="relative_abundance",title = paste(ifelse(is.null(var2_sub),ifelse(is.null(var1_sub),var1,var1_sub),var2_sub)),fill="taxonomy") +
        theme(legend.position = "right",axis.text.x=element_text(angle=0)) # change the angle of the text 
      
      
    }else if(tax=="Phylum"){
      ggplot(data2, aes(x = data2[,xvalue], y = value, fill = variable)) + 
        geom_bar(stat = "identity") + 
        scale_fill_manual(values = as.matrix(unique(tax_melt1$name_col)))+ # for phylum
        scale_x_continuous(limits = c(2,3))+
        labs(x=xvalue,y="relative_abundance",title = title_name,fill="taxonomy") +
        #scale_fill_manual(values = as.character(unique(tax_melt_new$name))) + 
        theme(legend.position = "right", axis.text.x=element_text(angle=0)) # change the angle of the text 
      
    }
  }else{
    if(tax=="Genus"){
      if(isTRUE(color_label)){
        ggplot(data, aes(x = data[,xvalue], y = value, fill = variable)) + 
          facet_grid(cols = vars(data[,facet]),scale="free",space="free")+
          geom_bar(stat = "identity") + 
          # scale_fill_manual(values = c(distinctive_colors(color_g-7),"purple","orange","brown","yellow","blue","green","red"))+
          scale_fill_manual(values = genus_color)+
          #scale_fill_manual(values = c(gray.colors(color_g-7,),"purple","orange","brown","yellow","blue","green","red"))+
          # labs(x=xvalue,title = paste(ifelse(is.null(var2_sub),ifelse(is.null(var1_sub),var1,var1_sub),var2_sub)),
          # fill="taxonomy") +
          labs(x=xvalue,y="relative abundance",title = title_name, fill="taxonomy") +
          theme_bw() + theme(panel.grid = element_blank(), legend.text = element_text(size = 15),
                             axis.text = element_text(size=15),axis.title=element_text(size=20),
                             strip.text=element_text(size=20)) +
          #scale_fill_manual(values = as.character(unique(tax_melt_new$name))) +
          # for family
          #scale_fill_manual(values = c('red','blue', '#00FF00', 'burlywood1', 'purple', 'orange', 'hotpink', '#FF00FF', 'yellow', 'green', 'cyan', '#995599','gray')) + 
          theme(legend.position = "right",axis.text.x=element_text(angle=0),
                panel.border = element_blank()) # change the angle of the text 
        #ggplot(df_melt, aes(x = id, y = value, fill = variable)) + 
        #geom_bar(stat = "identity") + scale_fill_manual(values = distinctive_colors(55)) + theme(legend.position = "top",axis.text.x=element_text(angle=90)) 
        
      }else(
        ggplot(data, aes(x = data[,xvalue], y = value, fill = variable)) +
          facet_grid(cols = vars(data[,facet]),scale="free",space = "free")+
          geom_bar(stat = "identity") +
          scale_fill_manual(values = distinctive_colors(74))+
          labs(x=xvalue,y="relative_abundance",title = title_name,
               fill="taxonomy") +
          guides(fill= "none") +
          theme_bw() + theme(panel.grid = element_blank(), legend.text = element_text(size = 15),
                             axis.text = element_text(size=15),axis.title=element_text(size=20),
                             strip.text=element_text(size=20)) +
          # scale_fill_manual(values = as.character(unique(tax_melt_new$name))) +
          # for family
          #scale_fill_manual(values = c('red','blue', '#00FF00', 'burlywood1', 'purple', 'orange', 'hotpink', '#FF00FF', 'yellow', 'green', 'cyan', '#995599','gray')) + 
          theme(legend.position = "right",axis.text.x=element_text(angle=0),
                panel.border = element_blank()) # change the angle of the text 
        #ggplot(df_melt, aes(x = id, y = value, fill = variable)) + 
        #geom_bar(stat = "identity") + scale_fill_manual(values = distinctive_colors(55)) + theme(legend.position = "top",axis.text.x=element_text(angle=90)) 
        
      )
      
    }else if(tax=="Phylum"){
      if(!is.null(facet2)){
        data$facet2 <- paste(data[,facet],data[,facet2],sep = "_")
        data[,facet] <- data$facet2
      }
      # data$variable <- factor(data$variable,levels = unique(data$variable))
      data$name_col <- factor(data$name_col,levels= unique(as.character(data$name_col)))
      #   levels(data$name_col)
      unique(data$variable)
      ggplot(data, aes(x = data[,xvalue], y = value, fill = variable)) + 
        facet_grid(cols = vars(data[,facet]),scale="free",space = "free")+
        geom_bar(stat = "identity") + 
        #margin(5,0,5,0)+
        #geom_histogram(binwidth = 0.3)+
        #scale_y_reverse() +
        # scale_fill_manual(values =levels(tax_melt1$name_col))+ # for phylum
        scale_fill_manual(values =as.matrix(unique(data$name_col))) +
        labs(x=xvalue,y="relative_abundance",title =title_name,fill="taxonomy") +
        #scale_fill_manual(values = as.character(unique(tax_melt_new$name))) + 
        theme(legend.position = "right", axis.text.x=element_text(angle=0),
              plot.margin = margin(0, 1.5, 0, 1.5, "cm"),
              panel.border = element_blank(), ###delete the border
              panel.grid.major = element_blank()
        ) # change the angle of the text 
      
    }
  }
  
  
}



## line
#### for line plot 
mean_line_get= function(dis_melt0,time,group,value){
  dis_melt<- dis_melt0
  dis_melt[,group] <- as.factor(dis_melt[,group])
  dis_melt[,value] <- as.numeric(dis_melt[,value])
  dis_melt <- na.omit(dis_melt)
  dis_melt$mix <- paste(dis_melt[,group],dis_melt[,time],sep = "_")
  # dis_melt$variable <- as.numeric(as.matrix(dis_melt$variable))
  ###through rowsum to calculate the mean 
  dis_mean <- rowsum(dis_melt[,value],group=dis_melt$mix) /summary(as.factor(dis_melt$mix))
  colnames(dis_mean) <- "mean"
  dis_melt <- merge(dis_melt,dis_mean,by.x = "mix",by.y = "row.names")
  return(dis_melt)
  
}

line_plot2 =function(tax_melt,var1={},var1_sub={},var2={},var2_sub={},Value="value",
                     group,time,tax_num=7,mad=TRUE,title_name="",
                     mad_all=FALSE,color_level={},
                     scale_color=c( distinctive_colors(tax_num-7),'black', 'grey', 'purple', 'cyan', 'orange', 
                                    'red', 'yellow', 'green','brown','blue' )){
  
  ##### select some var1
  if(!is.null(var1_sub)){
    tax1 <- tax_melt[tax_melt[,var1]==var1_sub, ]
  }else{
    tax1 <- tax_melt
  }
  
  tax1=sub_select(tax_melt,var1,var1_sub,var2,var2_sub)
  
  tax2 <- mean_line_get(dis_melt0=tax1,group=group,time=time,
                        value=Value)
  ### get the top mad taxonomy ,h2 for row mad, row for taxonomy and col for time, select the large changes,
  if(isTRUE(mad_all)){
    count=0
    for (i in unique(tax_melt[,var1])) {
      h1 <- tax_melt[tax_melt[,var1]==i, ]
      
      h2 <- demelt(df=h1,value=Value,var1=group)
      ##### based on row taxonomy to select the large changes
      h_mad <- data.frame(variable=rownames(h2),mad_value=apply(h2, 1,mad,constant=1 ))
      h_mad <- h_mad[order(h_mad$mad_value,decreasing = TRUE),]
      h_mad$id=h_mad$variable
      if(count==0){
        top_mad <- as.character(h_mad$id[1:tax_num])
        count=count+1
      }else{
        
        top_mad <- append(top_mad,as.character(h_mad$id[1:tax_num]))
      }
    }
    top_mad <- unique(top_mad)
    tax2 <- subset(tax2,variable %in% top_mad)
    
  }else if (isTRUE(mad)){
    tax2_mad <- demelt(df=tax2,value=Value,var1="variable")
    tax2_mad2 <- data.frame(id=rownames(tax2_mad),mad_value=apply(tax2_mad, 1,mad,constant=1 ))
    tax2_mad2 <- tax2_mad2[order(tax2_mad2$mad_value,decreasing = TRUE),]
    tax2 <- subset(tax2,variable %in% tax2_mad2$id[1:tax_num])
  }
  
  ### set the same color as the related abundance
  if(!is.null(color_level)){
    color_level1 <- subset(color_level,color_level$taxonomy %in% unique(tax2$variable))
    scale_color <- as.character(color_level1$genus_color)
  }
  
  
  unique(tax2$variable)
  
  tax2[,time] <- as.numeric(as.matrix(tax2[,time]))
  
  
  
  
  ggplot(tax2,aes(x = tax2[,time], y =  tax2[,Value], group = tax2[,group])) +
    geom_line(aes(y=mean,color=tax2[,group]),size=2)+
    geom_point(aes(color=tax2[,group]),size=1,position = position_dodge2(width = 0.2))+
    labs(x = 'treatment_time', y = 'Relative Abundance', color = "taxonomy",
         title = paste(var1_sub,sep = ""))+
    scale_colour_manual(values = scale_color) +
    ylim(0,1)+
    theme( panel.grid.minor = element_blank(),
           panel.background = element_blank(),panel.border = element_blank(),
           axis.line = element_line(colour = "black",linetype="solid",size = 1),
           panel.grid.major=element_line(colour=NA),axis.text=element_text(size = 15, face= "bold",family = ""),
           legend.text = element_text(size = 10),legend.title = element_text(size=15),axis.title=element_text(size=15))
  
  
}


######################## PCoA-----------------

##### PCoA
### could be used with data.frame to calculate the bray crutis distance
### rowSums could not be zero
pcoa_product = function(Data, meta, metric, trefile){
  if (metric %in% c("unifrac", "wunifrac")){
    if (is.null(trefile)){
      stop("Phylogenic tree should be specified for Unifrac or Weighted-Unifrac distance")
    }else{
      library(phyloseq)
      
      sample_data = sample_data(meta)
      otu_table = otu_table(Data, taxa_are_rows = TRUE)
      
      if(metric=="unifrac"){
        Distance = rbiom::unifrac(as.matrix(t(otu_table)),weighted=FALSE,phy_tree(Data))  %>% as.matrix()
      }else{
        Distance = rbiom::unifrac(as.matrix(t(otu_table)),weighted=TRUE,phy_tree(trefile))  %>% as.matrix()
      }
      pcoa = cmdscale(Distance, k = 3, eig = T)
    }
    
  }else{
    if(class(Data)=="data.frame"){ ### for tax_melt , the dataframe will be turned to sample as row
      Distance <- vegan::vegdist(t(Data), method = metric)
    }else{
      Distance <- vegan::vegdist(otu_table(Data), method = metric)
    }
    pcoa = cmdscale(Distance, k = 3, eig = T)  # cal the distance and get the top 3 eigen value
  }
  
  
  mds1 = as.data.frame(pcoa$points)
  colnames(mds1) = c("pc1","pc2","pc3")
  data_ord=merge(mds1, meta, by = "row.names")
  row.names(data_ord) = data_ord$Row.names
  
  result = list(data_ord=data_ord,pcoa=pcoa,Distance=Distance)
  return(result)
}

##### for a mean base distance, all the time points distance to the mean base distance,and remove the base time value
##### filter after mean, if two group needed to mean ,set a mean_second in group column
#### only used for donor data,gropu_sub is OD or YD,used with mean_second
pcoa_distance_product2 = function(data, meta, metric, trefile,time="timeDSS",remove_baseline=TRUE,
                                  time_base=0,group="donor2acceptor",mean_second={},group_sub={}){
  if(is.data.frame(data)){
    otu_table=data
  }else{
    library(phyloseq)
    otu_table = otu_table(data, taxa_are_rows = TRUE)%>% as.data.frame()
    
  }
  ### add mean value
  meta1 <- meta[rownames(otu_table),]
  meta1$mix <- paste(meta1[,group],meta1[,time],sep = "_")
  #### first mean
  for (i in unique(na.omit(str_match(meta1$mix,".*_0")))) {
    
    otu_mean <- apply( otu_table[rownames(meta1[meta1$mix==i,]),],2,mean) 
    otu_table[paste(i,"mean",sep = "~"),]<- otu_mean
    meta1[paste(i,"mean",sep = "~"),]<- paste(i,"mean",sep = "~")
    meta1[paste(i,"mean",sep = "~"),time] <- str_split(i,"_")[[1]][2]
    meta1[paste(i,"mean",sep = "~"),group] <- str_split(i,"_")[[1]][1]
    if(isTRUE(remove_baseline)){
      otu_table <- otu_table[-which(rownames(otu_table) %in% rownames(meta1[meta1$mix==i,])),]
      meta1 <- meta1[-which(rownames(meta1) %in% rownames(meta1[meta1$mix==i,])),]
    }
    
  }
  ####second mean
  if(!is.null(mean_second)){
    for (j in mean_second) {
      name= paste(j,"_0~mean",sep = "")
      meta2=meta1[rownames(subset(meta1,rownames(meta1) %in% name),j),]
      otu_mean2 <- t(data.frame(apply( otu_table[rownames(meta2),],2,mean)))
      if(str_split(name[1],"")[[1]][2]=="Y"){
        new_name <- "YD"
        
      }else if (str_split(name[1],"")[[1]][2]=="O"){
        new_name <- "OD"
      }
      rownames(otu_mean2)<- new_name
      otu_table <- rbind(otu_table,otu_mean2)
      meta1[ new_name,]<-  new_name
      meta1[new_name,time] <- 0
      
    }
    if(isTRUE(remove_baseline)){
      otu_table <- otu_table[-which(!is.na(str_match(rownames(otu_table),"mean"))),]
      meta1 <- meta1[-which(!is.na(str_match(rownames(meta1),"mean"))),]
    }
    
    meta1 <- subset(meta1,meta1[,group] %in% group_sub)
    otu_table <- otu_table[rownames(meta1),]
  }
  
  
  ### distance calculation
  if(is.data.frame(data)){
    Distance <- vegan::vegdist(otu_table, method = metric) %>% as.matrix()
    
    
  }else{
    phylo1 = phyloseq::phyloseq(otu_table(otu_table, taxa_are_rows = F),sample_data(meta1), phy_tree(trefile))
    if (metric %in% c("unifrac", "wunifrac")){
      if (is.null(trefile)){
        stop("Phylogenic tree should be specified for Unifrac or Weighted-Unifrac distance")
      }else{
        if(metric=="unifrac"){
          Distance = rbiom::unifrac(as.matrix(t(otu_table)),weighted=FALSE,phy_tree(data))  %>% as.matrix()
        }else{
          Distance = rbiom::unifrac(as.matrix(t(otu_table)),weighted=TRUE,phy_tree(trefile))  %>% as.matrix()
        }
      }
      
    }else if (metric %in% c("bray","euclidean")){
      Distance <- vegan::vegdist(otu_table(phylo1), method = metric) %>% as.matrix()
      
    }
  }
  
  
  pcoa = cmdscale(Distance, k = 3, eig = T)
  
  mds1 = as.data.frame(pcoa$points)
  colnames(mds1) = c("pc1","pc2","pc3")
  data_ord=merge(mds1, meta1, by = "row.names")
  row.names(data_ord) = data_ord$Row.names
  
  
  results = list(Distance=Distance,phylo1=phylo1,data_ord=data_ord,pcoa=pcoa)
  return(results)
}



#### plot time based on the size ,shape should be   hollow,if use time to gradient color, no factor
##plot mean value of each time each group
pcoa_plot_size=function(df, pcoa, metric, text,color,
                        shape,shape_sub={},var1={},var1_sub={},var2={},var2_sub={},var3={},var3_sub={},
                        color_g=8,title_name={},size={}){
  
  df1=sub_select(df,var1 = var1,var1_sub = var1_sub,var2 = var2,var2_sub = var2_sub)
  
  
  
  df1[,shape] <- as.factor(df1[,shape])
  
  if(!is.null(shape_sub)){
    shape_value=shape_sub
  }else{
    shape_value=seq(21,21+length(unique(df1[,shape])))
  }
  
  
  ### mean and order based on time
  df1= dd_rowmean(data1=df1,var1 = "pc1",var2 = "pc2",
                  group1 = color,group2=size)
  
  df1 <- df1[order(as.numeric(as.matrix(df1[,size]))),]
  
  
  
  df1[,color] <- as.factor(df1[,color])
  
  df1[,size]=as.numeric(df1[,size])
  
  # ggplot(df1, aes( pc1, pc2)) +
  ggplot(df1, aes( pc1_mean, pc2_mean)) +
    geom_path(aes(pc1_mean,pc2_mean,group=df1[,color],color=df1[,color]),
              size=1.5,alpha=0.5)+
    # geom_path(aes(pc1_mean,pc2_mean,group=df1[,color]),
    #           color="black",size=1.5,alpha=0.5)+
    geom_point(aes(fill =df1[,  color],size=as.numeric(as.factor(df1[,size]))),
               color="black",shape=21) +
    annotate("text",x = df1$pc1_mean, y =  df1$pc2_mean,
             label=df1[,text],size=6)+
    
    scale_size_continuous(range=c(1,12))+
    scale_shape_manual(values = shape_value) + 
    labs(color = color,shape= shape,fill=color,size=size) + # change the name of the legend's title
    # guides(size = FALSE) + ## omit the legend of size
    xlab(paste("PCoA1 ", round(100*as.numeric(pcoa$eig[1] / sum(pcoa$eig)), 1), "%", sep = "")) +
    ylab(paste("PCoA2 ", round(100*as.numeric(pcoa$eig[2] / sum(pcoa$eig)), 1), "%", sep = "")) +
    ggtitle(paste(metric,title_name,sep = " ")) + theme(plot.title = element_text(hjust = 1, size = 10)) +
    guides(fill = guide_legend(override.aes = list(size = 4)))  +
    
    theme_bw() + theme(panel.grid = element_blank(), legend.text = element_text(size = 15),legend.title = element_text(size=20),
                       axis.text = element_text(size=15),axis.title=element_text(size=20)) 
  
  
  
}


### get the day0 distance to average day0 based on group , only used for donor2acceptor and id , 
####### group_sub could not contain values not in group
########  otu_table rows refers to sampleid
distance_base_product = function(data, meta,var1,var1_sub, group,group_sub,metric, trefile,
                                 mouse_id="label_id",time=var1){
  if(is.data.frame(data)){
    #3# col names for sample id, add the mean row before, should not contain the mean cols
    meta2=sub_select(meta,var1=var1,var1_sub = var1_sub)
    otu_table=data[,rownames(meta2)]
    
  }else{
    library(phyloseq)
    
    phylo2 <- phylo1_get(phylo_ls = data,meta,var1,var1_sub,
                         var2=group,var2_sub=group_sub)
    
    otu_table = t(otu_table(phylo2, taxa_are_rows = TRUE)) %>% as.data.frame() ### col for sample id
    meta2 <- meta[colnames(otu_table),]
    
    
  }
  res1= add_group_mean(counts=otu_table,meta2,time,group)
  otu_table=res1$counts %>% as.data.frame()
  
  colSums(otu_table)
  
  meta3=res1$meta
  
  
  meta3$mix <- paste(meta3$label_id,meta3[,group],meta3[,time],sep = "_")
  base_mean = meta3[is.na(str_match(rownames(meta3),"sample")),]  ##### select the mean row
  
  
  if(metric=="unifrac"){  ### otu_table should be col for sample id here 
    Distance = rbiom::unifrac(as.matrix(otu_table),weighted=FALSE,phy_tree(data))  %>% as.matrix()
    #col for  sample id 
  }else if(metric=="wunifrac"){
    Distance = rbiom::unifrac(as.matrix(otu_table),weighted=TRUE,phy_tree(data))  %>% as.matrix()
  }else if (metric=="bray"){ # row for sample id 
    Distance <- vegan::vegdist(t(otu_table), method = metric)  %>% as.matrix()
    # pcoa = cmdscale(Distance, k = 3, eig = T)
  }
  
  ### select other sample and the mean sample, mean as cols, and sample for row
  Distance2 <- Distance[!is.na(str_match(rownames(Distance),"sample")),
                        is.na(str_match(colnames(Distance),"sample"))]  #%>% as.matrix()
  
  ###change the colnames and the rownames matching the meta2$mix
  if(is.vector(Distance2)){## colnames for mean, rownames for sample id 
    Distance2 = Distance2%>% as.data.frame()
    colnames(Distance2) = rownames(base_mean)
  }
  
  
  # rownames(Distance2) <- meta3[rownames(Distance2),"mix"]
  
  count=0
  #####
  for (i in group_sub) {  ## colnames for mean, rownames for sample id 
    h <- Distance2[!is.na(str_match(rownames(Distance2),"sample")),
                   !is.na(str_match(colnames(Distance2),i))] %>% as.data.frame()
    
    colnames(h) <- "value"
    rownames(h) <- meta3[rownames(Distance2),"mix"]
    count=count+1
    
    if(count==1){
      Distance3 <- h
    }else{
      Distance3 <- rbind(Distance3,h)
    }
    
  }
  
  
  
  return(Distance3)
}


### time_pos for time position ,alway the third position in the character ,all the time to base time for distance 
### each mouse time point distance to its own
distance_refine = function(pcoa_distance,meta,mouse_id,group,time,time_pos=3,phylo1={},
                           time_base=0){
  
  # meta1 <- meta[rownames(otu_table(phylo1)),]
  meta1 <- meta[rownames(pcoa_distance),]
  meta1$mix <- paste(meta1[,mouse_id],meta1[,group],meta1[,time],sep = "_")
  
  rownames(pcoa_distance) <- meta1$mix
  colnames(pcoa_distance) <- meta1$mix
  
  count=0
  for (i in unique(meta1[,mouse_id])) {
    # dis1 <- pcoa_distance[na.omit(str_match(colnames(pcoa_distance),paste("^",i,".*",sep = ""))),
    #                       na.omit(str_match(colnames(pcoa_distance),paste("^",i,".*",sep = "")))]
    dis1 <- pcoa_distance[na.omit(str_match(colnames(pcoa_distance),paste("^",i,"_.*",sep = ""))),
                          na.omit(str_match(colnames(pcoa_distance),paste("^",i,"_.*",sep = "")))]
    # rownames(dis1) <- as.numeric(as.matrix(rownames(h)))
    if(class(dis1)=="numeric")
      next
    ### change to time
    rownames(dis1) <- sapply(str_split(rownames(dis1),"_" ),"[",time_pos)
    ## clear the col to day 0, remain the row with different time
    h <- dis1[,na.omit(str_match(colnames(dis1),paste(".*_",time_base,sep = "")))] %>%as.data.frame()
    rownames(h)=rownames(dis1)
    colnames(h) <- na.omit(str_match(colnames(dis1),paste(".*_",time_base,sep = "")))
    ###some of the smaple has fewer time , deal with it in advanced
    # if(isTRUE(i=="J3")){
    #   h["14",] <- NA
    # }
    h$order <- as.numeric(as.matrix(rownames(h))) ## sort the time
    h <- h[order(h$order),]
    
    if(count==0){
      dis_all <- h
      count = count +1
    }else {
      dis_all <- merge(dis_all,h,by="order",all = TRUE) # use merge, not cbind, which easily make error, some will miss time 
    }
    
    
  }
  rownames(dis_all) <- dis_all$order
  dis_all <- dis_all[,-1]
  dis_all <- t(dis_all) %>% as.data.frame()
  dis_all$id <- rownames(dis_all)
  dis_all[,group]<- sapply(str_split(rownames(dis_all),"_" ),"[",2)
  
  dis_melt <- melt(dis_all,id.vars = c("id",group))
  
  
  return(dis_melt)
}



########## each mouse time point distance to many group's mean value,setting the O2,OO,YO base line to OD
###########  get the distance O2,OO,YO to OD, YY,OY to YD
distance_refine3 = function(pcoa_distance,meta,mouse_id,group,time,time_pos=3,
                            time_base=0){
  
  # meta1 <- meta[rownames(otu_table(phylo1)),]
  meta1 <- meta[rownames(pcoa_distance),]
  meta1[is.na(str_match(meta1$mix,"mean")),"mix"] <- 
    paste(meta1[is.na(str_match(meta1$mix,"mean")),mouse_id],meta1[is.na(str_match(meta1$mix,"mean")),group],meta1[is.na(str_match(meta1$mix,"mean")),time],sep = "_")
  
  rownames(pcoa_distance) <- meta1$mix
  colnames(pcoa_distance) <- meta1$mix
  
  count=0
  for (i in na.omit(unique(meta1[,group]))) {
    if(i %in% c("O2","OO","YO")){
      base_group="OD"
    }else if (i %in% c("YY","OY")){
      base_group="YD"
    }
    # dis1 <- pcoa_distance[na.omit(str_match(colnames(pcoa_distance),paste("^",i,".*",sep = ""))),
    #                       na.omit(str_match(colnames(pcoa_distance),paste("^",i,".*",sep = "")))]
    dis1 <- data.frame(value=pcoa_distance[na.omit(str_match(colnames(pcoa_distance),paste(".*",i,"_.*",sep = ""))),
                                           na.omit(str_match(colnames(pcoa_distance),paste(".*",base_group,"_.*",sep = "")))])
    dis1[,group]<- i 
    dis1[,"distanceTo"] = base_group
    
    # rownames(dis1) <- as.numeric(as.matrix(rownames(h)))
    if(class(dis1)=="numeric")
      next
    
    dis1[,time] <-  sapply(str_split(rownames(dis1),"_"), '[',time_pos)  
    dis1 <- na.omit(dis1)
    
    dis1$order <- as.numeric(dis1[,time]) ## sort the time
    dis1 <- dis1[order(dis1$order),]
    
    if(count==0){
      dis_all <- dis1
      count = count +1
    }else {
      dis_all <- rbind(dis_all,dis1) # use merge, not cbind, which easily make error, some will miss time 
    }
    
    
  }
  dis_all$id <- rownames(dis_all)
  
  
  
  return(dis_all)
}

####  O2,OO,YO distance to YD
distance_refine3_2 = function(pcoa_distance,phylo1,meta,mouse_id,group,time,time_pos=3,
                              time_base=0){
  
  # meta1 <- meta[rownames(otu_table(phylo1)),]
  meta1 <- meta[rownames(pcoa_distance),]
  meta1[is.na(str_match(meta1$mix,"mean")),"mix"] <- 
    paste(meta1[is.na(str_match(meta1$mix,"mean")),mouse_id],meta1[is.na(str_match(meta1$mix,"mean")),group],meta1[is.na(str_match(meta1$mix,"mean")),time],sep = "_")
  
  rownames(pcoa_distance) <- meta1$mix
  colnames(pcoa_distance) <- meta1$mix
  
  count=0
  for (i in na.omit(unique(meta1[,group]))) {
    if(i %in% c("O2","OO","YO")){
      base_group="YD"
    }else if (i %in% c("YY","OY")){
      base_group="OD"
    }
    # dis1 <- pcoa_distance[na.omit(str_match(colnames(pcoa_distance),paste("^",i,".*",sep = ""))),
    #                       na.omit(str_match(colnames(pcoa_distance),paste("^",i,".*",sep = "")))]
    dis1 <- data.frame(value=pcoa_distance[na.omit(str_match(colnames(pcoa_distance),paste(".*",i,"_.*",sep = ""))),
                                           na.omit(str_match(colnames(pcoa_distance),paste(".*",base_group,"_.*",sep = "")))])
    dis1[,group]<- i 
    dis1[,"distanceTo"] = base_group
    
    # rownames(dis1) <- as.numeric(as.matrix(rownames(h)))
    if(class(dis1)=="numeric")
      next
    
    dis1[,time] <-  sapply(str_split(rownames(dis1),"_"), '[',time_pos)  
    dis1 <- na.omit(dis1)
    
    dis1$order <- as.numeric(dis1[,time]) ## sort the time
    dis1 <- dis1[order(dis1$order),]
    
    if(count==0){
      dis_all <- dis1
      count = count +1
    }else {
      dis_all <- rbind(dis_all,dis1) # use merge, not cbind, which easily make error, some will miss time 
    }
    
    
  }
  dis_all$id <- rownames(dis_all)
  
  
  
  return(dis_all)
}


###### pheatmap -----------
##### containing mix in meta1, mix could be used for ,counts could be only DEG counts 
#### no cluster when average, mix contain groups and time
#### DEG is the gene names wich written in the row,meta1 contain the mix for time and group
##### no need to  select counts, but based on meta1, counts should col names for sample
pheatmap_plot = function(counts,meta,DEG={},cluster_cols=TRUE,cluster_rows=TRUE,title_name="",
                         average=FALSE,normalized=TRUE,scale="row",nchar=50,Annotation=FALSE,
                         row_colors=c("red","yellow","blue","purple","orange","black","pink","brown")){
  ## whether need to 
  if("mix" %in% colnames(meta) ){
    meta$mix=as.factor(meta$mix)
  }else{
    print("ERR, no mix var in the meta dataframe")
  }
  if(isTRUE(normalized)){
    counts_heatmap <- round(sweep(counts,2,colSums(counts)/min(colSums(counts)),"/"),digits = 0)
  }else{
    counts_heatmap <- counts
  }
  # colnames(counts_heatmap)==rownames(meta)
  counts_heatmap = counts_heatmap[,rownames(meta)]
  ### keep the order 
  if(is.null(DEG)){
    counts_heatmap1 <- counts_heatmap 
  }else{
    counts_heatmap1 <-  na.omit(counts_heatmap[rownames(DEG),rownames(meta)])
  }
  if(isTRUE(Annotation)){ ### annotation should be the same order as counts data,used for cols
    
    Anno = data.frame(row_anno=factor(meta$mix)) ## rownames should be the same as counts data
    rownames(Anno)=rownames(meta)
    row_colors=row_colors
    names(row_colors)=levels(Anno$row_anno) 
    ann_colors=list(row_anno=row_colors) ## names of colors should the same as Anno, like row_anno
    
  }else{
    Anno=NULL
    ann_colors=NULL
  }
  
  
  if(isFALSE(average)){
    gaps_col = matrix(summary(meta$mix),nrow=length(levels(meta$mix)),ncol=2)
    
    for (i in 1:nrow(gaps_col)) {
      gaps_col[i,2] =sum(gaps_col[1:i,1])
      
    }
    gaps_col <- gaps_col[1:(nrow(gaps_col)-1),2]
    # gaps_col=c(5,10,14,19)
    
    meta=meta[order(meta$mix),]
    
    rownames(counts_heatmap1) <- str_sub(rownames(counts_heatmap1),1,nchar)
    counts_heatmap1=counts_heatmap1[,rownames(meta)]
    
    if(scale=="row"){
      scale_rows = function(x){
        m = apply(x, 1, mean, na.rm = T)
        s = apply(x, 1, sd, na.rm = T)
        return((x - m) / s)
      }
      
      max_value= max(scale_rows(counts_heatmap1))
      min_value=min(scale_rows(counts_heatmap1))
      legend_breaks=c(min_value,max_value)
    }else{
      legend_breaks={}
    }
    
    
    
    pheatmap::pheatmap(as.matrix(counts_heatmap1), color = colorRampPalette(c("green", "black", "red"))(50),
                       cluster_cols = cluster_cols,  cluster_rows = cluster_rows,scale = scale,angle_col=0,
                       main = title_name,fontsize_col = 10,gaps_col = gaps_col,
                       annotation_col = Anno,annotation_colors = ann_colors,legend_breaks = legend_breaks,
                       legend_labels = c("min","max"),
                       display_numbers = FALSE)
    
  }else{  ##for average
    counts_diff1 <- t(counts_heatmap1) 
    counts_diff1 <- as.data.frame(apply(counts_diff1, 2,as.numeric))
    
    counts_diff1 <- cbind(counts_diff1,meta$mix) %>% as.data.frame()
    
    #### mean each group and plot heatmap
    h <- rowsum(counts_diff1[,1:(ncol(counts_diff1)-1)],group = counts_diff1[,ncol(counts_diff1)])
    summary(as.factor(counts_diff1[,ncol(counts_diff1)]))
    h1 <- h/tabulate(as.factor(counts_diff1[,ncol(counts_diff1)]))
    
    pheatmap::pheatmap(as.matrix(na.omit(t(h1))), color = colorRampPalette(c("green", "black", "red"))(50),
                       cluster_cols = FALSE,  cluster_rows = TRUE,scale = scale,angle_col=0,main = title_name,
                       fontsize_col = 20)
  }
  
  
}


#######  ########################## RNA analsyis ---------------------------

DESeq2_analysis = function(counts,meta,condition="condition",
                           level_detail=c("ctrl","abx"),level_ref="ctrl"){
  ## level_ref for reference
  ddsMat <- DESeq2::DESeqDataSetFromMatrix(countData =counts,
                                           colData = meta ,design = formula(paste("~",condition)))
  #### filter less than 10 counts in rows
  # ddsMat <- ddsMat[rowSums(counts(ddsMat)) >= 10,]
  
  colData(ddsMat)[,condition] = factor(colData(ddsMat)[,condition],levels = level_detail)
  colData(ddsMat)[,condition] <- relevel(colData(ddsMat)[,condition],ref = level_ref)
  
  ddsMat <- DESeq2::DESeq(ddsMat)
  return(ddsMat)
}

###get the data from ddsMat
ddsMat_diff_get = function(ddsMat,keep_signif=TRUE,pvalue=0.05){
  dresult <-DESeq2::results(ddsMat)
  ##### select some condition to find different genes
  # dresult <- results(ddsMat,name = "condition_ctrl_vs_abx")
  # dresult <- results(ddsMat,contrast=c("condition","age"))
  ddsMat_data <- data.frame(dresult@listData,row.names = dresult@rownames)  #%>% na.omit()  
  ddsMat_data$threshold = as.factor(ifelse(ddsMat_data$pvalue < pvalue & abs(ddsMat_data$log2FoldChange) >= 1,
                                           ifelse(ddsMat_data$log2FoldChange>1 ,'Up','Down'),'NoSignifi'))
  
  if(isTRUE(keep_signif)){
    ddsMat_diff <-  ddsMat_data[ddsMat_data$threshold != "NoSignifi",]
  }else{
    ddsMat_diff <-  ddsMat_data
  }
  ddsMat_diff=ddsMat_diff[!is.na(ddsMat_diff$threshold),]
  
  ddsMat_diff$id <- rownames(ddsMat_diff)
  return(ddsMat_diff)
}


#### kegg  ,diff_all with id data

kegg_get=function(diff_all,DEG_threshold={},threshold=0.05,
                  count_threshold=1,enrich=TRUE){
  ##### get pathway levels of KEGG
  kegg_level <- read.csv("F:/R file/DaiLab/pathway/kegg_refine2020-08.csv")
  kegg_level <- unique(kegg_level[,c("level1","level2","level3","path_ko")])
  
  if(is.null(DEG_threshold)){
    diff_all1 <- diff_all
  }else{
    diff_all1 <- diff_all[diff_all$threshold.DESeq2==DEG_threshold,]
    # diff_all1 <- diff_all[diff_all$threshold.DESeq==DEG_threshold|
    #                         diff_all$threshold.DESeq2==DEG_threshold,]
  }
  ensids<- as.character(diff_all1$id)
  # ensids<- as.character(rownames(counts[counts$rate>0.0001,]))
  ensids2 <- bitr(ensids,fromType = "SYMBOL",toType = "ENTREZID",
                  OrgDb = org.Mm.eg.db)
  if(isTRUE(enrich)){
    kegg2 <- enrichKEGG(ensids2$ENTREZID, organism = "mmu",
                        pvalueCutoff = threshold)
    
    kegg_rich <- summary(kegg2)
    #http://www.genome.jp/kegg/catalog/org_list.html 
    kegg_rich2 <- kegg_rich[kegg_rich$Count>count_threshold,]
    kegg_rich2 <- unique(merge(kegg_rich2,kegg_level[c("level1","level2","level3","path_ko")],
                               by.x="Description",by.y="level3",all.x=TRUE))
    
    # intersect(kegg_rich2$Description,kegg_level$level3)
    
    
    
    
    return(kegg_rich2)
  }else{
    ensids2 <- merge(ensids2,kegg_level,by.x="ENTREZID",by.y="entrezid",all.x=TRUE)
    
    
    return(ensids2)
  }
  
}



#### volcano

plot_volcano = function(df,x="log2FoldChange",y="padj",y_threshold=0.05,
                        colour="threshold",title_name={},log2FoldChange_threshold=3){
  ## with all data ,but plot the threshold , log2 threshold to label with text
  library(ggrepel)
  ddsMat_diff <- df
  ddsMat_diff$id <- rownames(ddsMat_diff)
  
  # up=ddsMat_diff[ddsMat_diff$threshold=="Up",c("padj","log2FoldChange")]
  # down=ddsMat_diff[ddsMat_diff$threshold=="Down",c("padj","log2FoldChange")]
  
  
  
  #### FDR p < y_threshold,0.05
  ddsMat_diff=ddsMat_diff[!is.na(ddsMat_diff[,y]),]  # remove the NA in the y
  ddsMat_diff[ddsMat_diff[,y]>y_threshold,colour]="NoSignifi"
  
  
  # ddsMat_diff2 <- subset(ddsMat_diff, 
  #                        ddsMat_diff[,y] < y_threshold & abs(ddsMat_diff[,x]) >= log2FoldChange_threshold)
  # 
  # 
  ggplot(data = ddsMat_diff, 
         aes(x = ddsMat_diff[,x], y = -log10(ddsMat_diff[,y]),
             colour=ddsMat_diff[,colour],label =ddsMat_diff$id)) +
    geom_point(alpha=0.4, size=3.5) +
    scale_color_manual(values=c("blue", "gray","red")) +
    labs(title = title_name,x=x,y="-log10 pvalue",colour=colour) +
    ###annotate some of the genes
    # ggrepel::geom_text_repel(data = ddsMat_diff2,aes(x = ddsMat_diff2[,x], y = -log10(ddsMat_diff2[,y]),
    #                                                  colour=ddsMat_diff2[,colour],
    #                                                  label = ddsMat_diff2$id),size = 3,
    #                          box.padding = unit(0.5, "lines"),
    #                          point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = TRUE )+
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0) +
    
    # geom_hline(yintercept = -log10(y_threshold),color="red",size=1.5)+
    # geom_vline(xintercept = 1,color="blue",size=1.5)+
    # geom_vline(xintercept = -1,color="blue",size=1.5)+
    # geom_curve(aes(ymin=2,ymax=6,xmin=1,xmax=3))+
    # geom_polygon(data=data.frame(x=c(min(up$log2FoldChange),max(up$log2FoldChange),
    #                                  max(up$log2FoldChange),min(up$log2FoldChange)),
    #                              y=-log10(c(min(up$padj),min(up$padj),
    #                                  max(up$padj),max(up$padj)))),aes(x=x,y=y,label=y),color="red",alpha=0.3,fill="red")+
    # geom_polygon(data=data.frame(x=c(min(down$log2FoldChange),max(down$log2FoldChange),
    #                                  max(down$log2FoldChange),min(down$log2FoldChange)),
  #                              y=-log10(c(min(down$padj),min(down$padj),
  #                                         max(down$padj),max(down$padj)))),aes(x=x,y=y,label=y),color="blue",alpha=0.3,fill="blue")+
  theme_bw() + theme(panel.grid = element_blank(), legend.text = element_text(size = 15),legend.title = element_text(size=20),
                     axis.text = element_text(size=15),axis.title=element_text(size=20)) 
  
  
}



### for upset
upset_get=function(DEG_res,pvar="padj",p_threshold=1){
  for (i in 1:length(DEG_res)) {
    df1=DEG_res[[i]] 
    df1=df1[df1[,pvar]<p_threshold,]
    
    df1[,names(DEG_res[i])]=1
    df1=df1[,c("id",names(DEG_res[i]))]
    if(i==1){
      df_all=df1
    }else{
      df_all=full_join(df_all,df1)
    }
    
    
  }
  rownames(df_all)=df_all$id
  df_all=subset(df_all,select=-id) %>% as.matrix()
  df_all[is.na(df_all)]=0
  df_all=as.data.frame(df_all)
  return(df_all)
}
