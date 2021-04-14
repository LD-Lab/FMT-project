######  all figure in a paper order
setwd("F:/R file/DaiLab/wangying/FMT_project")
####### change the tree file , and use rbiom::unifrac to replace the phyloseq
## change color from red to orange.

file_name="figure_20210415"

#### rna---------

# setwd("F:/R file/DaiLab/wangying/rna_seq")
library("GenomicAlignments")
library(dplyr)

### import data for rna ------
counts <- read.csv("F:/R file/DaiLab/wangying/FMT_project/rna_seq/all_counts.csv",row.names = 1)
colSums(counts)

## filter some of the gene id, which is not necessary,
filter_counts =read.csv("F:/R file/DaiLab/wangying/FMT_project/rna_seq/rna_filter.csv")


counts=counts[setdiff(rownames(counts),filter_counts$id),]
counts=counts[is.na(str_match(rownames(counts),".*Rik|Gm.*|.*-ps.*|Ig*")),]

colSums(counts)

######## filter the genes with less than 5 reads in at least one sample
counts=counts[rowSums(counts>5)>1,]

# rowSums(counts>5)

colnames(counts)<- sapply(str_split(colnames(counts),"_"),"[",1 )

meta <- read.csv("F:/R file/DaiLab/wangying/FMT_project/rna_seq/metadata.csv",row.names = 1)

meta$mix <- paste(meta$treatment,meta$treatment_time,sep = "-")

meta=meta[-which(rownames(meta) %in% c("ZL30","ZL28")),]




###pvalue =0.05, make DEG---------
deseq2_res=function(meta,counts,var1,var1_sub,var2={},var2_sub={},var3={},var3_sub={},pvalue= 0.05,
                    keep_diff=TRUE,condition="treatment",level_detail=c("ctrl","high"),level_ref="ctrl"){
  
  meta1 <- sub_select(meta,var1 = var1,var1_sub = var1_sub,var2 = var2,
                      var2_sub = var2_sub,var3=var3,var3_sub=var3_sub)
  
  counts1 <- counts[,rownames(meta1)]
  
  ddsMat= DESeq2_analysis(counts=counts1,meta=meta1,condition=condition,
                          level_detail=level_detail,level_ref=level_ref)
  
  
  
  ddsMat_diff <- ddsMat_diff_get(ddsMat,keep_signif = keep_diff,pvalue=pvalue) ## padj =0.05
  # ddsMat_diff1 <- ddsMat_diff[ddsMat_diff$padj< Padj,]
  return(ddsMat_diff)
}


# h=rowSums(counts1) %>% as.matrix()


DEG_res=list()
## gut 
DEG_res[["gut-18-1"]]=deseq2_res(meta=meta,counts=counts,var1 = "age",var1_sub = 18,var2 = "mix",
                                 var2_sub = c("ctrl-1","high-1"))



DEG_res[["gut-18-56"]]=deseq2_res(meta=meta,counts=counts,var1 = "age",var1_sub = 18,var2 = "mix",
                                  var2_sub = c("ctrl-1","high-56"),condition = "mix",
                                  level_detail = c("ctrl-1","high-56"),level_ref = "ctrl-1")

DEG_res[["gut-18-OO"]]=deseq2_res(meta=meta,counts=counts,var1 = "age",var1_sub = 18,var2 = "mix",
                                  var2_sub = c("ctrl-1","OO-56"),condition = "mix",
                                  level_detail = c("ctrl-1","OO-56"),level_ref = "ctrl-1")

DEG_res[["gut-18-YO"]]=deseq2_res(meta=meta,counts=counts,var1 = "age",var1_sub = 18,var2 = "mix",
                                  var2_sub = c("ctrl-1","YO-56"),condition = "mix",
                                  level_detail = c("ctrl-1","YO-56"),level_ref = "ctrl-1")



save(DEG_res,file=paste(file_name,"/DEG_data0415.RData",sep = ""))

# h=DEG_res$`gut-18-56`
# 
# h=h[setdiff(rownames(h),filter_counts$id),]
# h=h[is.na(str_match(rownames(h),".*Rik|Gm.*|.*-ps.*|Ig*")),]

# write.csv(h,file = paste(file_name,"/gut_20 months day1 to ctrl DEGs.csv",sep = ""))

#### all data , padj set to 0.05 for Up or Down description
DEG_res_all=list()

DEG_res_all[["gut-18-1"]]=deseq2_res(meta=meta,counts=counts,var1 = "age",var1_sub = 18,var2 = "mix",
                                 var2_sub = c("ctrl-1","high-1"),keep_diff = FALSE)

DEG_res_all[["gut-18-56"]]=deseq2_res(meta=meta,counts=counts,var1 = "age",var1_sub = 18,var2 = "mix",
                                  var2_sub = c("ctrl-1","high-56"),condition = "mix",
                                  level_detail = c("ctrl-1","high-56"),level_ref = "ctrl-1",keep_diff = FALSE)

DEG_res_all[["gut-18-OO"]]=deseq2_res(meta=meta,counts=counts,var1 = "age",var1_sub = 18,var2 = "mix",
                                  var2_sub = c("ctrl-1","OO-56"),condition = "mix",
                                  level_detail = c("ctrl-1","OO-56"),level_ref = "ctrl-1",keep_diff = FALSE)
DEG_res_all[["gut-18-YO"]]=deseq2_res(meta=meta,counts=counts,var1 = "age",var1_sub = 18,var2 = "mix",
                                  var2_sub = c("ctrl-1","YO-56"),condition = "mix",
                                  level_detail = c("ctrl-1","YO-56"),level_ref = "ctrl-1",keep_diff = FALSE)




save(DEG_res_all,file=paste(file_name,"/DEG_data_all0415.RData",sep = ""))




# h=DEG_res_all$`gut-18-56`
# h1=counts[,rownames(meta[meta$mix=="high-56"|meta$mix=="ctrl-1",])]
# 
# h2=merge(h,h1,by="row.names")
# write.csv(h2,"ctrl1-high56.csv")

### DEGs upset
# gene overlap upset---------------------------



###### based on pvalue, not the padj

load(paste(file_name,"/DEG_data0415.RData",sep = ""))

h=DEG_res$`gut-18-56`


p_threshold=0.05

upset_df=upset_get(DEG_res,pvar = "pvalue",p_threshold = p_threshold)

upset_df1=upset_df[,!is.na(str_match(colnames(upset_df),"18"))]
upset_df1=upset_df1[,c("gut-18-1","gut-18-56","gut-18-OO","gut-18-YO")]

p=  UpSetR::upset(upset_df1,nsets = 10,nintersects = 60,
                  keep.order = T,sets = rev(c("gut-18-1","gut-18-56","gut-18-OO","gut-18-YO")),
                  mainbar.y.label = paste("18-month DEG Intersection Size padj <",p_threshold),
                  sets.x.label = "DEG Set Size", set_size.show = T,
                  text.scale = c(2,2,2,2,2,2),
                  set_size.scale_max = max(colSums(upset_df1))+200,
                  decreasing = c(FALSE,FALSE))
p


# jpeg(paste(file_name,"/padj ",p_adj,"-18 month only gut DEG overlap.png",sep = ""),width=600,height = 600)
# print(p)
# dev.off()



########* gene venn ---------

library(VennDiagram)

load(paste(file_name,"/DEG_data0415.RData",sep = ""))
### set padj for small numbers of DEGs
p_threshold=0.05
pvar = "pvalue"
upset_df=upset_get(DEG_res,pvar = pvar,p_threshold = p_threshold)

# write.csv(upset_df,paste(file_name,"/venn_data.csv",sep = ""))

venn_df=apply(upset_df,2,function(x) rep(names(x),x) )







venn.diagram(
  x=list("gut-20-56"=as.vector(venn_df$`gut-18-56`),
         "gut-20-OO"=as.vector(venn_df$`gut-18-OO`),
         "gut-20-YO"=as.vector(venn_df$`gut-18-YO`)),
  fill=c(color_db[c("OO_DSS","OO","YO"),"value"]),
  alpha=0.5,
  filename=paste(file_name,"/20 months gut only three venn,",pvar,"-",p_threshold,".png",sep = ""))



### DEGs analysis -----------

load(paste(file_name,"/DEG_data0414.RData",sep = ""))


### select pvalue <0.05
h=DEG_res$`gut-18-1`

rna_high1_ctrl=DEG_res$`gut-18-1` 
sum(rna_high1_ctrl$pvalue)
# rna_high1_ctrl1 <- rna_high1_ctrl[rna_high1_ctrl$pvalue<0.05,]

rna_high1_ctrl1 <- arrange(rna_high1_ctrl,pvalue) %>% head(.,n=30)

## fold change
rna_high1_ctrl1 <- arrange(rna_high1_ctrl,log2FoldChange)

rna_high1_ctrl1= rbind(head(rna_high1_ctrl1,n=20),
                       tail(rna_high1_ctrl1,n=20))

meta_1=meta
meta_1$mix = paste(meta$donor2acceptor,meta$treatment_time,sep = "-")

meta2=sub_select(meta_1,var1="mix",var1_sub = c("O0-1","O2-1","O2-56","OO-56","YO-56"))







#######* plot the heatmap based on DEG from ctrl and high day 1 in all groups-------

p=pheatmap_plot(counts=counts,meta2,DEG=rna_high1_ctrl1,cluster_cols = FALSE, 
                # scale="none",
                title_name = paste("pvalue < 0.05","RNA DEG based on 18_ctrl1-high1"),Annotation = TRUE,
                row_colors = color_db[c("O0","O2","O3","OO","YO"),"value"])

p

# pdf(file=paste(file_name,"20 months RNA DEG heatmap.pdf",sep="/"), height = 8, width = 12) #for pdf
# print(p)
# dev.off()








######## volcano ---------
###* FDR  p < 0.01 as the heatmap--------
load(paste(file_name,"/DEG_data_all0415.RData",sep = ""))



padj=0.01

for (i in 1:length(DEG_res_all)) {
  DEG_df=DEG_res_all[[i]]
  
  DEG_df1=DEG_df[is.na(str_match(DEG_df$id,"^[0-9]+")),]
  
  p=plot_volcano(df=DEG_df1,y = "padj",y_threshold=padj,log2FoldChange_threshold=3 )+
    labs(title=names(DEG_res_all[i]))+xlim(-6,6)
  

  ggsave(paste(file_name,"/",names(DEG_res_all[i]),"-",padj,"--canon.png",sep = ""),p,width = 6,height = 6)
  # pdf(file=paste(file_name,"/",names(DEG_res[i]),"--canon.pdf",sep = ""), height = 8, width = 8) #for pdf
  # print(p)
  # dev.off()
  # 
}


