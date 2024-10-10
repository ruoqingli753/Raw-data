if (T) {
  dir.create("scripts")
  dir.create("files")
  dir.create("figures")
  dir.create("00_origin_datas/GEO",recursive = T)
  dir.create("00_origin_datas/TCGA")
  dir.create("00_pre.data",recursive = T)
  dir.create("data",recursive = T)
  
}

library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)
options(expressions = 5e5)

mg_merge_plot <- function(...,ncol = NULL, nrow = NULL,
         labels = NULL, label.x = 0, label.y = 1, hjust = -0.5,
         vjust = 1.5, font.label = list(size = 14, color = "black", face ="bold"
                                        , family = NULL), align = c("none", "h", "v", "hv")
         ,widths = 1, heights = 1, legend = NULL, common.legend = FALSE){
  ml=list(...)
  if(length(ml)==1){
    if(is.list(ml[[1]])){
      ml=ml[[1]]
    }
  }
  gal=ggpubr::ggarrange(plotlist = ml, ncol = ncol, nrow = nrow
                        ,labels = labels,label.x = label.x, label.y = label.y
                        , hjust = hjust,
                        vjust = vjust, font.label = font.label, align = align
                        ,widths = widths, heights = heights, legend = legend
                        , common.legend = common.legend
  )
  return(gal)
}
mg_surv_pROC <- function(time,status,score,mks=c(1,3,5)){
  #time=g.os
  #status=g.ev
  #score=as.numeric(cpm.score)
  #cx=coxRun(data.frame(time,status,score))
  #if(cx[1]<=1){
  #  score=-1*score
  #}
  roc.tm=mg_surv_pROC(time,status,score,mks)
  print('roc.tm')
  print((roc.tm))
  library(survival)
  library(ggplot2)
  mks=mg_predict_time_ymd(time,mks)
  print(mks)  
  ROC.DSST=timeROC::timeROC(T=time,
                            delta=status
                            ,marker=score,
                            cause=1,weighting="marginal",
                            times=mks,
                            iid=TRUE)
  print(ROC.DSST)
  mks=mks[which(!is.na(ROC.DSST$AUC)&ROC.DSST$AUC>0)]
  print(mks)
  if(length(mks)>0){
    if(max(ROC.DSST$AUC)<0.5){
      score=-1*score
    }
    ROC.DSST=timeROC::timeROC(T=time,
                              delta=status
                              ,marker=score,
                              cause=1,weighting="marginal",
                              times=mks,
                              iid=TRUE)
    print(ROC.DSST$times)
    if(max(ROC.DSST$times)<20){
      lb=paste0(ROC.DSST$times,'-Years')
    }else if(max(ROC.DSST$times)<365){
      lb=paste0(round(ROC.DSST$times/12,0),'-Years')
    }else{
      lb=paste0(round(ROC.DSST$times/365,0),'-Years')
    }
    
    lbs=paste0(lb,',AUC=',round(ROC.DSST$AUC,2),',95%CI(',paste0(round(confint(ROC.DSST,level = 0.95,na.rm=T)$CI_AUC[,1]/100,2),'-',
                                                                 round(confint(ROC.DSST,level = 0.95,na.rm=T)$CI_AUC[,2]/100,2)),')')
    #roc.tm=ROC.DSST$times[which(ROC.DSST$times>0)]
    
    #p.dat=rbind()
    #for(i in which(ROC.DSST$times>0)){
    #los=lowess(ROC.DSST$FP[,i], y=ROC.DSST$TP[,i], f = 1/3, iter = 100)
    #los$x=c(0,los$x,1)
    #los$y=c(0,los$y,1)
    # p.dat=rbind(p.dat,data.frame(los$x, y=los$y,rep(lbs[i],length(los$y)),stringsAsFactors = F))
    #}
    
    p.dat=rbind()
    print(length(roc.tm))
    for(i in 1:length(roc.tm)){
      #print(i)
      r1=roc.tm[[i]]
      x1=1-r1$specificities
      y1=r1$sensitivities
      #print(cbind(1-r1$specificities,r1$sensitivities))
      nx1=unique(x1)
      ny1=c()
      for(x in unique(x1)){
        x.inds=which(x1==x)
        if(length(x.inds)>0&x<0.5){
          ny1=c(ny1,min(y1[x.inds]))
        }else if(length(x.inds)>0){
          ny1=c(ny1,max(y1[x.inds]))
        }else{
          ny1=c(ny1,y1[x.inds][1])
        }
      }
      #print(cbind(nx1,ny1))
      p.dat=rbind(p.dat,data.frame(x=nx1, y=ny1,rep(lbs[i],length(nx1)),stringsAsFactors = F))
    }
    colnames(p.dat)=c('V1','V2','Type')
    p.dat=as.data.frame(p.dat)
    
    p1=ggplot(p.dat, aes(x=V1,y=V2, fill=Type))
    p1=p1+geom_line(aes(colour=Type),lwd=1.1)+theme_bw()+xlab('False positive fraction')+ylab('True positive fraction') 
    #p1=p1+stat_smooth(aes(colour=Type),se = FALSE, size = 1)+theme_bw()+xlab('False positive fraction')+ylab('True positive fraction') 
    
    # p1=p1+theme(axis.text.y=element_text(family="Times",face="plain"),axis.text.x=element_text(family="Times",face="plain")
    #             ,axis.title.x=element_text(family="Times",face="plain"),axis.title.y=element_text(family="Times",face="plain")
    #             ,plot.title=element_blank()
    #             ,plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "inches")
    #             ,legend.position=c(1,0)
    #             ,legend.justification=c(1,0)
    #             ,legend.background = element_rect(fill = NA, colour = NA)
    #             ,legend.title = element_text(family="Times",face="plain")
    #             ,legend.text = element_text(family="Times",face="plain"))
    return(p1)
  }else{
    return(mg_getplot_bank('No data plot by ROC!'))
  }
}
bioForest=function(rt=null,col){
  #
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  
  #
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}
my_volcano=function(dat,p_cutoff=0.05,fc_cutoff=1,col=c("red","blue","black"),
                    ylab='-log10 (adj.PVal)',xlab='log2 (FoldChange)',leg.pos='right'){
  degs_dat=dat$DEG
  degs_dat$type=factor(ifelse(degs_dat$adj.P.Val<p_cutoff & abs(degs_dat$logFC) > fc_cutoff, 
                              ifelse(degs_dat$logFC> fc_cutoff ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))
  p=ggplot(degs_dat,aes(x=logFC,y=-log10(adj.P.Val),color=type))+
    geom_point()+
    scale_color_manual(values=col)+#
    # geom_text_repel(
    #   data = tcga.diff$DEG[tcga.diff$DEG$adj.P.Val<p_fit & abs(tcga.diff$DEG$logFC)>fc_fit,],
    #   #aes(label = Gene),
    #   size = 3,
    #   segment.color = "black", show.legend = FALSE )+#
    theme_bw()+#
    theme(
      legend.title = element_blank(),#
      legend.position = leg.pos,
      text = element_text(family = 'Times')
    )+
    ylab(ylab)+#
    xlab(xlab)+#
    geom_vline(xintercept=c(-fc_cutoff,fc_cutoff),lty=3,col="black",lwd=0.5) +#|FoldChange|>2
    geom_hline(yintercept = -log10(p_cutoff),lty=3,col="black",lwd=0.5)#padj<0.05
  return(p)
}

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(text = element_text(family = 'Times'),panel.grid = element_blank())
}
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  
  ## 
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 
  return(geneList)
}
my_boxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),#test_method='kruskal.test',
                    fill= "Group",label=c("p.format",'p.signif')[1],notch=F,size=10,
                    xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat=ARGs.score[tcga.subtype$Samples,'score']
  # group=tcga.subtype$Cluste
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(data, aes(x=Group, y=value,fill=Group)) +
    geom_boxplot(notch = notch) +
    scale_fill_manual(values = group_cols)+   #
    # if(length(names(table(group)))>2){
    #   test_method=''
    # }
    ggpubr::stat_compare_means(aes(group=Group), label = label, method = test_method)+
    labs(x=xlab, y = ylab, fill = fill) +
    theme_bw()+
    theme(legend.position =legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          text = element_text(family = 'Times',size = size),
          axis.text.x = element_text(size = x.size),
          axis.text.y = element_text(size = y.size)) # 
  return(p)
}
my_violin=function(dat,group,group_cols=ggsci::pal_aaas()(10),#test_method='kruskal.test',
                   fill= "Group",label=c("p.format",'p.signif')[1],
                   xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat = tcga.b.cell$GeneSet,
  # group = tcga.subtype$Cluster
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='t.test'
  }
  p=ggplot(data,aes(x=Group, y=value,fill=Group)) +
    geom_violin(trim = F)+  
    geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
    scale_fill_manual(values = group_cols)+
    theme_classic(base_size = 20)+labs(x=xlab,y=ylab,title = title)+
    ggpubr::stat_compare_means(aes(group=Group), label = label, method =test_method)+
    theme(legend.position = legend.position,axis.text = element_text(color = 'black'),
          title = element_text(size = 12),text = element_text(family = 'Times'),
          axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))
  return(p)
}

my_mutiboxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                        #test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],
                        bw=T,xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,
                        legend.position='top',fill='group',notch=F){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method=c('wilcox.test','t.test')[1]
  }
  p=dat.melt %>%
    ggplot(aes(x=type, y=value,fill=Group)) +
    geom_boxplot(notch = notch) +  
    scale_fill_manual(values =group_cols)+   #
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    #theme_light()+
    # theme_bw()+
    theme_classic()+
    theme(legend.position = legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          text = element_text(family = 'Times',size = size),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust)) # 
  return(p)
}
my_mutiviolin=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                       #test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],
                       bw=T,xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,
                       legend.position='top',fill='group'){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggviolin(dat.melt, x = "type", y = "value", fill = "Group",
             add = "boxplot",palette = group_cols)+
    #scale_fill_manual(values =group_cols)+   #
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    theme(text = element_text(family = 'Times',size = size),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust))
  return(p)
}
my_riskplot=function(cli_dat,cols=c("red","blue"),xlab='Samples',
                     a.ylab="Risk score",b.labs="Survival time(year)",cutoff=0,labs=c('A','B')){
  # cli_dat=tcga.risktype.cli
  cli.dat.order=cli_dat[order(cli_dat$Riskscore),c('OS.time','Status','Riskscore','Risktype')]
  fp_dat=data.frame(Samples=1:nrow(cli_dat),cli.dat.order)
  p1=ggplot(fp_dat,aes(x=Samples,y=Riskscore))+geom_point(aes(color=Risktype))+
    scale_colour_manual(values =cols)+
    theme_bw()+labs(x=xlab,y=a.ylab)+
    geom_hline(yintercept=cutoff,colour="black", linetype="dotted",size=0.8)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p1
  p2=ggplot(fp_dat,aes(x=Samples,y=OS.time))+geom_point(aes(col=Status))+theme_bw()+
    scale_colour_manual(values =cols)+
    labs(x=xlab,y=b.labs)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p2
  p.merge=mg_merge_plot(p1,p2,nrow=2,ncol=1,labels = labs)
  return(p.merge)
}
coxFun <- function(dat){
  library(survival)
  colnames(dat)=c('time','status','gene')
  fmla=as.formula("Surv(time,status)~gene")
  cox=coxph(fmla,data=dat)
  p=summary(cox)[[7]][5]
  result=c(p,summary(cox)[[8]][1],summary(cox)[[8]][3],summary(cox)[[8]][4])
  return(result)
}
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  library(glmnet)
  set.seed(2021)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Partial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=names(sig.coef),lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}

##TCGA#####
#
genecode=read.delim('data/GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]
head(mrna_genecode)

tcga.pancancer.cli=read.xlsx('data/TCGA_pancancer_cli_PMID_29625055.xlsx')
head(tcga.pancancer.cli)
tcga.cli=tcga.pancancer.cli[which(tcga.pancancer.cli$type=='STAD'),]
head(tcga.cli)
tcga.cli=data.frame(Samples=paste0(tcga.cli$bcr_patient_barcode,'-01'),
                    Age=tcga.cli$age_at_initial_pathologic_diagnosis,
                    Gender=tcga.cli$gender,
                    AJCC_stage=tcga.cli$ajcc_pathologic_tumor_stage,
                    # Clinical_Stage=tcga.cli$clinical_stage,
                    Grade=tcga.cli$histological_grade,
                    tcga.cli[,c('OS','OS.time','DSS','DSS.time','DFI','DFI.time','PFI','PFI.time')])
rownames(tcga.cli)=tcga.cli$Samples
head(tcga.cli)
tcga.cli$OS.time
tcga.cli=tcga.cli %>% drop_na(OS.time)
tcga.cli=tcga.cli[tcga.cli$OS.time>0,]
dim(tcga.cli)
#417  13
head(tcga.cli)
fivenum(tcga.cli$Age)
tcga.cli$Age1=ifelse(tcga.cli$Age>67,'>67','<=67')
table(tcga.cli$Gender)
table(tcga.cli$AJCC_stage)
tcga.cli$AJCC_stage[tcga.cli$AJCC_stage=='[Discrepancy]'|tcga.cli$AJCC_stage=='[Not Available]']=NA
tcga.cli$AJCC_stage=gsub('[ABC]','',tcga.cli$AJCC_stage)
tcga.cli$AJCC_stage=gsub('Stage ','',tcga.cli$AJCC_stage)

table(tcga.cli$Grade)
tcga.cli$Grade[tcga.cli$Grade=='GX']=NA


tcga.data=read.delim('00_origin_datas/TCGA/Merge_RNA_seq_TPM.txt',row.names = 1,check.names = F)
tcga.data[1:4,1:4]
table(substr(colnames(tcga.data),14,15))
dim(tcga.data)

sample_T=colnames(tcga.data)[which(as.numeric(substr(colnames(tcga.data),14,15))==1)]#
sample_T=intersect(sample_T,tcga.cli$Samples)
sample_N=colnames(tcga.data)[which(as.numeric(substr(colnames(tcga.data),14,15))==11)]#
tcga_type=data.frame(Samples=c(sample_T,sample_N),Type=rep(c('Tumor','Normal'),c(length(sample_T),length(sample_N))))
rownames(tcga_type)=tcga_type$Samples
table(tcga_type$Type)


length(sample_T)
#353

range(tcga.data)
tcga.data=log2(tcga.data+1)
range(tcga.data)
tcga.data[1:5,1:5]


tcga.exp=tcga.data[intersect(mrna_genecode$SYMBOL,rownames(tcga.data)),sample_T]
dim(tcga.exp)
tcga.cli=tcga.cli[sample_T,]
tcga.exp.all=tcga.data[intersect(mrna_genecode$SYMBOL,rownames(tcga.data)),tcga_type$Samples]
###################################

#GSE6143 HP infection#####
# cagA , babA2 , and vacAs1 triple-positive samples 
GSE6143 <- getGEOExpData('GSE6143')
GSE6143.cli=GSE6143$Sample
GSE6143.df <- GSE6143$Exp$GPL193_862_Data_col3
dim(GSE6143.df)
range(GSE6143.df)
GSE6143.df <- log2(GSE6143.df+1)
GSE6143.exp <- exp_probe2symbol_v2(datExpr = GSE6143.df,GPL = 'GPL193')
GSE6143.cli <- data.frame(samples=GSE6143.cli$Acc,type=GSE6143.cli$Title)
GSE6143.cli$group <- ifelse(grepl("Hp", GSE6143.cli$type), "Control", "HP")
table(GSE6143.cli$group)
GSE6143.exp <- GSE6143.exp[,GSE6143.cli$samples]
dim(GSE6143.exp)
# 855  24

#GSE5081HP infection######
#
GSE5081 <- getGEOExpData('GSE5081')
GSE5081.cli=GSE5081$Sample
GSE5081.df <- GSE5081$Exp$GPL570_54675_Data_col1
dim(GSE5081.df)
range(GSE5081.df)
GSE5081.df <- log2(GSE5081.df+1)
GSE5081.exp <- exp_probe2symbol_v2(datExpr = GSE5081.df,GPL = 'GPL570')
GSE5081.cli <- data.frame(samples=GSE5081.cli$Acc,type=GSE5081.cli$Title)
table(GSE5081.cli$type)
GSE5081.cli$'group'[grepl("Gastric biopsy HP- ER-,",GSE5081.cli$type)] <- "Control" 
GSE5081.cli$'group'[grepl("Gastric biopsy HP\\+", GSE5081.cli$type)] <- "HP"
GSE5081.cli <- subset(GSE5081.cli, group %in% c("Control", "HP"))
table(GSE5081.cli$group)
GSE5081.exp <- GSE5081.exp[,GSE5081.cli$samples]
dim(GSE5081.exp)
# 20549    24

#GSE27411 HP infection##############
GSE27411 <- getGEOExpData('GSE27411')
GSE27411.cli=GSE27411$Sample
GSE27411.df <- GSE27411$Exp$GPL6255_20589_Data_col1
dim(GSE27411.df)
range(GSE27411.df)
#GSE27411.df <- log2(GSE27411.df+1)
GSE27411.exp <- exp_probe2symbol_v2(datExpr = GSE27411.df,GPL = 'GPL6255')
GSE27411.cli <- data.frame(samples=GSE27411.cli$Acc,type=GSE27411.cli$`disease state`)
table(GSE27411.cli$type)
# atrophy     infected non-infected 
# 6            6            6 
GSE27411.cli <- GSE27411.cli[GSE27411.cli$type %in% c('infected','non-infected'),]
GSE27411.cli$group <- ifelse(GSE27411.cli$type=='infected','HP',"Control")
GSE27411.exp <- GSE27411.exp[,GSE27411.cli$samples]
dim(GSE27411.exp)
#16261    12
#GSE60662 HP infection###############
#
GSE60662 <- getGEOExpData('GSE60662')
GSE60662.cli=GSE60662$Sample
GSE60662.df <- GSE60662$Exp$GPL13497_34127_Data_col1
dim(GSE60662.df)
range(GSE60662.df)
GSE60662.df <- log2(GSE60662.df+1)
range(GSE60662.exp)
GSE60662.exp <- exp_probe2symbol_v2(datExpr = GSE60662.df,GPL = 'GPL13497')
GSE60662.cli <- data.frame(samples=GSE60662.cli$Acc,type=GSE60662.cli$Title)
table(GSE60662.cli$type)
GSE60662.cli$group <- "HP"
GSE60662.cli$group[grepl("control",GSE60662.cli$type)] <- "Control"
table(GSE60662.cli$group)
# Control      HP 
# 4      12 
GSE60662.exp <- GSE60662.exp[,GSE60662.cli$samples]
dim(GSE60662.exp)
#21840    16
# GSE28541 ##########
library(stringr)
GSE28541_cli <- openxlsx::read.xlsx('00_origin_datas/GEO/GSE91061/clinic_datas.xlsx', sheet = 5)
# View(GSE28541_cli)
# rm(GSE28541_cli)
GSE28541_cli <- GSE28541_cli[, c("GEO.ID", "OS.m", "Death.(1=yes,.0=no)")]
colnames(GSE28541_cli) <- c("samples", "OS.time", "OS")
GSE28541_cli$OS.time <- GSE28541_cli$OS.time * 30
GSE28541_cli <- GSE28541_cli[GSE28541_cli$OS.time > 0, ]

GSE28541_cli1 <- getGEOSampleData('GSE28541')
GSE28541_cli1 <- GSE28541_cli1[, c("Acc", "sample_id")]


# GSE28541_cli1$Title <- str_split_fixed(GSE28541_cli1$Title, '_', 5)[, 4]

GSE28541_cli <- merge(GSE28541_cli, GSE28541_cli1,
                      by.x = 'samples', by.y = 'sample_id')
GSE28541_cli$samples <- GSE28541_cli$Acc

GSE28541 <- getGEOExpData('GSE28541')


GSE28541_anno <- GSE28541$Anno$GPL13376

GSE28541_exp <- GSE28541$Exp$GPL13376_48701_Data_col1
dim(GSE28541_exp)
GSE28541_exp[1:5, 1:5]
boxplot(GSE28541_exp[, 1:5])

GSE28541_exp <- exp_probe2symbol_v2(GSE28541_exp,
                                    GPL = 'GPL13376')
# GSE28541_exp <- log2(GSE28541_exp + 1)
boxplot(GSE28541_exp[, 1:15])


intersect(GSE28541_cli$samples, colnames(GSE28541_exp))

rownames(GSE28541_cli) <- GSE28541_cli$samples



GSE28541_exp <- GSE28541_exp[, intersect(colnames(GSE28541_exp), GSE28541_cli$samples)]
GSE28541_cli <- GSE28541_cli[intersect(colnames(GSE28541_exp), GSE28541_cli$samples), ]


# GSE13861 ##########
library(stringr)
GSE13861_cli <- openxlsx::read.xlsx('00_origin_datas/GEO/GSE91061/clinic_datas.xlsx', sheet = 2)
View(GSE13861_cli)
# rm(GSE13861_cli)
GSE13861_cli <- GSE13861_cli[, c("Patients_ID", "OS.m", "Death.(1=yes,.0=no)")]
colnames(GSE13861_cli) <- c("samples", "OS.time", "OS")
GSE13861_cli$OS.time <- GSE13861_cli$OS.time * 30
GSE13861_cli <- GSE13861_cli[GSE13861_cli$OS.time > 0, ]

GSE13861_cli1 <- getGEOSampleData('GSE13861')
GSE13861_cli1 <- GSE13861_cli1[, c("Acc", "Title")]


GSE13861_cli1$Title <- str_split_fixed(GSE13861_cli1$Title, '_', 5)[, 4]

GSE13861_cli <- merge(GSE13861_cli, GSE13861_cli1,
                      by.x = 'samples', by.y = 'Title')
GSE13861_cli$samples <- GSE13861_cli$Acc

GSE13861 <- getGEOExpData('GSE13861')

GSE13861_anno <- GSE13861$Anno$GPL6884

GSE13861_exp <- GSE13861$Exp$GPL6884_48803_Data_col1
dim(GSE13861_exp)
GSE13861_exp[1:5, 1:5]
boxplot(GSE13861_exp[, 1:5])

GSE13861_exp <- exp_probe2symbol_v2(GSE13861_exp,GPL ='GPL6884' )
range(GSE13861_exp)
# GSE13861_exp <- log2(GSE13861_exp + 1)
boxplot(GSE13861_exp[, 1:15])


intersect(GSE13861_cli$samples, colnames(GSE13861_exp))

rownames(GSE13861_cli) <- GSE13861_cli$samples

GSE13861_exp <- GSE13861_exp[, intersect(colnames(GSE13861_exp), GSE13861_cli$samples)]
GSE13861_cli <- GSE13861_cli[intersect(colnames(GSE13861_exp), GSE13861_cli$samples), ]

dim(GSE13861_cli)
dim(GSE13861_exp)
##GSE38749#####
GSE38749 <- getGEOExpData('GSE38749' )
GSE38749.exp <- GSE38749$Exp$GPL570_54675_Data_col1
GSE38749.cli <- GSE38749$Sample
GSE38749.cli <- data.frame(Samples)

GSE38749.cli=data.frame(Samples=GSE38749.cli$Acc,
                        Age=GSE38749.cli$`age (y)`,
                        Gender=GSE38749.cli$gender,
                        Stage=GSE38749.cli$`tumor stage (ajcc)`,
                        OS.time=GSE38749.cli$`time (months) overall survival`/12*365,
                        OS=GSE38749.cli$status)


head(GSE38749.cli)
rownames(GSE38749.cli)=GSE38749.cli$Samples
GSE38749.cli$OS.time
GSE38749.cli=GSE38749.cli%>%drop_na(OS.time)
GSE38749.cli=GSE38749.cli[GSE38749.cli$OS.time>0,]
GSE38749.cli$OS=ifelse(GSE38749.cli$OS=='alive',0,1	)
dim(GSE38749.cli)

GSE38749.exp <- exp_probe2symbol_v2(datExpr = GSE38749.exp,GPL = 'GPL570')
range(GSE38749.exp)

##GSE15459############
GSE15459=readRDS('00_pre.data/GSE15459.rds')
GSE15459.cli=readxl::read_xls('00_pre.data/GSE15459_outcome.xls')
head(GSE15459.cli)
GSE15459.cli=data.frame(Samples=GSE15459.cli$`GSM ID`,
                        Age=GSE15459.cli$Age_at_surgery,
                        Gender=GSE15459.cli$Gender,
                        Stage=GSE15459.cli$Stage,
                        OS.time=GSE15459.cli$`Overall.Survival (Months)**`/12*365,
                        OS=GSE15459.cli$`Outcome (1=dead)`)
head(GSE15459.cli)
rownames(GSE15459.cli)=GSE15459.cli$Samples
GSE15459.cli$OS.time
GSE15459.cli=GSE15459.cli%>%drop_na(OS.time)
GSE15459.cli=GSE15459.cli[GSE15459.cli$OS.time>0,]
GSE15459.cli$OS
dim(GSE15459.cli)


GSE15459.df=exprs(GSE15459)
GSE15459.df[1:5,1:5]
GSE15459.exp=exp_probe2symbol_v2(datExpr = GSE15459.df,GPL = 'GPL570')
rm(GSE15459.df)
range(GSE15459.exp)
GSE15459.exp=log2(GSE15459.exp+1)
GSE15459.exp=GSE15459.exp[,GSE15459.cli$Samples]
dim(GSE15459.exp)



##GSE62254########
GSE62254=readRDS('00_pre.data/GSE62254.rds')
GSE62254.cli=read.xlsx('00_pre.data/GSE66229_cli.xlsx')
head(GSE62254.cli)
GSE62254.cli=data.frame(Samples=GSE62254.cli$GEO_ID,
                        Age=GSE62254.cli$age,Gender=GSE62254.cli$sex,
                        Stage=GSE62254.cli$Stage,
                        OS=GSE62254.cli$Death,OS.time=GSE62254.cli$OS.m/12*365)
rownames(GSE62254.cli)=GSE62254.cli$Samples
GSE62254.cli$OS.time


GSE62254.df=exprs(GSE62254)
GSE62254.df[1:5,1:5]
GSE62254.exp=exp_probe2symbol_v2(datExpr = GSE62254.df,GPL = 'GPL570')
range(GSE62254.exp)
GSE62254.exp=GSE62254.exp[,GSE62254.cli$Samples]
dim(GSE62254.exp)

##GSE66229#######
GSE66229=readRDS('00_pre.data/GSE66229.rds')

GSE66229.pheno=pData(GSE66229)
GSE66229.pheno=data.frame(Samples=GSE66229.pheno$geo_accession,tissue=GSE66229.pheno$`tissue:ch1`)
rownames(GSE66229.pheno)=GSE66229.pheno$Samples
head(GSE66229.pheno)

GSE66229.df=exprs(GSE66229)
GSE66229.df[1:5,1:5]
GSE66229.exp=exp_probe2symbol_v2(datExpr = GSE66229.df,GPL = 'GPL570')
rm(GSE66229.df)
range(GSE66229.exp)
dim(GSE66229.exp)

GSE66229.cli=read.xlsx('00_pre.data/GSE66229_cli.xlsx')
head(GSE66229.cli)
GSE66229.cli=data.frame(Samples=GSE66229.cli$GEO_ID,
                        Age=GSE66229.cli$age,Gender=GSE66229.cli$sex,
                        Stage=GSE66229.cli$Stage,
                        OS=GSE66229.cli$Death,OS.time=GSE66229.cli$OS.m/12*365)
rownames(GSE66229.cli)=GSE66229.cli$Samples
GSE66229.cli$OS.time
dim(GSE66229.cli)



##GSE84437#####
GSE84437=readRDS('00_pre.data/GSE84437.rds')
GSE84437.cli=pData(GSE84437)
GSE84437.cli=data.frame(Samples=GSE84437.cli$geo_accession,
                        Age=GSE84437.cli$`age:ch1`,Gender=GSE84437.cli$`Sex:ch1`,
                        OS=GSE84437.cli$`death:ch1`,OS.time=GSE84437.cli$`duration overall survival:ch1`)
rownames(GSE84437.cli)=GSE84437.cli$Samples
GSE84437.cli$OS.time
GSE84437.cli=GSE84437.cli%>%drop_na(OS.time)
GSE84437.cli=GSE84437.cli[GSE84437.cli$OS.time>0,]
GSE84437.cli=crbind2DataFrame(GSE84437.cli)
GSE84437.cli$OS.time=GSE84437.cli$OS.time/12*365
head(GSE84437.cli)
dim(GSE84437.cli)
#

GSE84437.df=exprs(GSE84437)
GSE84437.df[1:15,1:15]
dim(GSE84437.df)

GSE84437.exp=exp_probe2symbol_v2(datExpr = GSE84437.df,GPL = 'GPL6947')
rm(GSE84437.df)
GSE84437.exp[1:5,1:5]
boxplot(GSE84437.exp[,1:10])
GSE84437.exp=log2(GSE84437.exp+1)
GSE84437.exp=GSE84437.exp[,GSE84437.cli$Samples]
dim(GSE84437.exp)

#GSE26253####
# GSE26253 <- getGEOExpData('GSE26253')
 #save(GSE26253,file='00_pre.data/GEO/GSE26253.RData')
load('00_pre.data/GEO/GSE26253.RData')
GSE26253_cli=GSE26253$Sample
head(GSE26253_cli)
GSE26253_cli=data.frame(Sampels=GSE26253_cli$Acc,
                        Stage=GSE26253_cli$`pathological stage`,
                        OS=GSE26253_cli$`status (0=non-recurrence, 1=recurrence)`,
                        OS.time=GSE26253_cli$`recurrence free survival time (month)`*30)
rownames(GSE26253_cli)=GSE26253_cli$Sampels
GSE26253_exp=exp_probe2symbol_v2(datExpr =GSE26253$Exp$GPL8432_17418_Data_col1,
                                 GPL ='GPL8432')
GSE26253_exp=GSE26253_exp[,rownames(GSE26253_cli)]
range(GSE26253_exp)
dim(GSE26253_exp)
#13391   432

#GSE26942  #####
#GSE26942 <- getGEOExpData('GSE26942')
#save(GSE26942,file='00_pre.data/GEO/GSE26942.RData')
load('00_pre.data/GEO/GSE26942.RData')
GSE26942_cli=read.delim('00_pre.data/GEO/GSE26942_cli.txt',header = T,sep='\t')
GSE26942_sample=read.table('00_pre.data/GEO/GSE26942_samples_code.txt',header = T,sep = '\t')
GSE26942_cli=merge(GSE26942_cli,GSE26942_sample,by.x = 'Patients_ID',by.y = 'Title')
GSE26942_cli=data.frame(Samples=GSE26942_cli$Accession,
                        Sex=GSE26942_cli$Sex,
                        Age=GSE26942_cli$Age,
                        Stage=GSE26942_cli$AJCC.stage,
                        M.Stage=GSE26942_cli$M.stage,
                        OS=GSE26942_cli$Death..1.yes..0.no.,
                        OS.time=ceiling(GSE26942_cli$OS.m*30),
                        RFS=GSE26942_cli$Recurrence..1.yes..0.no.,
                        RFS.time=ceiling(GSE26942_cli$RFS.m*30))
rownames(GSE26942_cli)=GSE26942_cli$Samples
GSE26942_exp<-exp_probe2symbol_v2(GSE26942$Exp$GPL6947_48803_Data_col1,
                                  GSE26942$Anno$GPL6947[,c(1,14)],method = 'mean')
GSE26942_exp=GSE26942_exp[,GSE26942_cli$Samples]
GSE26942_exp <- GSE26942_exp[stringr::str_split_fixed(rownames(GSE26942_exp), ' /// ', 3)[, 2] == '', ]
GSE26942_exp=crbind2DataFrame(GSE26942_exp)
GSE26942_exp=na.omit(GSE26942_exp)
range(GSE26942_exp)
dim(GSE26942_exp)
#25127    93



# GSE66229##########
#GSE66229 <- getGEOExpData('GSE66229')
#save(GSE66229,file='00_pre.data/GEO/GSE66229.RData')
load('00_pre.data/GEO/GSE66229.RData')
GSE66229_cli=openxlsx::read.xlsx('00_pre.data/GEO/GSE66229_cli.xlsx',sheet = 1)
GSE66229_cli <- as.data.frame(GSE66229_cli)
GSE66229_cli=data.frame(Samples=GSE66229_cli$GEO_ID,
                        OS=GSE66229_cli$Death,
                        OS.time=ceiling(GSE66229_cli$OS.m*30),
                        Sex=GSE66229_cli$sex,
                        T.Stage=GSE66229_cli$T,
                        N.Stage=GSE66229_cli$N,
                        M.Stage=GSE66229_cli$M,
                        Stage=GSE66229_cli$Stage)
GSE66229_cli=crbind2DataFrame(GSE66229_cli)
rownames(GSE66229_cli)=GSE66229_cli$Samples
GSE66229_exp<-exp_probe2symbol_v2(GSE66229$Exp$GPL570_54675_Data_col1,
                                  GSE66229$Anno$GPL570[,c(1,11)],method = 'mean')
GSE66229_exp=GSE66229_exp[,GSE66229_cli$Samples]
GSE66229_exp <- GSE66229_exp[stringr::str_split_fixed(rownames(GSE66229_exp), ' /// ', 3)[, 2] == '', ]
head(GSE66229_cli)
range(GSE66229_exp)
dim(GSE66229_exp)
#21655   300

#01 HP Gene###########
###1.1HP_Gene#######
dir.create("01_HP_Gene")
HP_genes <- read.table("01_HP_Gene/HP_Gene_PMID_36248367.txt",header = T,sep = "\t",check.names = F)
HP_genes <- HP_genes$Genes
setdiff(HP_genes,rownames(tcga.exp))
HP_score <- ssGSEAScore_by_genes(tcga.exp.all,genes =HP_genes )
HP_score <- t(HP_score)
###HP score
p1a <- my_violin(dat =HP_score ,group =tcga_type$Type,ylab = "HP_score" ,group_cols = c('#7469B6','#FFDB5C'))
savePDF(filename = "01_HP_Gene/p1a.pdf",p1a,height = 6,width = 6)
####

###1.2 GEO  ######
####tcga DEG########
tcga.limma <- mg_limma_DEG(tcga.exp.all,group =tcga_type$Type,ulab = "Tumor" ,dlab = "Normal")
tcga.limma$Summary
#######
tcga.degs=tcga.limma$DEG
tcga.degs=tcga.degs[abs(tcga.degs$logFC)>log2(1.5) & tcga.degs$P.Value<0.05,]
dim(tcga.degs)
write.csv(tcga.degs,'01_HP_Gene/tcga.degs.csv')


####GSE6143 DEG########
GSE6143.limma <- mg_limma_DEG(GSE6143.exp,group =GSE6143.cli$group,ulab = "HP" ,dlab = "Control")
GSE6143.limma$Summary
#######
GSE6143.degs=GSE6143.limma$DEG
GSE6143.degs=GSE6143.degs[abs(GSE6143.degs$logFC)>log2(1.5) & GSE6143.degs$P.Value<0.05,]
dim(GSE6143.degs)
write.csv(GSE6143.degs,'01_HP_Gene/GSE6143.degs.csv')


####GSE5081 DEG########
GSE5081.limma <- mg_limma_DEG(GSE5081.exp,group =GSE5081.cli$group,ulab = "HP" ,dlab = "Control")
GSE5081.limma$Summary
# 1.2-fold  1.3-fold  1.5-fold  2-fold   
# p<0.05   "686|824" "610|703" "418|457" "141|181"
# p<0.01   "211|136" "208|123" "173|75"  "84|28"  
# FDR<0.05 "4|0"     "4|0"     "4|0"     "4|0"    
# FDR<0.01 "0|0"     "0|0"     "0|0"     "0|0" 
#######
GSE5081.degs=GSE5081.limma$DEG
GSE5081.degs=GSE5081.degs[abs(GSE5081.degs$logFC)>log2(1.5) & GSE5081.degs$P.Value<0.05,]
dim(GSE5081.degs)
write.csv(GSE5081.degs,'01_HP_Gene/GSE5081.degs.csv')

####GSE27411 DEG########
GSE27411.limma <- mg_limma_DEG(GSE27411.exp,group =GSE27411.cli$group,ulab = "HP" ,dlab = "Control")
GSE27411.limma$Summary
#          1.2-fold   1.3-fold  1.5-fold  2-fold  
# p<0.05   "1111|977" "859|618" "473|215" "144|56"
# p<0.01   "483|289"  "459|248" "326|113" "119|33"
# FDR<0.05 "70|20"    "70|20"   "69|18"   "49|12" 
# FDR<0.01 "16|1"     "16|1"    "16|1"    "15|1"  
#######
GSE27411.degs=GSE27411.limma$DEG
GSE27411.degs=GSE27411.degs[abs(GSE27411.degs$logFC)>log2(1.5) & GSE27411.degs$P.Value<0.05,]
dim(GSE27411.degs)
write.csv(GSE27411.degs,'01_HP_Gene/GSE27411.degs.csv')

####GSE60662 DEG########
GSE60662.limma <- mg_limma_DEG(GSE60662.exp,group =GSE60662.cli$group,ulab = "HP" ,dlab = "Control")
GSE60662.limma$Summary
#         1.2-fold    1.3-fold    1.5-fold    2-fold   
# p<0.05   "1718|2083" "1660|1856" "1480|1212" "954|392"
# p<0.01   "670|631"   "670|614"   "653|467"   "543|178"
# FDR<0.05 "8|6"       "8|6"       "8|6"       "8|6"    
# FDR<0.01 "0|1"       "0|1"       "0|1"       "0|1"  
#######
GSE60662.degs=GSE60662.limma$DEG
GSE60662.degs=GSE60662.degs[abs(GSE60662.degs$logFC)>log2(1.5) & GSE60662.degs$P.Value<0.05,]
dim(GSE60662.degs)
write.csv(GSE60662.degs,'01_HP_Gene/GSE60662.degs.csv')




###1.4 #######
HP.col <- c('#FFA62F','#ACD793')
HP_score.all <- ssGSEAScore_by_genes(tcga.exp,genes =HP_genes )
HP_score.all <- as.data.frame(t(HP_score.all))
HP.cli <- data.frame(Sample=tcga.cli$Samples,HP_score.all=HP_score[tcga.cli$Samples,'GeneSet'],OS=tcga.cli$OS,
                     OS.time=tcga.cli$OS.time)
# HP.data.point <- surv_cutpoint(HP.cli, time = "OS.time", event = "OS",
#                                 variables = 'Riskscore')
# HP.cutoff <- as.numeric(summary(HP.data.point)[1])
# HP.cutoff
HP.cli$type=ifelse(HP.cli$HP_score>median(HP.cli$HP_score),'High','Low')
HP.km=ggplotKMCox(data.frame(HP.cli$OS.time/365,
                             HP.cli$OS,
                             HP.cli$type),
                  palette = HP.col,show_confint = F,title = 'HP cohort')
HP.km

#02.WGCNA####
dir.create('02_WGCNA')
library(WGCNA)
allowWGCNAThreads(nThreads = 36)#
enableWGCNAThreads(nThreads = 36)# 
my_mad <- function(x){mad(x,na.rm = TRUE)} #mad
wgcna_exp=t(tcga.exp)
m.mad <- apply(wgcna_exp,2,my_mad)
dim(tcga.exp)
# 19503   353
#
tpm_T2 <- wgcna_exp[,which(m.mad >max( quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01))]
# dim(tpm_T2)
#  60 13836
#tpm_T2=tcga.exp[which(apply(tcga.exp,1,sd)>0.5),]
#tpm_T2=(2^tpm_T2-1)
range(tpm_T2)

pdf('02_WGCNA/1.pdf',width = 10,height = 10)
tpm_T2.power=mg_wgcna_get_power(tpm_T2)
dev.off()

# minModuleSize = 30,    ##
# mergeCutHeight = 0.25, ##
tpm_T2.power$cutPower
tpm_T2.module=mg_WGCNA_getModule(tpm_T2,
                                 power = tpm_T2.power$cutPower,
                                 deepSplit=2,
                                 mergeCutHeight=0.3,
                                 minModuleSize=50)


table(tpm_T2.module$Modules[,2])
length(table(tpm_T2.module$Modules[,2]))
#9
write.csv(tpm_T2.module$Modules,file = "02_WGCNA/WGCNA_Modules.csv",col.names = T,row.names = T)
pdf('02_WGCNA/2.pdf',height = 5,width = 6)
plotDendroAndColors(tpm_T2.module$Tree, tpm_T2.module$Modules,
                    c("Dynamic Module",'Merged Module'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
#

pdf('02_WGCNA/3.pdf',height = 6,width = 6)
mg_barplot_point(labels = names(table(tpm_T2.module$Modules[,2]))
                 ,values = as.numeric(table(tpm_T2.module$Modules[,2]))
                 ,point_sizes = 2
                 ,point_cols = names(table(tpm_T2.module$Modules[,2]))
                 ,xlab = 'Number of Genes',legend.pos = NULL)
dev.off()

#### 
# Calculate eigengenes
MEs = tpm_T2.module$MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf('02_WGCNA/4.pdf',height = 6,width = 12,onefile = T)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
dev.off()




##choose 
tcga_cli_use <-cbind.data.frame(tcga.cli[tcga.cli$Samples,],HP_score=HP_score[tcga.cli$Samples,])
head(tcga_cli_use)
tcga_cli_use=tcga_cli_use[,-1]
colnames(tcga_cli_use)
tcga_cli_use.part=data.frame(HP_score=tcga_cli_use[tcga.cli$Samples,c(14)])
str(tcga_cli_use.part)
tcga_cli_use.part=sapply(tcga_cli_use.part, function(x)as.numeric(as.factor(x)))


spms=tcga_cli_use.part
MEs_col<-tpm_T2.module$MEs
dim(MEs_col)
modTraitCor = cor(MEs_col[,rownames(MEDiss)[METree$order]]
                  , spms
                  ,use = 'pairwise.complete.obs')
modTraitP = corPvalueStudent(modTraitCor, dim(spms)[1])
textMatrix = paste(signif(modTraitCor, 2), " (", format(modTraitP,scientific =TRUE,digits = 3), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
dim(textMatrix)

rownames(modTraitCor)=gsub("ME","",rownames(modTraitCor))
rownames(textMatrix)=gsub("ME","",rownames(textMatrix))
colnames(modTraitCor)

pdf('02_WGCNA/5.pdf',width = 4,height =6)
labeledHeatmap(Matrix = data.frame(modTraitCor),
               xLabels = colnames(modTraitCor),
               yLabels = rownames(modTraitCor),
               cex.lab = 1,
               ySymbols = colnames(t(modTraitCor)), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = data.frame(textMatrix),
               setStdMargins = FALSE,
               cex.text = 0.8, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#
geneModuleMembership <- signedKME(tpm_T2
                                  , data.frame(tpm_T2.module$MEs)
                                  , outputColumnName = "")
head(geneModuleMembership)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership)
                                           , nrow(tpm_T2.module$MEs)))
#
geneTraitSignificance <- as.data.frame(cor(tpm_T2
                                           , spms
                                           , use = 'pairwise.complete.obs'))

GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance)
                                           , nrow(spms)))

modNames<-colnames(geneModuleMembership)
modNames

module = "purple"
column = match(module, modNames)
column
module1 = "royalblue"
column1 = match(module1, modNames)
column1

moduleGenes<- (tpm_T2.module$Modules[,'mergedColors']==module)
moduleGenes1 <- (tpm_T2.module$Modules[,'mergedColors']==module1)
tcga.wgcna.gene=c(names(which(moduleGenes)),names(which(moduleGenes1)))
length(tcga.wgcna.gene)
#1103
#3######
dir.create('03_GO_KEGG')
wgcna.enrichment=mg_clusterProfiler(genes = tcga.wgcna.gene)
p1=barplot(wgcna.enrichment$KEGG)+ggtitle('KEGG')+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))+theme(text=element_text(family = 'Times'))
p2=barplot(wgcna.enrichment$GO_BP)+ggtitle('Biological Process')+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))+theme(text=element_text(family = 'Times'))
p3=barplot(wgcna.enrichment$GO_CC)+ggtitle('Cell Component')+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))+theme(text=element_text(family = 'Times'))
p4=barplot(wgcna.enrichment$GO_MF)+ggtitle('Molecular Function')+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))+theme(text=element_text(family = 'Times'))

pdf('03_GO_KEGG/enrichment.pdf',height = 15,width = 15,onefile = F)
mg_merge_plot(
mg_merge_plot(p1,p2,ncol=2,widths = c(1.3,1),labels = LETTERS[1:2]),
mg_merge_plot(p3,p4,ncol=2,widths =c(1,1.2),labels = LETTERS[3:4]),nrow = 2)
dev.off()

write.xlsx(list(KEGG=wgcna.enrichment$KEGG@result,BP=wgcna.enrichment$GO_BP@result,
                CC=wgcna.enrichment$GO_CC@result,MF=wgcna.enrichment$GO_MF@result),
           '03_GO_KEGG/wgcna_enrichment_res.xlsx',overwrite = T)

##3.1venny####

######upset######
install.packages("UpSetR")
library(UpSetR)
venny.data <- list(
  GSE60662 = rownames(GSE60662.degs),
  GSE6143 = rownames(GSE6143.degs),
  GSE27411 = rownames(GSE27411.degs),
  GSE5081 = rownames(GSE5081.degs),
  wgcna.module.gene= tcga.wgcna.gene
)

gene_matrix <- stack(venny.data) %>%
  mutate(ind = as.factor(ind), values = as.factor(values)) %>%table() %>%as.matrix()

binary_matrix <- gene_matrix>0
binary_matrix <- 1 * binary_matrix  # TRUE 为 1, FALSE 为 0
binary_matrix <- as.data.frame(binary_matrix )
pdf("03_venny/upset_HP_genes.pdf",height = 4,width =10)
upset(
  binary_matrix,
  nset = 5,  # You have 5 sets in venny.data
  #nintersects = 20,
  order.by = c('degree', 'freq'),
  decreasing = c(TRUE, TRUE),
  mb.ratio = c(0.7, 0.3),
  point.size = 1.8,
  line.size = 1,
  mainbar.y.label = "Intersection size",
  sets.x.label = "Set Size",
  main.bar.color = "#2a83a2",
  sets.bar.color = "#3b7960"
  # queries = list(list(query = intersects,
  #                      params = list("GSE60662", "GSE6143", "GSE27411", "GSE5081", "HP_genes"),
  #                      active = TRUE,
  #                      color = "#d66a35",
  #                      query.name = "Complex Intersection"))
)
dev.off()
######venny#####
library(VennDiagram)
library(RColorBrewer)
library(limma)

venn_ploy <-venn.diagram(
  venny.data,
  filename = NULL,
  fill = brewer.pal(5, "Set1")
)
grid.draw(venn_ploy)
ggsave('03_venny/VENNY_HP_gene.pdf',height = 6,width = 6)


##########
# 
pairwise_intersects <- list(
  intersect(rownames(GSE60662.degs), tcga.wgcna.gene),
  intersect(rownames(GSE6143.degs), tcga.wgcna.gene),
  intersect(rownames(GSE27411.degs), tcga.wgcna.gene),
  intersect(rownames(GSE5081.degs), tcga.wgcna.gene)
)


venny_gene_tcga_wgcna <- Reduce(union,pairwise_intersects)
length(venny_gene_tcga_wgcna)
# 532

#######
dir.create('files/model_select')
select_gene_zscore <- function(dat1, dat2 = NULL, dat3 = NULL, 
                               a = 1, n = 100, 
                               ratio = 0.5, cut_p = 0.05, 
                               years = c(1,3,5)){
  library(timeROC)
  library(survival)
  library(glmnet)
  # library(mosaic)
  sample.index <- data.frame()
  SampleingTime <- c()
  SamplingSeed <- c()
  tra.auc <- c()
  test.auc <- c()
  geo1.auc <- c()
  geo2.auc <- c()
  all.auc <- c()
  tra.p <- c()
  test.p <- c()
  geo1.p <- c()
  geo2.p <- c()
  all.p <- c()
  tcga.dat <- dat1
  geo1.dat <- dat2
  geo2.dat <- dat3
  gene.list <- c()
  GeneNums <- c()
  for (i in a:n) {
    # i=2
    set.seed(i)
    par(mfrow = c(2, 3))
    myd.index <- seq(1,nrow(tcga.dat),1)
    tra.index <- sample(myd.index,size = round(nrow(tcga.dat)*ratio))
    tra.dat <- tcga.dat[tra.index,]
    test.dat <- tcga.dat[-tra.index,]
    write.table(rownames(tra.dat), file = paste0('files/model_select/tra.dat_zscore_', i, '.txt'), sep = '\t', quote = F, row.names = F)
    write.table(rownames(test.dat), file = paste0('files/model_select/test.dat_zscore_', i, '.txt'), sep = '\t', quote = F, row.names = F)
    tra.cox <- t(apply(tra.dat[,3:c(ncol(tra.dat))],2,function(x){
      vl=as.numeric(x)
      tm=tra.dat$OS.time
      ev=tra.dat$OS
      #ev=ifelse(ev=='Alive',0,1)
      dat=data.frame(tm,ev,vl)[which(tm > 0 & !is.na(vl)),]
      return(coxFun(dat))
    }))
    colnames(tra.cox) <- c('p.value','HR','Low 95%CI','High 95%CI')
    tra.cox <- as.data.frame(tra.cox)
    write.table(tra.cox, file = paste0('files/model_select/sig_cox_zscore_', i, '.txt'), sep = '\t', quote = F)
    cut.cox=tra.cox[which(tra.cox[,1]<cut_p),] 
    if (nrow(cut.cox) > 1) {
      print(paste("Processing....: ",i," resampling",sep=""))
      flush.console()
      geneid=rownames(cut.cox)
      geneid=geneid[which(!is.na(geneid))]
      if (length(geneid)>1) {
        sample.index<-c(sample.index,data.frame(tra.index))
        set.seed(i)
        cv.fit <- cv.glmnet(as.matrix(tra.dat[,geneid]), cbind(time=tra.dat$OS.time, 
                                                               status=tra.dat$OS)
                            ,family="cox", nlambda=100, alpha=1) 
        sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
        write.table(sig.coef, file = paste0('files/model_select/sig.coef_zscore_', i, '.txt'), sep = '\t', quote = F)
        if (length(names(sig.coef)>1)) {
          dat1 <- cbind(time=tra.dat$OS.time,
                        status=tra.dat$OS,
                        tra.dat[,match(names(sig.coef),
                                       colnames(tra.dat))])
          fmla <- as.formula(paste0("Surv(time, status) ~"
                                    ,paste0(names(sig.coef),collapse = '+')))
          cox <- coxph(fmla, data =as.data.frame(dat1))
          # lan=coef(cox)
          # print(lan)
          cox1=step(cox, trace = 0)
          lan=coef(cox1)
          write.table(lan, file = paste0('files/model_select/lan_zscore_', i, '.txt'), sep = '\t', quote = F)
          # lan=sig.coef
          final_gene=names(cox1$coefficients)
          # final_gene=names(cox$coefficients)
          GeneNums <-c(GeneNums, length(final_gene))
          
          risk.tra=as.numeric(lan%*%as.matrix(t(tra.dat[,final_gene])))
          # risk.tra=zscore(risk.tra)
          ROC.DSST <- timeROC(T=tra.dat$OS.time,
                              delta=tra.dat$OS,
                              marker=risk.tra,
                              cause=1,
                              weighting="marginal",
                              times=c(365*years[1],365*years[2],365*years[3]),
                              iid=TRUE)
          tra.auc=c(tra.auc,max(ROC.DSST$AUC))
          trap=plotKMCox(data.frame(tra.dat$OS.time,tra.dat$OS,ifelse(risk.tra>=median(risk.tra),'H','L')))
          tra.p=c(tra.p,trap)
          risk.test=as.numeric(lan%*%as.matrix(t(test.dat[,final_gene])))
          # risk.test=zscore(risk.test)
          ROC.DSST1=timeROC(T=test.dat$OS.time,delta=test.dat$OS,marker=risk.test,cause=1,weighting="marginal",
                            times=c(365*years[1],365*years[2],365*years[3]),iid=TRUE)
          test.auc=c(test.auc,max(ROC.DSST1$AUC))
          testp=plotKMCox(data.frame(test.dat$OS.time,test.dat$OS,ifelse(risk.test>=median(risk.test),'H','L')))
          test.p=c(test.p,testp)
          risk.all=as.numeric(lan%*%as.matrix(t(tcga.dat[,final_gene])))
          # risk.all=zscore(risk.all)
          ROC.DSST3=timeROC(T=tcga.dat$OS.time,delta=tcga.dat$OS,marker=risk.all,cause=1,weighting="marginal",
                            times=c(365*years[1],365*years[2],365*years[3]),iid=TRUE)
          all.auc=c(all.auc,max(ROC.DSST3$AUC))
          allp=plotKMCox(data.frame(tcga.dat$OS.time,tcga.dat$OS,ifelse(risk.all>=median(risk.all),'H','L')))
          all.p=c(all.p,allp)
          # final_gene1=as.character(gene.type[final_gene,1])
          
          if (length(intersect(final_gene,colnames(geo1.dat)))==length(final_gene)) {
            risk.geo1=as.numeric(lan%*%as.matrix(t(geo1.dat[,intersect(final_gene,colnames(geo1.dat))])))
            # risk.geo1=zscore(risk.geo1)
            ROC.DSST2=timeROC(T=geo1.dat$OS.time,delta=geo1.dat$OS,marker=risk.geo1,cause=1,weighting="marginal",
                              times=c(365*years[1],365*years[2],365*years[3]),iid=TRUE)
            geo1.auc=c(geo1.auc,max(ROC.DSST2$AUC))
            geop=plotKMCox(data.frame(geo1.dat$OS.time,geo1.dat$OS,ifelse(risk.geo1>=median(risk.geo1),'H','L')))
            geo1.p=c(geo1.p,geop)
          }else{
            geo1.auc=c(geo1.auc,NA)
            geo1.p=c(geo1.p,NA)
          }
          
          if (length(intersect(final_gene,colnames(geo2.dat)))==length(final_gene)) {
            risk.geo2=as.numeric(lan%*%as.matrix(t(geo2.dat[,intersect(final_gene,colnames(geo2.dat))])))
            # risk.geo2=zscore(risk.geo2)
            ROC.DSST4=timeROC(T=geo2.dat$OS.time,delta=geo2.dat$OS,marker=risk.geo2,cause=1,weighting="marginal",
                              times=c(365*years[1],365*years[2],365*years[3]),iid=TRUE)
            geo2.auc=c(geo2.auc,max(ROC.DSST4$AUC))
            geop=plotKMCox(data.frame(geo2.dat$OS.time,geo2.dat$OS,ifelse(risk.geo2>=median(risk.geo2),'H','L')))
            geo2.p=c(geo2.p,geop)
          }else{
            geo2.auc=c(geo2.auc,NA)
            geo2.p=c(geo2.p,NA)
          }
          print(c('AUC',max(ROC.DSST$AUC),max(ROC.DSST1$AUC),max(ROC.DSST3$AUC)))
          print(c('P',trap,testp,allp))
          print(c('num',length(final_gene)))
          print("........")
          SampleingTime=c(SampleingTime,i)
          gene.list=c(gene.list,as.character(final_gene))
          SamplingSeed=c(SamplingSeed,rep(i,length(final_gene)))
        }
      }
    }
  }
  myd.clustering.df=data.frame("SamplingTime"=SampleingTime,
                               "TrainRiskP"=tra.p,"TrainRiskAUC"=tra.auc,
                               "TestRiskP"=test.p,"TestRiskAUC"=test.auc,
                               "TCGARiskP"=all.p,"TCGARiskAUC"=all.auc,
                               "GEO1RiskP"=geo1.p,"GEO1RiskAUC"=geo1.auc,
                               "GEO2RiskP"=geo2.p,"GEO2RiskAUC"=geo2.auc,
                               "GeneNums"=GeneNums)
  sample.index=as.data.frame(sample.index)
  colnames(sample.index)=paste("seed",SampleingTime,sep="")
  myd.clustering.genes=data.frame("SamplingSeed"=SamplingSeed,"Genes"=gene.list)
  return(list(myd.clustering.df,sample.index,myd.clustering.genes))
}

tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[venny_gene_tcga_wgcna, tcga.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))

myd_exp_resampling <- select_gene_zscore(dat1 = tcga_model_data,
                                         # dat2 = silu_model_data,
                                         # dat3 = GSE84437_model_data,
                                         a = 1,
                                         n = 3000,
                                         ratio = 0.7,
                                         cut_p = 0.05,
                                         years = c(1:5))
myd_exp_resampling[[1]]

write.csv(myd_exp_resampling[[1]],file = "04_model/myd_exp_resampling.csv")
#04.#######
dir.create('04_model')
num <- 2297#OK
tra.samples <- rownames(read.delim(paste0('files/model_select/tra.dat_zscore_',num,'.txt'), 
                                   header = T, row.names = 1, stringsAsFactors = F))
test.samples <- rownames(read.delim(paste0('files/model_select/test.dat_zscore_',num,'.txt'),
                                    header = T, row.names = 1, stringsAsFactors = F))
tra.data <- tcga_model_data[tra.samples, ]
dim(tra.data)
test.data <- tcga_model_data[test.samples, ]
dim(test.data)

# write.csv(data.frame(cohort=rep(c('Train','Test'),c(length(tra.samples),length(test.samples))),
#                      tcga.cli.merge[c(tra.samples,test.samples),
#                                     c('Age','Gender','pathologic_T','pathologic_N','pathologic_M','pathologic_stage','OS','OS.time')]),
#           'results/TCGA_clinical.csv',quote = F,row.names = F)


tra.cox=cox_batch(dat = tcga.exp[venny_gene_tcga_wgcna,tra.samples],
                  time = tcga.cli[tra.samples,]$OS.time,event = tcga.cli[tra.samples,]$OS)
tra.cox=na.omit(tra.cox)
head(tra.cox)

rownames(tra.cox)=gsub('-','__',rownames(tra.cox))
p_cutoff=0.05
table(tra.cox$p.value<p_cutoff)
# FALSE  TRUE 
# 515    17 
tra.cox.fit=tra.cox[tra.cox$p.value<p_cutoff,]
tra.cox.fit
write.csv(tra.cox.fit,file = "04_model/tra.cox.fit.csv")

#########lasso
library(glmnet)
set.seed(num)
fit1=glmnet(as.matrix(tra.data[,rownames(tra.cox.fit)])
            ,cbind(time=tra.data$OS.time,
                   status=tra.data$OS)
            ,family="cox",nlambda=100, alpha=1) 

cv.fit<-cv.glmnet(as.matrix(tra.data[,rownames(tra.cox.fit)])
                  ,cbind(time=tra.data$OS.time,
                         status=tra.data$OS)
                  ,family="cox",nlambda=100, alpha=1)
sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
cv.fit$lambda.min
names(sig.coef)

pdf('04_model/LASSO.pdf',height = 5,width = 10,onefile = F)
par(mfrow=c(1,2))
plot(fit1)
plot(cv.fit)
dev.off()

fmla <- as.formula(paste0("Surv(OS.time, OS) ~",paste0(names(sig.coef),collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tra.data))
cox=step(cox)

module.coxforest=ggforest(cox, data = tra.data, 
                          main = "Hazardratio", fontsize =1.0, 
                          noDigits = 2)
module.coxforest
ggsave('04_model/gene_forest.pdf',module.coxforest,height = 4,width = 9)

lan <- coef(cox)
lan
paste0(round(lan, 3), '*', names(lan),collapse = '+')
"0.142*EMB+0.134*CPVL+-0.436*CTLA4+0.329*FAM241A+0.259*CXCR4"
write.csv(data.frame(gene=names(lan),coef=as.numeric(lan)),'04_model/gene_coef.csv',row.names = F)
####
###########
lan.dataframe <- as.data.frame(lan)
lan.dataframe$gene <- rownames(lan.dataframe) 
lan.dataframe$gene <- factor(lan.dataframe$gene,levels = rownames(lan.dataframe)[order(lan.dataframe$lan)])
# 
lan.dataframe$color_group <- ifelse(lan.dataframe$lan > 0, "Positive", "Negative")
library(ggplot2)
# 
p <- ggplot(lan.dataframe, aes(x=gene, y=lan,fill=color_group)) +
  geom_bar(stat="identity") +
  xlab("Gene Name") +
  ylab("Coefficient") +
  ggtitle("Gene Coefficients") +
  coord_flip() +
  scale_fill_manual(values = c("Positive" = "#FC4100", "Negative" = "#00215E")) +
  theme_bw()+
  guides(fill=FALSE)
p1 <- p+geom_text(aes(label=sprintf("%.3f", lan)), hjust=-0.2, size=3, color="black")
ggsave('04_model/gene_ Coefficients.pdf',p1,height = 4,width = 9)

##4.1########
risktype.col=c('#FFDA78',"#003285")
risk.tra=as.numeric(lan%*%as.matrix(t(tra.data[tra.samples,names(lan)])))
tra.risktype.cli=data.frame(tcga.cli[tra.samples,],Riskscore=risk.tra)

#######
tra.data.point <- surv_cutpoint(tra.risktype.cli, time = "OS.time", event = "OS",
                                variables = 'Riskscore')
tra.cutoff <- as.numeric(summary(tra.data.point)[1])
tra.cutoff
tra.risktype.cli$Risktype=ifelse(tra.risktype.cli$Riskscore>tra.cutoff,'High','Low')
#tra.risktype.cli$Risktype=ifelse(tra.risktype.cli$Riskscore>median(risk.tra),'High','Low')
tra.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                               data = tra.risktype.cli),
                  data=tra.risktype.cli,
                  conf.int = T,pval = T,fun = "pct",risk.table =T, size = 0.7,
                  surv.median.line = 'hv',title='Train cohort',
                  linetype = c("solid", "dashed","strata")[1],
                  palette = risktype.col,ggtheme = custom_theme(),
                  legend = c(0.8,0.75), # 
                  legend.title = "Risktype",legend.labs=c('High','Low'))
tra.km=mg_merge_plot(tra.km$plot,tra.km$table,nrow=2,heights = c(2.5,1),align = 'v')
tra.km


tra.roc=ggplotTimeROC(tra.risktype.cli$OS.time,
                      tra.risktype.cli$OS,
                      tra.risktype.cli$Riskscore,mks = c(1,2,3,4,5))
tra.roc



########4.2#######
risk.test=as.numeric(lan%*%as.matrix(t(test.data[test.samples,names(lan)])))
test.risktype.cli=data.frame(tcga.cli[test.samples,],Riskscore=risk.test)
###
test.data.point <- surv_cutpoint(test.risktype.cli, time = "OS.time", event = "OS",
                                 variables = 'Riskscore')
test.cutoff <- as.numeric(summary(test.data.point)[1])
test.cutoff
test.risktype.cli$Risktype=ifelse(test.risktype.cli$Riskscore>test.cutoff,'High','Low')
#test.risktype.cli$Risktype=ifelse(test.risktype.cli$Riskscore>median(risk.test),'High','Low')
test.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = test.risktype.cli),
                   data=test.risktype.cli,
                   conf.int = T,pval = T,fun = "pct",risk.table =T, size = 0.7,
                   surv.median.line = 'hv',title='Test cohort',
                   linetype = c("solid", "dashed","strata")[1],
                   palette = risktype.col,ggtheme = custom_theme(),
                   legend = c(0.8,0.85), # 
                   legend.title = "Risktype",legend.labs=c('High','Low'))
test.km=mg_merge_plot(test.km$plot,test.km$table,nrow=2,heights = c(2.5,1),align = 'v')
test.km


test.roc=ggplotTimeROC(test.risktype.cli$OS.time,
                       test.risktype.cli$OS,
                       test.risktype.cli$Riskscore,mks = c(1,2,3,4,5))
test.roc









##4.3TCGA########
risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga.cli$Samples,names(lan)])))
tcga.risktype.cli=data.frame(tcga.cli,Riskscore=risk.tcga)
########
tcga.data.point <- surv_cutpoint(tcga.risktype.cli, time = "OS.time", event = "OS",
                                 variables = 'Riskscore')
tcga.cutoff <- as.numeric(summary(tcga.data.point)[1])
tcga.cutoff
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>tcga.cutoff,'High','Low')
#tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>median(risk.tcga),'High','Low')
tcga.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga.risktype.cli),
                   data=tcga.risktype.cli,
                   conf.int =T,pval = T,fun = "pct",risk.table =T, size = 0.7,
                   surv.median.line = 'hv',title='TCGA',
                   linetype = c("solid", "dashed","strata")[1],
                   palette = risktype.col,ggtheme = custom_theme(),
                   legend = c(0.8,0.85), # 
                   legend.title = "Risktype",legend.labs=c('High','Low'))
tcga.km=mg_merge_plot(tcga.km$plot,tcga.km$table,nrow=2,heights = c(2.5,1),align = 'v')
tcga.km

# pdf("04_model/all_tcga_roc.pdf",height = 8,width = 8)
# mg_surv_pROC(tcga.risktype.cli$OS.time,
#               tcga.risktype.cli$OS,
#               tcga.risktype.cli$Riskscore,mks = c(1,2,3,4,5))
# dev.off()

tcga.roc=ggplotTimeROC(tcga.risktype.cli$OS.time,
                       tcga.risktype.cli$OS,
                       tcga.risktype.cli$Riskscore,mks = c(1,2,3,4,5))
tcga.roc
# ########4.4 GSE38749验证集##########
# GSE38749_model_data <- data.frame(GSE38749.cli[, c("OS.time", "OS")],
#                                   t(GSE38749.exp[intersect(names(lan),rownames(GSE38749.exp)), GSE38749.cli$Samples]))
# colnames(GSE38749_model_data) <- gsub('-', '_', colnames(GSE38749_model_data))
# 
# risk.GSE38749=as.numeric(lan%*%as.matrix(t(GSE38749_model_data[GSE38749.cli$Samples,names(lan)])))
# 
# GSE38749.risktype.cli=data.frame(GSE38749.cli,Riskscore=risk.GSE38749)
# #GSE38749.risktype.cli$Risktype=ifelse(GSE38749.risktype.cli$Riskscore>median(risk.GSE38749),'High','Low')
# GSE38749.data.point <- surv_cutpoint(GSE38749.risktype.cli, time = "OS.time", event = "OS",
#                                      variables = 'Riskscore')
# GSE38749.cutoff <- as.numeric(summary(GSE38749.data.point)[1])
# GSE38749.cutoff
# GSE38749.risktype.cli$Risktype=ifelse(GSE38749.risktype.cli$Riskscore>GSE38749.cutoff,'High','Low')
# GSE38749.roc=ggplotTimeROC(GSE38749.risktype.cli$OS.time,
#                            GSE38749.risktype.cli$OS,
#                            GSE38749.risktype.cli$Riskscore,mks = c(1,2,3,4))
# GSE38749.roc
# GSE38749.km=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
#                                    data = GSE38749.risktype.cli),
#                        data=GSE38749.risktype.cli,
#                        conf.int = T,pval = T,risk.table = T, 
#                        fun = "pct",size = 1,surv.median.line = 'hv',
#                        title='GSE38749',legend.title='Risktype',
#                        legend.labs = c('High','Low'),
#                        linetype = c("solid", "dashed","strata")[1],
#                        palette = risktype.col,ylab='Overall Survival(OS)',
#                        legend=c(0.85,0.25),#
#                        ggtheme = theme_bw(base_size = 12))
# GSE38749.km=mg_merge_plot(GSE38749.km$plot,GSE38749.km$table,nrow=2,heights = c(3,1),align = 'v')
# GSE38749.km

####TCGA#####
my_mutibarplot=function(df,xlab='group',leg.title='',cols=pal_d3()(10)[7:9]){
  prop.pval=round(chisq.test(df)$p.value,2)#round(-log10(chisq.test(df)$p.value),2)
  if( prop.pval<0.001)
    prop.pval='<0.001'
  df.prop=prop.table(df,margin=2)
  df.prop=reshape2::melt(df.prop)
  colnames(df.prop)<-c("type","group","Percentage")
  df.prop$Percentage<-round(df.prop$Percentage,digits=2)
  p=ggplot(df.prop,aes(x=group,y=Percentage,fill=type))+
    geom_bar(position = "fill",stat="identity")+
    scale_fill_manual(values = cols)+
    xlab(xlab)+labs(fill = leg.title,title = 'Chi-Squared Test',subtitle  =  paste0('pvalue  ',prop.pval))+
    theme_bw()+theme(text=element_text(family = 'Times'),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
  p
  return(p)
}
tcga.risktype.cli$Status <- ifelse(tcga.risktype.cli$OS==0,"alive","dead")
tcga.barplot=my_mutibarplot(df=table(tcga.risktype.cli$Status,tcga.risktype.cli$Risktype),xlab = 'Risktype',leg.title = 'Status')
#####合并图#####
p4=mg_merge_plot(mg_merge_plot(tra.km,test.km,tcga.km,common.legend = T,labels = c('C','D','E'),ncol=3),
                 mg_merge_plot(tra.roc,test.roc,tcga.roc,labels = c('F','G','H'),ncol=3),
                mg_merge_plot(tcga.barplot,GSE38749.km,GSE38749.roc,labels = c('I','J','K'),ncol=3),
                 nrow = 3)
 
  
ggsave('04_model/P4c-k.pdf',height = 15,width = 15)


#05.######
dir.create('05_immu')
###mcp###
tcga.mcp <- immu_MCPcounter(exp =tcga.exp,isTCGA = T)
saveRDS(taga.mcp,file ='05_immu/taga.mcp.RDS')
mg_PlotMutiBoxplot(data = taga.mcp[tcga.risktype.cli$Samples,],group =tcga.risktype.cli$Risktype,
                   legend.pos = 'top',group_cols = risktype.col,add = 'boxplot',test_method = 'wilcox.test',ylab = 'ssgsea Immune Score')
library(tidyverse)
library(ggcor)
library(vegan)
cr=psych::corr.test(x=tcga.risktype.cli[,'Riskscore'],
                    y=tcga.mcp[tcga.risktype.cli$Samples,]
                    ,method = 'spearman')
df=t(rbind(cr$r,cr$p))
colnames(df)=c('r','p.value')
df=data.frame(Riskscore='Riskscore',MCP_count=rownames(df),df)
df
df <- df %>%
  mutate(lty = cut(r, breaks = c(-1, 0, 1),
                   labels = c("r <= 0", "r > 0")),
         col = cut(p.value, breaks = c(0, 0.01, 0.05, 1),
                   labels = c("< 0.01", "< 0.05", ">= 0.05"),
                   right = FALSE, include.lowest = TRUE))
head(df)
corrmat.color=colorRampPalette(c('blue', 'white','red'))(100)

p5a<-quickcor(tcga.mcp[tcga.risktype.cli$Samples,], cor.test = TRUE,type = "lower") + #upper
  geom_square(data = get_data(p.value < 0.05, type = "lower")) + 
  anno_link(df, mapping = aes(colour = col,
                              size = abs(r),
                              linetype = lty),diag.label = TRUE) +
  scale_fill_gradient2n(colours = corrmat.color) +
  remove_x_axis()+
  scale_size_area(max_size = 2) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  guides(
    fill = guide_colourbar(title = "corr", order = 3),
    colour = guide_legend(title = "spearman's p", order = 2),
    size = guide_legend(title = "spearman's r", order = 1),
    linetype = "none")
p5a
ggsave('05_immu/p5a.pdf',height = 6,width = 6)
# ICGs####
tcga.icgs=immu_ICGs(tcga.exp)
colnames(tcga.icgs)
icg.dat.RS=cbind(tcga.risktype.cli$Riskscore
                 ,tcga_model_data[tcga.risktype.cli$Samples,names(lan)]
                 ,tcga.icgs[tcga.risktype.cli$Samples,c('CTLA4','PDCD1','PDCD1LG2','LGALS9','CD80','CD28','HAVCR2')])
#c('CTLA4','PDCD1','PDCD1LG2','LGALS9','CD80','CD28','HAVCR2')
colnames(icg.dat.RS)[1]='Riskcsore'

icg_cor_res <- Hmisc::rcorr(as.matrix(icg.dat.RS),type = 'spearman')
icg_cor_res$P[is.na(icg_cor_res$P)] <- 0
icg_cor_res.p=icg_cor_res$P
icg_cor_res.p[1:5,1:5]
icg_cor_res.p<-ifelse(icg_cor_res.p<0.0001,'****',
                      ifelse(icg_cor_res.p<0.001,'***', 
                             ifelse(icg_cor_res.p<0.01,'**',
                                    ifelse(icg_cor_res.p<0.05,'*',''))))

pdf('05_immu/p5b.pdf',height = 6,width = 7,onefile = F)
p5b=pheatmap(icg_cor_res$r[-c(1:6),c(names(lan),'Riskcsore')],
         color = circlize::colorRamp2(c(-1, 0, 1), c('#3B4992FF', 'white', '#EE0000FF')),
         main="Heatmap", # 
         display_numbers = icg_cor_res.p[-c(1:6),c(names(lan),'Riskcsore')], # 
         cluster_cols = F, # 
         cluster_rows = F,
         show_rownames = T, #
         show_colnames = T,
         fontsize_row = 12, # 
         fontsize_col = 16)
dev.off()


ggsave('05_Risktype.immune/p5ab.pdf',p5ab,height = 6,width = 15)


#06####
HP.cox=cox_batch(dat = tcga.exp[HP_genes,tcga.cli$Samples],
                 time = tcga.cli$OS.time,event = tcga.cli$OS)
HP.cox=na.omit(HP.cox)
head(HP.cox)

rownames(HP.cox)=gsub('-','__',rownames(HP.cox))
p_cutoff=0.05
table(HP.cox$p.value<p_cutoff)
# FALSE  TRUE 
#34     3 
HP.cox.fit=HP.cox[HP.cox$p.value<p_cutoff,]
HP.cox.fit
write.csv(HP.cox.fit,file = "singlecell/HP.cox.fit.csv")

