library("keras")
library("GenomicRanges")
library("rtracklayer")
library("data.table")
library("openxlsx")
library("pheatmap")
library("RColorBrewer")
library("grid")
library('ggplotify')
library("cowplot")
library("pROC")
library("ggplot2")
library("ggrepel")
library("plyr")
library("dplyr")

format_sample_table=function(file, main.class="Disease", sub.class="IDH1", read.col="Filter_reads", sort=F){
    if(is.data.frame(file)){
        sample.df=file
    }else{
        sample.df=read.xlsx(file)
    }
    sample.df$Label=paste(sample.df[, main.class],sep="_")
    sample.df$Label=factor(sample.df$Label, levels=sort(unique(sample.df$Label)))

    sample.df$Detail_Label=paste(sample.df[, main.class],
        ifelse(is.na(sample.df[, sub.class]),"",paste("_",sub.class, sample.df[,sub.class],sep="")),
        sep="")
    sample.df$Detail_Label=factor(sample.df$Detail_Label, levels=sort(unique(sample.df$Detail_Label)))

    sample.df$ReadsNum=get_read_number(sample.df[, read.col])
    if(sort){
      sample.df=sample.df[order(sample.df$Detail_Label, -sample.df$ReadsNum),]
    }

    return(sample.df)
}

get_read_number=function(mystr){
  mystr=gsub(",","",mystr)
  myvalue=sub("^(\\d+).*","\\1",mystr)
  return(as.numeric(myvalue))
}

get_read_percent=function(mystr){
  mystr=gsub(",","",mystr)
  mypercent=sub(".*\\((.+)\\%\\).*","\\1",mystr)
  return(as.numeric(mypercent)/100)
}

get_file_list=function(info.list, bw.norm.dir, bw.bias.dir, bias.suffix="_partition.bw"){
    data.list=list()
    for(n in 1:length(info.list)){
       
        dens.bw=paste(bw.norm.dir,"/",info.list[[n]]$Sample,"/",info.list[[n]]$Sample,".bw",sep="")
        names(dens.bw)=info.list[[n]]$Sample
       
        bias.bw=paste(bw.bias.dir,"/",info.list[[n]]$Sample, bias.suffix,sep="")
        names(bias.bw)=info.list[[n]]$Sample
        
        input.label=info.list[[n]]$Label
        names(input.label)=info.list[[n]]$Sample
        
        input.name=sub("^([^_]+)_(.+)$","\\1",info.list[[n]]$Sample)
        input.name=paste(input.name,"_",info.list[[n]]$PatientID,"_",round(info.list[[n]]$ReadsNum/1e6),"M",sep="")
        names(input.name)=info.list[[n]]$Sample
        
       
        if("Detail_Label" %in% colnames(info.list[[n]])){
          detail.label=info.list[[n]]$Detail_Label
          names(detail.label)=info.list[[n]]$Sample

          data.list[[n]]=list(Dens=dens.bw, Bias=bias.bw, Label=input.label, Detail_Label=detail.label, Name=input.name)
        }else{
          data.list[[n]]=list(Dens=dens.bw, Bias=bias.bw, Label=input.label, Name=input.name)
        }

    }
    names(data.list)=names(info.list)

    return(data.list)
}


sample_pred=function(my.sample, input.label=NULL, pred.level, 
  one.model.dir, out.dir=NULL,mid.dir=NULL,
  value.type="abs", transform=F, N_model=10, model.name="GLMNET_model", N_thread=12, force=F, rm.ref=T, hp.color=NA){
    result.list=list()
    for(j in 1:length(model.name)){
        cat("Predict by", model.name[j],"\n")
        if(!is.null(out.dir)){
          pred.result.dir=paste(out.dir,"/",model.name[j],"_test",sep="")
          if(!dir.exists(pred.result.dir)){dir.create(pred.result.dir, recursive = T)}
        }else{
          pred.result.dir=NULL
        }
        
        model.list=read_my_model(model.dir=one.model.dir, model.name=model.name[j], N_model=N_model)
        pred.list=list() 
        ##simplify
        for (i in 1:length(model.list)){
            cat("Predict by",model.name[j],":",i,"model\n")
            if(rm.ref==T){
                
                model.input.file=paste(one.model.dir,"/Model_Input/","model",i,"_input.txt",sep="")
                model.input=read.table(model.input.file, header = T, sep="\t", quote="\"",stringsAsFactors = F)
                model.sample=colnames(model.input)
                keep.ind=!names(my.sample) %in% model.sample
            }else{
                keep.ind=1:length(my.sample)
            }
            keep.label=input.label[keep.ind]

            rand.bw=my.sample[keep.ind]
            rand.bed=paste(one.model.dir,"/feature/","rand",i,".bed",sep="")
            if(!is.null(out.dir)){
              pred.mat.file=paste(out.dir,"/Prediction/","prediction",i,"_input.txt",sep="")
              if(!dir.exists(dirname(pred.mat.file))){dir.create(dirname(pred.mat.file), recursive = T)}
            }else{
              pred.mat.file=NULL
            }
            

            pred.mat=get_bw_mat(bw.file=rand.bw, region.file=rand.bed, 
                      mat.file=pred.mat.file, mid.dir=mid.dir,
                      N_threads=N_thread, force=force)


            if(value.type=="abs"){
              pred.mat=abs(pred.mat)
            }
            train.input.file=paste(one.model.dir,"/Model_Input/","model",i,"_input.txt",sep="")
            train.mat=read.table(train.input.file, header = T, sep="\t", quote="\"",stringsAsFactors = F)
            if(value.type=="abs"){
              train.mat=abs(train.mat)
            }
            train.avg=apply(train.mat, 1, mean, na.rm=T)
            if(transform==T){
              pred.final.mat=(pred.mat - train.avg) / apply(train.mat, 1, sd, na.rm=T)
              pred.final.mat[is.na(pred.final.mat)]=0
            }else{
              pred.final.mat=pred.mat
              pred.final.mat[is.na(pred.final.mat)]=0
            }

            if(!is.null(out.dir)){
              cat("plot signal heatmap...\n")
              if(is.null(input.label[1])|is.na(input.label[1])){
                col.anno=NA
              }else{
                col.anno=data.frame(Label=input.label)
                rownames(col.anno)=names(my.sample)
              }
              pred.zmat=(pred.mat - train.avg) / apply(train.mat, 1, sd, na.rm=T)
              pred.hp.file=sub("\\.txt$",".pdf", pred.mat.file)
              pred.zhp.file=sub("\\.txt$","_zmat.pdf", pred.mat.file)
              hplot1=draw_heatmap(list(pred_mat=pred.final.mat), graph.file=pred.hp.file, 
                              smooth=NA, my.color=hp.color, ncolor=100,
                              fontsize=10,
                              sort=F,show_rownum=30,show_colnum=1,width=7,height=7,
                              annotation_col=col.anno,
                              cluster_rows=F,cluster_cols=F)
              hplot2=draw_heatmap(list(zmat=pred.zmat), graph.file=pred.zhp.file, 
                              smooth=NA, my.color=hp.color, ncolor=100,
                              zMin = -2,zMax=2,
                              fontsize=10,
                              sort=F,show_rownum=30,show_colnum=1,width=7,height=7,
                              annotation_col=col.anno,
                              cluster_rows=F,cluster_cols=F)
            }

            cat("Predict by",i, model.name[j],"\n")
            one.pred=my_model_pred(list(model.list[[i]]), mod.name=paste(model.name[j],i,sep="_"), input.mat=t(pred.final.mat), pred.level=pred.level, input.label=keep.label)
            pred.list[[i]]=one.pred[[1]]
            names(pred.list)[i]=names(model.list)[i]

            if(!is.null(pred.result.dir)){
              model.result.file=paste(pred.result.dir,"/", model.name[j],"_",i,"_result.txt",sep="")
              write.table(pred.list[[i]], model.result.file,quote=F,sep="\t",row.names = F,col.names = T)
            }
        }

        out.avg=pred_average(pred.list = pred.list)
        if(!is.null(pred.result.dir)){
          avg.result.file=paste(pred.result.dir,"/", model.name[j],"_avg","_result.txt",sep="")
          write.table(out.avg, avg.result.file,quote=F,sep="\t",row.names = F,col.names = T)
        }
        result.list[[j]]=pred.list
        names(result.list)[j]=model.name[j]
    }
    return(result.list)
}



read_my_model=function(model.dir, model.name="GLMNET_model",N_model=1){
  
  model.type.dir=paste(model.dir,"/",model.name,sep="")
  if(!is.null(N_model)){
    model.file=paste(model.type.dir,"/",model.name,"_",1:N_model,".rds",sep="")
    model.list=lapply(model.file, readRDS)
    names(model.list)=paste(model.name,"_",1:N_model,sep="")
  }else{
      
      model.file=list.files(path=model.type.dir, pattern=paste("^",model.name,".*\\.rds$",sep=""),full.names = T)
      rand.id=sub(paste("^",model.name,"_(\\d+).*",sep=""), "\\1",basename(model.file))
      cancer.id=sapply(1:length(model.file), function(x){
        sub.match=regexpr(paste("(?<=",model.name,"_",rand.id[x],"_)","\\w+","(?=\\.rds$)",sep=""),
                          basename(model.file)[x],
                          perl=T)
        sub.model=ifelse(sub.match==-1,NA,regmatches(basename(model.file)[x], sub.match))
      })

      if(length(unique(cancer.id))<2){
        model.list=lapply(model.file, readRDS)
      }else{
        rand.type=unique(rand.id)
        cancer.type=unique(cancer.id)
        cat("Read",length(rand.type),"of",length(cancer.type), "submodel:",cancer.type,"\n")
        model.list=list()
        for(i in 1:length(rand.type)){
          sub.model.file=model.file[grepl(paste("^",model.name,"_",rand.type[i],sep=""), basename(model.file))]
          sub.model.list=lapply(sub.model.file, readRDS)
          names(sub.model.list)=sub("\\..*","",basename(sub.model.file))
          model.list[[i]]=sub.model.list
          names(model.list)[i]=paste(model.name,"_",rand.type[i],sep="")
        }
      }
    
  }
  
  return(model.list)
}

get_bw_mat=function(bw.file, region.file, mat.file=NULL, mid.dir, N_threads=12, force=F, verbose=F){

  if(!is.null(mat.file) & !force){
    if(file.exists(mat.file)){
      cal.flag=T
    }else{
      cal.flag=T
    }
  }else if(is.null(mat.file) & !force){
    cal.flag=T
  }else{
    cal.flag=T
  }

  if(cal.flag==F){

  }else{
    cat("Cal input...\n")
    if(!dir.exists(mid.dir)){dir.create(mid.dir, recursive = T)}
    score.files=paste(mid.dir,"/",names(bw.file),"_score.txt",sep="")
    names(score.files)=names(bw.file)
    
    stepi=round(length(bw.file)/10)
    pb=txtProgressBar(min=0, max=10)
    for(k in 1:length(bw.file)){
      if(verbose){
        cat(k, bw.file[k],"...\n")
      }else{
        setTxtProgressBar(pb, round(k/stepi))
      }

    }
    close(pb)

    
    score.mat=build_bias_mat(score.files)
    score.mat=score.mat[,-c(1,2)]
    my.gr=import.bed(region.file)
    rownames(score.mat)=my.gr$name
    if(!is.null(mat.file)){
      if(!dir.exists(dirname(mat.file))){dir.create(dirname(mat.file), recursive = T)}
      write.table(score.mat, file=mat.file, quote=F, sep="\t", row.names = T, col.names = T)
    }
  }
  return(score.mat)
}


build_bias_mat=function(peak.bias.file){
  bias.mat=c()
  for(i in 1:length(peak.bias.file)){
    window.df=as.data.frame(fread(peak.bias.file[i], header = T, sep="\t", quote = "\"", stringsAsFactors = F))
    if(i==1){
      bias.mat=data.frame(region=paste(window.df[,1],":",window.df[,2],"-",window.df[,3],sep=""),
                          name=window.df$name,
                          score=window.df$score)
    }else{
      bias.mat=cbind(bias.mat, window.df$score)
    }
    colnames(bias.mat)[i+2]=names(peak.bias.file)[i]
  }
  
  return(bias.mat)
}


draw_heatmap=function(mat.file, graph.file, 
                      zMin=NA, zMax=NA,zRange=c(0,1),smooth=NA,
                      my.color=NA, ncolor=100, 
                      sort=FALSE,show_rownum=1,show_colnum=1,width=NA,height=NA,append=F,
                      hrow.min=100,hrow.max=1000,combine=F,...){
  
  
  obj.list=c()
  sub.label=c()
  sub.height=c()
  for(i in 1:length(mat.file)){
    cat(i,"loading data of ")
    if(is.matrix(mat.file[[i]])){
      plot.mat=mat.file[[i]]
      plot.name=names(mat.file)[[i]]
    }else if(is.data.frame(mat.file[[i]])){
      plot.mat=mat.file[[i]]
      plot.name=names(mat.file)[[i]]
    }else if(is.vector(mat.file[[i]]) & length(mat.file[[i]])==1){
      if(file.exists(mat.file[[i]])){
        if(grepl("\\.xlsx$", mat.file[[i]])){
          plot.mat=read.xlsx(mat.file[[i]],colNames = T, rowNames = T)
          plot.name=sub("(\\.mat)*\\.xlsx$", "", basename(mat.file[[i]]))
        }else{
          plot.mat=read.table(mat.file[[i]], header = T, sep="\t", quote = "\"", stringsAsFactors = F,check.names = F)
          plot.name=sub("\\..*$", "", basename(mat.file[[i]]))
        }
      }else{
        stop(paste(i,mat.file[[i]],"not exists!\n"))
      }
    }else{
      stop(paste(i,names(mat.file)[i],"is not a matrix!\n"))
    }
    cat(plot.name,"\n")
    
    remove.ind=which(apply(plot.mat, 1, function(x){all(is.nan(x))}))
    if(length(remove.ind)>0){
      cat(paste(rownames(plot.mat)[remove.ind], collapse = ","), "removed because of NaN\n")
      plot.mat=plot.mat[-c(remove.ind), ]
    }
    
    if(is.na(zMin) | is.na(zMax)){
      vmin=max(min(unlist(plot.mat),na.rm=T), mean(unlist(plot.mat),na.rm=T)-2*sd(unlist(plot.mat),na.rm=T))
      vmax=min(max(unlist(plot.mat),na.rm=T), mean(unlist(plot.mat),na.rm=T)+2*sd(unlist(plot.mat),na.rm=T))
      my.range=c(vmin,vmax)
    }
    if(is.na(zMin)){zMin=my.range[1]}
    if(is.na(zMax)){zMax=my.range[2]}
   
    if(!is.na(smooth)){ plot.mat=apply(plot.mat, 1, function(x){rollapply(x, width=2*smooth, function(y){mean(y,na.rm=T)})})}

    if(sort==T){
      cat("sorting the matrix...\n")
      w1=apply(abs(plot.mat), 1, function(x){mean(x,na.rm=T)})
      w2=apply(abs(plot.mat), 1, function(x){sd(x,na.rm = T)})
      w2[w2==0]=NA
      w2[is.na(w2)]=min(w2,na.rm = T)
      w=w1

      plot.mat=plot.mat[order(w, decreasing = T),]
    }
    plot.mat[plot.mat > zMax]=zMax
    plot.mat[plot.mat < zMin]=zMin

    if(is.na(my.color[1])){
      my.color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(ncolor)
    }else{
      if(all(grepl("^#",my.color))){

      }else{
        my.color=colorRampPalette(my.color)(ncolor)
      }
    }
    
    if((show_rownum >=0)&(show_rownum <= 1)){
      row.num=nrow(plot.mat)*show_rownum
    }else{
      row.num=show_rownum
    }
    if(row.num >1){
      cut=ceiling(nrow(plot.mat) / row.num)
      show.ind=seq(cut,nrow(plot.mat),by=cut)
      tmp.name=rownames(plot.mat)[show.ind]
      plot.mat=as.matrix(plot.mat)
      rownames(plot.mat)=rep("",nrow(plot.mat))
      rownames(plot.mat)[show.ind]=tmp.name
    }else if(row.num==0){
      plot.mat=as.matrix(plot.mat)
      rownames(plot.mat)=rep("",nrow(plot.mat))
    }
    
    if((show_colnum >=0)&(show_colnum <= 1)){
      col.num=ncol(plot.mat)*show_colnum
    }else{
      col.num=show_colnum
    }
    if(col.num >1){
      cut=ceiling(ncol(plot.mat)/ col.num)
      show.ind=seq(cut,ncol(plot.mat),by=cut)
      tmp.name=colnames(plot.mat)[show.ind]
      plot.mat=as.matrix(plot.mat)
      colnames(plot.mat)=rep("",ncol(plot.mat))
      colnames(plot.mat)[show.ind]=tmp.name
    }else if(col.num==0){
      colnames(plot.mat)=rep("",ncol(plot.mat))
    }
    
    my.break=seq(zMin, zMax, length.out = ncolor+1)
    
    if(nrow(plot.mat) <= hrow.min){
      p.height=2
    }else if(nrow(plot.mat) >= hrow.max){
      p.height=height
    }else{
      p.height=2+round((height-2) * (nrow(plot.mat)-hrow.min)/(hrow.max-hrow.min), 1)
    }
    sub.height[i]=p.height
    
    sub.label[i]=paste(nrow(plot.mat),"rows",sep=" ")
    p.obj=pheatmap(plot.mat,
                   color=my.color,
                   breaks=my.break,
                   border_color = "NA",
                   main =paste(plot.name,"\n",nrow(plot.mat),"rows", sep=" "),
                   width=width,height=p.height
                   ,...
                   )
    obj.list[[i]]=p.obj
    names(obj.list)[i]=plot.name
  }
  
  gg.list=lapply(obj.list, as.ggplot)
  names(gg.list)=names(obj.list)
  if(combine){
    

    merge.plot=plot_grid(plotlist=gg.list, align = "v",ncol=1,
                         rel_heights = sub.height
                        
    )
  }

  if(grepl("\\.pdf",basename(graph.file))){
    if(!dir.exists(dirname(graph.file))){dir.create(dirname(graph.file),recursive = T)}
    if(!combine){
      pdf(graph.file,width=width,height=height,onefile = T,bg="white")
      lapply(gg.list,function(x){
        
        grid.draw(x)
      })
      dev.off()
    }else{
      ggsave(graph.file,plot=merge.plot,width=width,height = height)
    }
   
  }else{
    if(!dir.exists(graph.file)){dir.create(graph.file,recursive = T)}
    if(!combine){
      for(i in 1:length(gg.list)){
        jpeg(filename =file.path(graph.file,paste(names(gg.list)[i],".jpg",sep="")),width=width,height=height,units="in",res=300)
       
        grid.draw(gg.list[[i]])
        
        dev.off()
      }
    }else{
      jpeg(filename=paste(graph.file,".jpg",sep=""),width=width,height=height,units="in",res=300)
      grid.newpage()
      
      grid.draw(merge.plot)
      
      dev.off()
    }
    
  }
  
  if(!combine){
    result=list(gg.list)
  }else{
    result=merge.plot
  }
  
  return(result)
}

my_model_pred=function(mod.list, mod.name, input.mat, 
                       pred.level=NULL, 
                       input.label=NULL){
  if(length(mod.name) != length(mod.list)){stop("model list should be same length as model names!\n")}
  
  pred.list=list()
  for(i in 1:length(mod.list)){
    
      my.pred=predict(mod.list[[i]], input.mat,type="prob")
      my.pred$prediction=predict(mod.list[[i]],input.mat)
    
    if(!is.null(input.label)){
      my.pred$Label=input.label[rownames(my.pred)]
    }
    

    my.pred=cbind(data.frame(Sample=rownames(input.mat)),my.pred)
    
    pred.list[[i]]=my.pred
  }
  names(pred.list)=mod.name
  
  return(pred.list)
}

pred_average=function(pred.list){
  all.pred.df=Reduce(rbind, pred.list)
  label.level=colnames(all.pred.df)[!colnames(all.pred.df) %in% c("Sample","prediction","Label")]

  pred.unique.list=list()
  
  for(x in label.level){
    avg.value=tapply(all.pred.df[,x], all.pred.df[,1], mean)
    avg.df=data.frame(names(avg.value), avg.value)
    colnames(avg.df)=c("Sample", x)
    pred.unique.list[[x]]=avg.df
  }
  pred.unique.df=Reduce(merge, pred.unique.list)
 
  pred.ind=apply(pred.unique.df[,-1,drop=F], 1, which.max)
  pred.unique.df[,"prediction"]=label.level[pred.ind]
  
  if(any(grepl("^Label$",colnames(all.pred.df)))){
    sample.label=unique(all.pred.df[,c("Sample","Label")])
    pred.unique.df=merge(pred.unique.df, sample.label)
  }
  

  return(pred.unique.df)
}


evaluate_pred=function(result.list, input.label=NA, input.name=NA, normal.cut=NA, out.dir, my.color=c("red","pink","orange","blue"), verbose=F){
    if(is.na(input.name)[1]){
      input.name=rownames(result.list[[1]][[1]])
    }
    if(!is.na(input.label[1])){
      info.df=data.frame(Sample=names(input.label), Label=input.label, Name=input.name)
    }else{
      info.df=data.frame(Sample=result.list[[1]][[1]]$Sample, Label=input.label, Name=input.name)
    }
    prob.mat=c()
    pred.col=colnames(result.list[[1]][[1]]) 
    pred.level=pred.col[!pred.col %in% c("Sample","prediction","Label") ]
    for(n in 1:length(result.list)){
      one.prob.mat=problist_to_mat(df.list=result.list[[n]], pred.level=pred.level,prefix=names(result.list)[n]) 
      prob.mat=rbind(prob.mat, one.prob.mat)
    }
    
    mm.name=paste(sub("_.*","", names(result.list)),collapse = "_")
    prob.hp.file=paste(out.dir,"/",mm.name,"_","result_prob_mat.pdf",sep="")
    if(!is.na(input.label[1])){
      col.anno=data.frame(Label=input.label)
      rownames(col.anno)=names(input.label)
    }else{
      col.anno=NA
    }
    
    hp.plot=draw_heatmap(list(prob_mat=prob.mat), graph.file=prob.hp.file, 
                            smooth=NA, my.color=NA, ncolor=100,
                            fontsize=10,
                            
                            sort=F,show_rownum=30,show_colnum=1,width=7,height=7,
                            
                            annotation_col=col.anno,
                            
                            cluster_rows=F,cluster_cols=F)

    for(j in 1:length(result.list)){
        pred.result.dir=paste(out.dir,"/", names(result.list)[j],"_test",sep="")
        if(!dir.exists(pred.result.dir)){dir.create(pred.result.dir, recursive = T)}
        pred.list=result.list[[j]]
        for(i in 1:length(pred.list)){
            pred.list[[i]]=merge(pred.list[[i]], info.df[,c("Sample","Label")], sort=F)
            my.pred=pred.list[[i]]
            model.graph.file=paste(pred.result.dir,"/", names(result.list)[j],"_",i,"_prediction.pptx",sep="")
            model.result.file=paste(pred.result.dir,"/", names(result.list)[j],"_",i,"_result.txt",sep="")
           
            if(verbose){
            save_prediction_results(my.pred=my.pred,
                                    normal.cut=normal.cut,
                                    info.df=info.df,
                                    title=paste(names(result.list)[j],nrow(my.pred),"Sample"), font.size = 10,
                                    model.graph.file=model.graph.file, model.result.file=model.result.file, model.confusion.file=model.confusion.file
                                    ,my.color=my.color
                                    )
            }
        }
        cat("AUROC of", names(result.list)[j],":",length(pred.list),"models\n")
        avg.graph.file=paste(pred.result.dir,"/", names(result.list)[j],"_avg","_prediction.pdf",sep="")
        avg.result.file=paste(pred.result.dir,"/", names(result.list)[j],"_avg","_result.txt",sep="")
        
        if(!is.na(input.label[1])){
          auc.list=lapply(pred.list, pred_AUC)
          all.auc=Reduce(c, auc.list)
          all.auc.df=data.frame(Class=names(all.auc), AUROC=all.auc, color="blue")
          box.plot=draw_boxplot_style(all.auc.df,yMin=0,yMax=1,y.breaks=0.2,
                                      title=paste(model.name[j],"test sample"),
                                      xlab="Class", ylab="AUROC",
                                      font.size=12,legend.ncol=3,
                                      mycolor="blue",dot=T,text.x.angle=90)
          box.graph=paste(pred.result.dir,"/","pred_boxplot.pdf",sep="")
          pdf(box.graph, width=3,height=5,onefile = T,bg="white")
          lapply(list(box.plot),function(x){
            grid.draw(x)
          })
          dev.off()
        }
       
        cat("Save average prediction of",names(result.list)[j],"...\n")
        pred.avg=pred_average(pred.list = pred.list)
        save_prediction_results(my.pred=pred.avg,
                                normal.cut=normal.cut,
                                info.df=info.df,
                                title=paste(names(result.list)[j],nrow(pred.avg),"Sample"), font.size = 10,
                                model.graph.file=avg.graph.file, 
                                model.result.file=avg.result.file, 
                               
                                append=T,my.color=my.color
                                )
    }
}

problist_to_mat=function(df.list, pred.level, prefix="Dens"){
    all.pred=c()
    all.long=c()
    for(m in 1:length(df.list)){
      one.df=df.list[[m]]
      keep.id=grep("Sample|Label",colnames(one.df), value=T)
      one.long=melt(one.df, id.vars=keep.id, measure.vars=pred.level, variable.name="Cancer_type", value.name = "Prob")
      one.long$source=paste(prefix,"_",one.long$Cancer_type,"_",m,sep="")
      all.long=rbind(all.long, one.long)
      if("Label" %in% colnames(all.long)){
        
      }else{
        
      }
      
    }
    all.pred=dcast(all.long, source~Sample, value.var="Prob")
    rownames(all.pred)=all.pred$source
    all.pred=all.pred[match(unique(all.long$source), rownames(all.pred)), unique(all.long$Sample)]

    return(all.pred)
}

combine_problist=function(result1.list, result2.list, prefix1="Dens", prefix2="Bias",train.level){
  prob.mat1=c()
  for(n in 1:length(result1.list)){
    one.prob.mat=problist_to_mat(result1.list[[n]], train.level, prefix=paste(prefix1,"_",names(result1.list)[n],sep=""))
    prob.mat1=rbind(prob.mat1, one.prob.mat)
  }
  prob.mat2=c()
  for(n in 1:length(result2.list)){
    one.prob.mat=problist_to_mat(result2.list[[n]], train.level, prefix=paste(prefix2,"_",names(result2.list)[n],sep=""))
    prob.mat2=rbind(prob.mat2, one.prob.mat)
  }
  all.prob.mat=rbind(prob.mat1, prob.mat2)

  return(all.prob.mat)
}


save_prediction_results=function(my.pred, test.pred=my.pred, normal.cut=NA,
                                 info.df, gr="group", sort.by="ReadsNum",
                                 title
                                 ,model.graph.file, model.result.file
                                 , model.confusion.file=NA
                                 ,font.size=12,legend.ncol=3, my.color=c("red","pink","orange","blue")
                                 ,xlab="plasma samples", ylab="Probability"
                                 ,append=F){
  if(!dir.exists(dirname(model.graph.file))){dir.create(dirname(model.graph.file), recursive = T)}
  if(!dir.exists(dirname(model.result.file))){dir.create(dirname(model.result.file), recursive = T)}

  my.pred=merge(info.df[,c("Sample","Label")],my.pred, sort=F)
  test.pred=my.pred
  write.table(my.pred, model.result.file,quote=F,sep="\t",row.names = F,col.names = T)
 
  if(!is.na(normal.cut)){
    cut.pred=my.pred
    cut.pred$prediction=NA
    normal.label=grep("normal",colnames(cut.pred),value=T)
    cut.pred[cut.pred[,normal.label]>normal.cut,"prediction"]=normal.label
    cancer.pred=cut.pred[cut.pred[,normal.label] <= normal.cut, ]
    cancer.mat=cancer.pred[,!colnames(cancer.pred) %in% c("Sample","prediction","Label", normal.label)]
    cancer.pred$prediction=colnames(cancer.mat)[apply(cancer.mat,1,which.max)]
    cut.pred=df_update(old.df=cut.pred, map.df=cancer.pred, key.col="Sample", info.col="prediction")
    write.table(cut.pred, 
      sub("\\.txt",paste("_normal.cut",normal.cut,".txt",sep=""),model.result.file),quote=F,sep="\t",row.names = F,col.names = T)
  }
  
  bar.list=draw_prediction_bar(my.pred=my.pred, 
                      info.df=info.df, 
                      gr=gr, sort.by=sort.by
                      ,font.size=font.size,legend.ncol=legend.ncol,barcolor=my.color,
                      xlab=xlab, ylab=ylab
                      )
  pdf(model.graph.file,onefile = T,bg="white")
  lapply(bar.list,function(x){
    grid.draw(x)
  })
  dev.off()

  if(nrow(test.pred)>0 & length(unique(test.pred$Label))>1){

    if(!is.na(normal.cut)){

      normal.pred=cut.pred
      normal.pred$Label=as.character(normal.pred$Label)
      normal.pred$Label[normal.pred$Label!=normal.label]="Cancer"
      normal.pred$prediction[normal.pred$prediction!=normal.label]="Cancer"
      normal.tb=as.data.frame.matrix(table(real=normal.pred$Label, pred=normal.pred$prediction))
      colnames(normal.tb)=paste("pred_",colnames(normal.tb),sep="")
      rownames(normal.tb)=paste("real_",rownames(normal.tb),sep="")
      
      cancer.pred=cut.pred
      cancer.pred=cancer.pred[cancer.pred$Label!=normal.label,]
      cancer.tb=as.data.frame.matrix(table(cancer.pred$Label, cancer.pred$prediction))
      colnames(cancer.tb)=paste("pred_",colnames(cancer.tb),sep="")
      rownames(cancer.tb)=paste("real_",rownames(cancer.tb),sep="")
    }


    pROC.plot=draw_pROC(my.pred=test.pred, model.title=title, my.color = my.color)
    auc.graph.file=sub("\\.pdf$", "_AUC.pdf",model.graph.file)
    pdf(auc.graph.file, onefile = T,bg="white")
    lapply(list(pROC.plot),function(x){
      grid.draw(x)
    })
    dev.off()

    
    my.pred.long=melt(test.pred, id.vars=c("Sample","Label","prediction"),
                      variable.name = "Prediction_type",
                      value.name = "Probability")
    my.pred.long$Label=paste("real_",my.pred.long$Label, sep="")
    my.pred.long$Prediction_type=paste(my.pred.long$Prediction_type, "_prob",sep="")
    my.pred.long$color=my.pred.long$Prediction_type
    bar.group.plot=draw_color_bar_style(my.pred.long[,c("Label","Probability","color","Prediction_type")],
                                        group="Prediction_type",
                                  interval=T, yMin=0,yMax=1,y.dashed=NULL,
                                  xlab="", ylab="Average Probability",
                                  title=paste(nrow(test.pred), "Sample Average"),
                                  font.size=font.size, legend.ncol=legend.ncol,
                                  barcolor=my.color, flip=F,hide.x=F,
                                  text.x.angle=90, x.fac=T)+
      facet_grid(Prediction_type~Label, scales = "free",space="free")+
      theme(strip.text = element_text(size = 10, face="bold"))
    avg.bar.file=paste(dirname(model.graph.file),"/avg_bar.pdf",sep="")
      #sub("\\.pdf$", "_avg_bar.pdf",model.graph.file)
    pdf(avg.bar.file, width=5, height=7, onefile = T,bg="white")
    lapply(list(bar.group.plot),function(x){
      grid.draw(x)
    })
    dev.off()

    test.group.file=paste(dirname(model.result.file),"/avg_group.txt",sep="")
      #sub("\\.txt$","_vali_group.txt", model.result.file)
    test.pred.list=split(test.pred[, !colnames(test.pred) %in% c("Sample","Label","prediction")], test.pred$Label)
    test.pred.group=c()
    for(i in 1:length(test.pred.list)){
      avg.group=apply(test.pred.list[[i]], 2, mean)
      test.pred.group=rbind(test.pred.group, avg.group)
    }
    colnames(test.pred.group)=paste(colnames(test.pred.list[[1]]),"_prob",sep="")
    rownames(test.pred.group)=paste("real_",names(test.pred.list), sep="")
    write.table(as.data.frame(t(test.pred.group)), test.group.file, quote=F,sep="\t",row.names = T,col.names = T)
    
  }
  
}

pred_AUC=function(my.pred){
  
  if(nrow(my.pred)<1){return(c())}
  if(!"Label" %in% colnames(my.pred)){stop("Label is not found!")}
  
  my.pred$Label=make.names(my.pred$Label)
  my.label=unique(make.names(my.pred$Label))
  pred.auc=c()
  if(length(unique(my.pred$Label))>1){
    for(n in 1:length(my.label)){
      my.respone=ifelse(my.pred$Label==my.label[n],1,0)
      one.roc=roc(my.respone,my.pred[,my.label[n]], direction = "<", quiet=T)
      pred.auc[n]=round(as.numeric(one.roc$auc),2)
      names(pred.auc)[n]=my.label[n]
    }
  }
  return(pred.auc)
}


draw_boxplot_style=function(data.df,yMin=NA,yMax=NA,y.breaks=NA,
                            y.dashed=NULL,xlab="x", ylab="y",title="boxplot",subtitle="",
                            font.size=12,legend.ncol=NA,mycolor=NA,dot=T,text.x.angle=0, flip=F){
  

  if(is.na(mycolor[1])){mycolor=c("black", "#0073C2FF","#EFC000FF","red","#009E73","orange","purple", "grey")}
  if(is.na(yMin)){
    ymin=range(as.numeric(data.df[,2]),na.rm=T)[1]
  }else{
    ymin=yMin
  }
  if(is.na(yMax)){
    ymax=range(as.numeric(data.df[,2]),na.rm=T)[2]
  }else{
    ymax=yMax
  }
  data.df[,3]=factor(data.df[,3], levels =sort(unique(data.df[,3])), ordered = T)
  
  box.plot=ggplot(data.df,aes(x = data.df[,1],
                              y=as.numeric(data.df[,2])
                              ))+
    geom_boxplot(aes(color=data.df[,3]),
                 position = position_dodge(0.8),outlier.shape=NA)
  if(dot){
    box.plot=box.plot+geom_dotplot(aes(fill=data.df[,3]),position = position_dodge(0.8)
                                   ,binaxis = "y"
                                   ,dotsize=0.4
                                   ,stackdir = "center"
                                   )
  }
  box.plot=box.plot+scale_color_manual(name="group",breaks=levels(data.df[,3]),values=mycolor)+
    scale_fill_manual(name="group", breaks=levels(data.df[,3]), values=mycolor)+
    theme_bw()+
    theme(
      axis.text=element_text(size=font.size, face="bold", colour = "black"),
      axis.title = element_text(size=font.size+4,face="bold",colour="black"),
      axis.title.x=element_text(size=font.size,margin=ggplot2::margin(t=0,r=0,b=0,l=0)),
      axis.title.y=element_text(size=font.size),
      axis.ticks = element_line(size = 1, colour = "black"),
      axis.text.y=element_text(hjust=1,margin=ggplot2::margin(t=0,r=5,b=0,l=0)),
      axis.text.x = element_text(angle=text.x.angle,vjust = 0.5,hjust=0.8),
      panel.grid=element_blank(), panel.background = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = font.size,face="bold"),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )+
    geom_hline(yintercept = y.dashed,linetype="dashed",color="black", size=1)+
    labs(title=title,x=xlab, y=ylab,subtitle = subtitle)
  
  if(!is.na(yMin)|!is.na(yMax)){
    box.plot=box.plot+scale_y_continuous(limits = c(ymin,ymax))
  }
  if(!is.na(legend.ncol)){
    box.plot=box.plot+guides(color=guide_legend(ncol=legend.ncol))
  }
  if(!is.na(y.breaks)){
    box.plot=box.plot+scale_y_continuous(limits = c(ymin,ymax),breaks = seq(ymin*10e4,ymax*10e4,y.breaks*10e4)/10e4)
  }
  if(flip){
    box.plot=box.plot+coord_flip()
  }
  
 return(box.plot)
}

df_update=function(old.df, map.df, key.col, info.col){
  if(any(!key.col %in% colnames(old.df))){
    stop("key not found in old dataframe:",key.col,"!\n")
  }
  if(any(!key.col %in% colnames(map.df))){
    stop("key not found in map dataframe:",key.col,"!\n")
  }
  new.df=old.df
  for(i in 1:length(info.col)){
    if(info.col[i] %in% colnames(old.df)){
      cat("Column",info.col[i],"will be updated\n")
    }else{
      new.df[, info.col[i]]=NA
      cat("Add column:",info.col[i],"\n")
    }
  }

  up.num=0
  up.df=c()
  for(i in 1:nrow(old.df)){
    ##each row
    match.id=match(old.df[i, key.col], map.df[,key.col])
    if(!is.na(match.id)){
      ##each info column
      for(j in 1:length(info.col)){
        new.df[i, info.col[j]]=paste(map.df[match.id, info.col[j]], collapse = ";")
      }
      up.df=rbind(up.df, new.df[i,])
      up.num=up.num+1
    }
  }
  cat("Update",up.num,"rows\n")

  return(new.df)
}

draw_prediction_bar=function(my.pred, info.df, gr="group",sort.by="ReadsNum",
                             font.size=12,legend.ncol=3,barcolor=c("red","pink","orange","blue"),
                             xlab="plasma samples", ylab="Probability"){
  if("PatientID" %in% colnames(info.df) ){
    info.df$Name=sub("^([^_]+)_(.+)$","\\1",info.df$Sample)
    info.df$Name=paste(info.df$Name, info.df[,"PatientID"], sep="_")
  }
  if("Condition" %in% colnames(info.df) ){
    info.df$Name=paste(info.df$Name, sub("(raw|clean).*","\\1", info.df[,"Condition"]), sep="_")
  }
  if("Input_ssDNA(ng)" %in% colnames(info.df) ){
    info.df$Name=paste(info.df$Name, paste(info.df[,"Input_ssDNA(ng)"],"ng",sep=""), sep="_")
  }
  if("SS(ng/ul)" %in% colnames(info.df) ){
    info.df$Name=paste(info.df$Name, paste(info.df[,"SS(ng/ul)"],"/ul",sep=""), sep="_")
  }
  if("Rate" %in% colnames(info.df) ){
    info.df$Name=paste(info.df$Name, paste(round(info.df[,"Rate"]*100),"%",sep=""), sep="_")
  }
  if("ReadsNum" %in% colnames(info.df) ){
    info.df$Name=paste(info.df$Name, paste(round(info.df[,"ReadsNum"]/1e6),"M",sep=""), sep="_")
  }
  if("library_type" %in% colnames(info.df)){
    info.df$Name=paste(info.df$Name, info.df[,"library_type"], sep="_")
  }
  
  
  my.pred.info=merge(info.df,my.pred, sort=F)

  if(all(sort.by %in% colnames(my.pred.info)) ){
    sort.ind=do.call(order, my.pred.info[,c("Label",sort.by)])
    my.pred.info=my.pred.info[sort.ind, ]

  }else{

  }
  
  my.pred.long=melt(my.pred.info, id.vars=unique(c(colnames(info.df), "Label","prediction","Name")),
                    variable.name = "Prediction_type",
                    value.name = "Probability")
  
  
  label.type=sort(unique(my.pred.long$Label))
  if(length(label.type)==0){
    label.type="unknown"
    my.pred.long$Label="unknown"
  }
  bar.list=lapply(1:length(label.type), function(n){
    sub.pred=my.pred.long[my.pred.long$Label==label.type[n], ]
    bar.plot=draw_bar_percent(sub.pred[,c("Name","Probability","Prediction_type")],y.dashed=NULL,
                              xlab=xlab, ylab=ylab,
                              title=paste(as.character(label.type[n]), length(unique(sub.pred$Sample)), "Sample"),
                              font.size=font.size,legend.ncol=legend.ncol,
                              barcolor=barcolor,
                              flip=F,text.x.angle=90)+
      theme(strip.text = element_text(size = 12, face="bold"))
  })

  return(bar.list)
}


draw_bar_percent=function(data.df,y.dashed=NULL,yMin=NA,yMax=NA,
                          xlab="x",ylab="y",title="percentage",subtitle="",y.breaks=NA,
                          font.size=12,legend.ncol=NA,barcolor=NA,flip=T,text.x.angle=90){
  
  if(is.na(barcolor[1])){barcolor=c("black", "#0073C2FF","#EFC000FF","red","#009E73","yellow","purple", "grey")}
  ymin=yMin
  ymax=yMax
  data.df[,3]=factor(data.df[,3], levels =sort(unique(data.df[,3])), ordered = T)
  
  barplot=ggplot(data.df,aes(x = factor(data.df[,1], levels =unique(data.df[,1]), ordered = T),
                             y=as.numeric(data.df[,2]),fill=data.df[,3]))+
    geom_bar(stat="identity",position = "fill")+
    scale_fill_manual(name="group", breaks=levels(data.df[,3]), values=barcolor)+
    theme_bw()+
    theme(
      axis.text=element_text(size=font.size, face="bold", colour = "black"),
      axis.title = element_text(size=font.size+4,face="bold",colour="black"),
      axis.title.x=element_text(size=font.size,margin=ggplot2::margin(t=0,r=0,b=0,l=0)),
      axis.title.y=element_text(size=font.size),
      axis.ticks = element_line(size = 1, colour = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(angle=text.x.angle, vjust = 1,hjust=1),
      axis.text.y=element_text(hjust=1,margin=ggplot2::margin(t=0,r=5,b=0,l=0)),
      panel.grid=element_blank(), panel.background = element_blank(), #remove grid
      axis.line.x = element_blank(),
      axis.line.y=element_line(colour="grey"),
      legend.position = "top",
      legend.text = element_text(size = font.size,face="bold"),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )+
    scale_y_continuous(limits = c(ymin,ymax),expand=expansion(mult=c(0,0.1)))+
    labs(x=xlab, y=ylab,title=title,subtitle=subtitle)
  if(is.numeric(y.dashed)){
    barplot=barplot+geom_hline(yintercept = y.dashed,linetype="dashed",color="black", size=1)
  }  
  if(!is.na(legend.ncol)){
    barplot=barplot+guides(fill=guide_legend(ncol=legend.ncol))
  }
  if(!is.na(y.breaks)){
    barplot=barplot+scale_y_continuous(limits = c(ymin,ymax),breaks = seq(ymin*10e4,ymax*10e4,y.breaks*10e4)/10e4)
  }
  if(flip){
    barplot=barplot+coord_flip()
  }
    
  return(barplot)  
}


draw_pROC=function(my.pred, model.title="My model",my.color=c("red","pink","orange","blue")){
  
  if(nrow(my.pred)<1){return()}
  
  my.pred$Label=make.names(my.pred$Label)
  my.label=unique(make.names(my.pred$Label))
  roc.df=c()
  cut.df=c()
  for(n in 1:length(my.label)){
    my.respone=ifelse(my.pred$Label==my.label[n],1,0)
    one.roc=roc(my.respone,my.pred[,my.label[n]], direction = "<", quiet=T)
    cut.roc=pROC::coords(one.roc, x="best", ret=c("threshold","specificity", "sensitivity"),best.weights=c(1,0.5))
    ci.roc=pROC::ci(one.roc)
    roc.df=rbind(roc.df, 
                 data.frame(Sensitivity=one.roc$sensitivities,
                            Specificity=one.roc$specificities,
                            group=paste(my.label[n],
                                        "(AUC=",round(as.numeric(one.roc$auc),2),
                                        ";95%CI=",round(ci.roc[1],2),"-",round(ci.roc[3],2),
                                        ")",
                                        sep="")
                            )
    )
    cut.df=rbind(cut.df,
                 data.frame(Threshold=cut.roc$threshold,
                            Sensitivity=cut.roc$sensitivity,
                            Specificity=cut.roc$specificity,
                            group=my.label[n]))
  }
  roc.df$Specificity=1-roc.df$Specificity
  cut.df$Specificity=1-cut.df$Specificity
  roc.plot=draw_line_style(data.df=roc.df[,c(2,1,3)],xlab="1-Specificity", ylab="Sensitivity",
                           title=paste(model.title,"pROC curve"),
                           font.size=12,legend.ncol=1,mycolor=my.color, dot.show=F)
  roc.plot=roc.plot+
    coord_cartesian(clip = "off")+
    geom_label_repel(data=cut.df,
                     aes(x=Specificity,
                         y=Sensitivity,
                         group=group,
                         label=paste(group,":",round(Threshold,2),
                          ";Spec=",round(1-Specificity, 3),";Sensi=",round(Sensitivity,3),
                          sep="")),
                         point.padding=1,ylim = c(-Inf, Inf))+
    geom_point(data=cut.df, aes(x=Specificity,y=Sensitivity,group=group),color="red",size=3)
  return(roc.plot)
}

draw_line_style=function(data.df,yMin=NA,yMax=NA,xMin=NA,xMax=NA,y.breaks=NA,x.breaks=NA,
                         y.dashed=NULL,x.dashed=NULL,xlab="x", ylab="y",title="line-plot",subtitle="",
                         font.size=12,legend.ncol=NA,mycolor=NA,linetype=NA,linesize=1,dot.show=T,dot.size=1,line="path"){
  
  if(is.na(mycolor[1])){mycolor=c("black", "#0073C2FF","#EFC000FF","red","#009E73","orange","purple", "grey")}
  if(is.na(linetype[1])){linetype=rep(c("solid","dashed"),4)}
  if(is.na(yMin)){ymin=range(as.numeric(data.df[,2]),na.rm=T)[1]}else{ymin=yMin}
  if(is.na(yMax)){ymax=range(as.numeric(data.df[,2]),na.rm=T)[2]}else{ymax=yMax}
  if(is.na(xMin)){xmin=range(as.numeric(data.df[,1]),na.rm=T)[1]}else{xmin=xMin}
  if(is.na(xMax)){xmax=range(as.numeric(data.df[,1]),na.rm=T)[2]}else{xmax=xMax}
  data.df[,3]=factor(data.df[,3], levels =sort(unique(data.df[,3])), ordered = T)
  
  line.plot=ggplot(data.df,aes(x =as.numeric(data.df[,1]),
                             y=as.numeric(data.df[,2]),group=data.df[,3]))
  if(grepl("path",line,ignore.case = T)){
    line.plot=line.plot+geom_path(aes(color=data.df[,3],linetype=data.df[,3]),size=linesize)
  }else if(grepl("smooth",line,ignore.case = T)){
    line.plot=line.plot+geom_smooth(method="loess",aes(color=data.df[,3],linetype=data.df[,3]),size=linesize)
  }
    
  line.plot=line.plot+scale_color_manual(name="group", breaks=levels(data.df[,3]), 
                       values=mycolor)+
    scale_linetype_manual(name="group",breaks=levels(data.df[,3]),values = linetype)+
    theme_bw()+
    theme(
      axis.text=element_text(size=font.size, face="bold", colour = "black"),
      axis.title = element_text(size=font.size+4,face="bold",colour="black"),
      axis.title.x=element_text(size=font.size,margin=ggplot2::margin(t=5,r=0,b=0,l=0)),
      axis.title.y=element_text(size=font.size),
      axis.ticks = element_line(size = 1, colour = "black"),
      axis.text.y=element_text(hjust=1,margin=ggplot2::margin(t=0,r=5,b=0,l=0)),
      panel.grid=element_blank(), panel.border=element_blank(), panel.background = element_blank(), 
      axis.line.y=element_line(colour="black", size=1),
      axis.line.x=element_line(colour="black", size=1),
      legend.position = "top",
      legend.text = element_text(size = font.size,face="bold"),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )+
    geom_vline(xintercept = x.dashed,linetype="dashed",color="black",size=1)+
    geom_hline(yintercept = y.dashed,linetype="dashed",color="black", size=1)
  if(!is.na(yMin) | !is.na(yMax)){
    line.plot=line.plot+scale_y_continuous(limits = c(ymin,ymax))
  }
  if(!is.na(xMin) | !is.na(xMax)){
    line.plot=line.plot+scale_x_continuous(limits=c(xmin,xmax))
  }
   
  line.plot=line.plot+labs(title=title,x=xlab, y=ylab,subtitle = subtitle)
  if(dot.show){
    line.plot=line.plot+geom_point(aes(color=data.df[,3]),alpha=0.3,size=dot.size)
  }
  if(!is.na(legend.ncol)){
    line.plot=line.plot+guides(color=guide_legend(ncol=legend.ncol))
  }
  if(!is.na(y.breaks)){
    line.plot=line.plot+scale_y_continuous(limits = c(ymin,ymax),breaks = seq(ymin*10e4,ymax*10e4,y.breaks*10e4)/10e4)
  }
  if(!is.na(x.breaks)){
    line.plot=line.plot+scale_x_continuous(limits = c(xmin,xmax),breaks = seq(xmin*10e4,xmax*10e4,x.breaks*10e4)/10e4)
  }
  
  return(line.plot)
}

draw_color_bar_style=function(data.df, group=NULL,interval=F, yMin=NA,yMax=NA,y.dashed=NULL,
                              xlab="x", ylab="y",title="barplot",subtitle="",y.breaks=NA,
                              font.size=12,legend.ncol=NA,barcolor=NA,flip=T,hide.x=F,
                              text.x.angle=0, x.fac=T){
  
  if(is.na(barcolor[1])){barcolor=c("black", "#0073C2FF","#EFC000FF","red","#009E73","orange","purple", "grey")}
  ymin=yMin
  ymax=yMax
  if(!is.null(group)){
    sum.df=summarySE(data.df, measurevar = colnames(data.df)[2], groupvars = c(colnames(data.df)[1],group), type="mean")
    sum.df=merge(data.df[,-2], sum.df, sort=F, by=c(colnames(data.df)[1],group))
    data.df=sum.df[, c(colnames(data.df)[1:3],"se",colnames(data.df[-c(1:3)]))]
    data.df=unique(data.df)
  }
  if(x.fac){
    data.df[,1]=factor(data.df[,1], levels =unique(data.df[,1]), ordered = T)
  }
  
  bar.plot=ggplot(data.df,aes(x = data.df[,1],
                              y=as.numeric(data.df[,2])))
    
  if(is.numeric(data.df[,3])){
    bar.plot=bar.plot+
      geom_bar(aes(fill=data.df[,3]),stat = 'identity',position = position_dodge(width=0.9),alpha=0.8)+
      scale_fill_gradient(name="group",low=barcolor[1],high=barcolor[2])
  }else{
    if(!is.factor(data.df[,3])){
      data.df[,3]=factor(data.df[,3], levels =sort(unique(data.df[,3])), ordered = T)
    }
    bar.plot=bar.plot+
      geom_bar(aes(fill=data.df[,3]),stat = 'identity',position = position_dodge(width=0.9),alpha=0.8)+
      scale_fill_manual(name="group", breaks=levels(data.df[,3]), values=barcolor)
  }
  if(interval==T){
    bar.plot=bar.plot+geom_errorbar(aes(x=data.df[,1],
                                        ymin=data.df[,2], ymax=data.df[,2]+data.df[,"se"],
                                        group=data.df[,3]), 
                                    width=.2, position=position_dodge(.9))
  }
  bar.plot=bar.plot+theme_bw()+
    theme(
      axis.text=element_text(size=font.size, face="bold", colour = "black"),
      axis.title = element_text(size=font.size+4,face="bold",colour="black"),
      axis.title.x=element_text(size=font.size,margin=ggplot2::margin(t=5,r=0,b=0,l=0)),
      axis.title.y=element_text(size=font.size,margin =ggplot2::margin(t=0,r=5,b=0,l=0)),
      axis.ticks = element_line(size = 1, colour = "black"),
      axis.text.x = element_text(angle=text.x.angle,vjust = 0.8,hjust=0.5),
      axis.text.y=element_text(hjust=1,vjust=0.5,margin=ggplot2::margin(t=0,r=5,b=0,l=0)),
      panel.grid=element_blank(), panel.background = element_blank(),
      axis.line.x=element_line(colour="black", size=1),
      axis.line.y=element_line(colour="black", size=1),
      legend.position = "top",
      legend.text = element_text(size = font.size,face="bold"),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )+
    geom_hline(yintercept = y.dashed,linetype="dashed",color="black", size=1)+
    scale_y_continuous(limits = c(ymin,ymax),expand=expansion(mult=c(0,0.1)))+
    labs(title=title,x=xlab, y=ylab,subtitle = subtitle)
  if(hide.x){
    bar.plot=bar.plot+theme(axis.ticks.x = element_blank(),axis.text.x=element_blank())
  }
  if(!is.na(legend.ncol)){
    bar.plot=bar.plot+guides(fill=guide_legend(ncol=legend.ncol))
  }
  if(!is.na(y.breaks)){
    bar.plot=bar.plot+scale_y_continuous(limits = c(ymin,ymax),breaks = seq(ymin*10e4,ymax*10e4,y.breaks*10e4)/10e4)
  }
  if(flip){
    bar.plot=bar.plot+coord_flip()
  }
  
  return(bar.plot)
}

summarySE=function(data=NULL, measurevar, groupvars=NULL, na.rm=T,
                   conf.interval=0.95, .drop=TRUE,type="median") {
  
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  if(type=="median"){
    datac=ddply(data, groupvars, .drop=.drop,
                .fun = function(xx, col) {
                  c(N    = length2(xx[[col]], na.rm=na.rm),
                    median = median   (xx[[col]], na.rm=na.rm),
                    sd   = sd     (xx[[col]], na.rm=na.rm)
                  )
                },
                measurevar
    )
  }else if(type=="mean"){
    datac=ddply(data, groupvars, .drop=.drop,
                .fun = function(xx, col) {
                  c(N    = length2(xx[[col]], na.rm=na.rm),
                    mean = mean   (xx[[col]], na.rm=na.rm),
                    sd   = sd     (xx[[col]], na.rm=na.rm)
                  )
                },
                measurevar
    )
  }
    
  my.replace=measurevar
  names(my.replace)=type
  datac=plyr::rename(datac, my.replace)
  datac$se=datac$sd / sqrt(datac$N)  
  ciMult = ifelse(datac$N>1, qt(conf.interval/2 + .5, datac$N-1), 0)
  datac$ci = datac$se * ciMult
  
  return(datac)
}