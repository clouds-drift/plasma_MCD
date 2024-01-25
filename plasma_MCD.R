library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("ggplotify")
library("openxlsx")
library("GenomicRanges")
library("rtracklayer")
library("data.table")
library("BSgenome.Hsapiens.UCSC.hg19")
library("reshape2")
library("dplyr")
source("MCD_fun.r")
##input
N_thread=12
N_model=10
genome="hg19"
bw.norm.dir="./bw_norm_count/"

p.name="pvalue"
p.cut=c(0.01)
lfc.cut=1
imp.cut=0
ntop=c(Inf)
diff.dir="./QSEA_diff/215_training_set_other_nature_block_QSEA"
cutoff1.tb=expand.grid(p.name,p.cut,lfc.cut,ntop,imp.cut)
#feature1.saf="./QSEA_diff/215_training_set_other_nature_block_QSEA/all_sample/pvalue0.01_LFC1/total200.saf"
feature1.bed="./QSEA_diff/215_training_set_other_nature_block_QSEA/all_sample/pvalue0.01_LFC1/total200.bed"
if(any(!file.exists(feature1.bed))){stop(feature1.bed,"not found!\n")}
feature1.name="pvalue0.01_LFC1_total200"

bw.bias.dir="./bw_bias"
cov.cut=1 ##RPM cut
    #2 ##read cut
names(cov.cut)="RPM"
    #"read"
miss.percent=1
p.name="pvalue"
p.cut=0.01
delta.cut=0.3
base.cut=0.3
imp.cut=0
ntop=Inf
value.type=c("raw")
bias.trans=F
diff.dir="./data/DMBR_diff/215_training_set_other_nature_block_RPM1_each0.3_miss1"
    #"./DMBR_diff/175_training_set_other_nature_block_read2_miss1"
cutoff2.tb=expand.grid(p.name, p.cut, delta.cut, base.cut, ntop, value.type, imp.cut)
feature2.saf="./DMBR_diff/215_training_set_other_nature_block_RPM1_each0.3_miss1/all_sample/pvalue0.01_delta0.3_base0.3/total200.saf"
feature2.bed="./DMBR_diff/215_training_set_other_nature_block_RPM1_each0.3_miss1/all_sample/pvalue0.01_delta0.3_base0.3/total200.bed"
if(any(!file.exists(feature2.bed))){stop(feature2.bed,"not found!\n")}
feature2.name="RPM1_each0.3_pvalue0.01_delta0.3_base0.3_total200"
model.name="GLMNET_model"
cali.method="GLMNET_model"
force=F
##output
DMR.model.dir="./DMR_model"
DHMR.model.dir="./DHMR_model"
Calibration.model.dir="./Calibration_model"
Calibration.pred.dir="./Calibration_prediction"





cat("Prepare main model samples\n")
sample.df=read.xlsx("./215_training_set.xlsx")
discover.df=format_sample_table(file=sample.df[,c("Sample","PatientID","Disease","IDH1","Filter_reads","stage")])
sample.df=read.xlsx("./215_training_set.xlsx")
test.df=format_sample_table(file=sample.df[,c("Sample","PatientID","Disease","IDH1","Filter_reads")])
sample.df=read.xlsx("./56_validation_set.xlsx")
vali.df=format_sample_table(file=sample.df[,c("Sample","PatientID","Disease","IDH1","Filter_reads","stage")])
main.level=sort(unique(discover.df$Label))
detail.level=sort(unique(discover.df$Detail_Label))
info.list=list(train=discover.df, test=test.df, vali=vali.df)
data.list=get_file_list(info.list, bw.norm.dir, bw.bias.dir)



##main DMR model
project.name=paste(N_model,"model_",nrow(discover.df),"_DMR_",feature1.name,sep="")
main.dmr.model.dir=paste(DMR.model.dir,"/",project.name,"/Train_by",nrow(discover.df), sep="")

##main DHMR model
project.name=paste(N_model,"model_",nrow(discover.df),"_DHMR_",feature2.name,sep="") 
main.dhmr.model.dir=paste(DHMR.model.dir,"/",project.name,"/Train_by",nrow(discover.df), sep="")





cali.model.dir=paste(Calibration.model.dir,"/", basename(dirname(main.dmr.model.dir)), 
    "_and_", basename(dirname(main.dhmr.model.dir)),
    "/Calibrate_model_by",nrow(discover.df),"_",
    sub("_.*","",cali.method),
    "_",value.type,sep="")
pred.dir=paste(Calibration.pred.dir,"/",basename(dirname(cali.model.dir)), sep="")

cat("Evaluate DMR vali...\n")
vali1.sample=names(data.list$vali$Name)
names(vali1.sample)=names(data.list$vali$Name)
value.type="raw"
dens.trans=F
dmr.pred.dir=paste(pred.dir,
    "/DMR_model_pred",length(vali1.sample),
    "_",value.type,
    ifelse(dens.trans,"_zdens",""),
    sep="")
model1=main.dmr.model.dir
vali.label=data.list$vali$Label
vali.name=data.list$vali$Name
pred.level=main.level
main.dmr.mid.dir=paste(dirname(model1),"/score",sep="")
dmr.vali.list=sample_pred(my.sample=vali1.sample, input.label=vali.label, pred.level=pred.level, 
    one.model.dir=model1, 
    out.dir=dmr.pred.dir, 
    mid.dir=main.dmr.mid.dir, 
    value.type = value.type, transform=dens.trans, 
    N_model=10, model.name=model.name, N_thread=N_thread, force=F, rm.ref=F)
cat("\n\nShow DMR vali performance\n")
evaluate_pred(result.list=dmr.vali.list, input.label=vali.label, input.name=vali.name, 
    normal.cut=NA,out.dir=dmr.pred.dir, my.color=c("red","orange","blue"))


cat("Evaluate DHMR vali...\n")
vali2.sample=names(data.list$vali$Name)
names(vali2.sample)=names(data.list$vali$Name)
bias.trans=F
dhmr.pred.dir=paste(pred.dir,
    "/DHMR_model_pred",length(vali2.sample),
    "_",value.type,
    ifelse(bias.trans,"_zbias", ""), 
    sep="")
model2=main.dhmr.model.dir
main.dhmr.mid.dir=paste(dirname(model2),"/score",sep="")
dhmr.vali.list=sample_pred(my.sample=vali2.sample, input.label=vali.label, pred.level=pred.level, 
    one.model.dir=model2, 
    out.dir=dhmr.pred.dir, 
    mid.dir=main.dhmr.mid.dir, 
    value.type = value.type, transform=bias.trans, 
    N_model=10, model.name=model.name, N_thread=N_thread, force=F, rm.ref=F)
cat("\n\nShow DHMR vali performance\n")
evaluate_pred(result.list=dhmr.vali.list, 
    input.label=vali.label, input.name=vali.name, 
    normal.cut=NA, out.dir=dhmr.pred.dir, my.color=c("red","orange","blue"))


cat("Evaluate calibration vali...\n")
cali.pred.dir=paste(pred.dir,
    "/Calibrate_model_vali", length(vali1.sample),
    "_",paste(sub("_.*","",model.name),collapse = "_"),
    "_",value.type,
    ifelse(dens.trans,"_zdens",""),
    ifelse(bias.trans,"_zbias",""),
    sep="")
vali.prob.mat=combine_problist(dmr.vali.list, dhmr.vali.list, prefix1="Dens",prefix2="Bias",train.level = pred.level)
cali.mod.list=read_my_model(cali.model.dir,model.name=cali.method,N_model=NULL)
cali.vali.list=my_model_pred(cali.mod.list, mod.name=cali.method, input.mat=t(vali.prob.mat), pred.level=pred.level, input.label=vali.label)
cat("\n\nShow calibration vali performance\n")
evaluate_pred(result.list=list(calibration=cali.vali.list), 
        input.label=vali.label, input.name=vali.name, 
        normal.cut=NA,out.dir=cali.pred.dir, my.color=c("red","orange","blue"))

