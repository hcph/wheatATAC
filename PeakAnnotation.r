
#peakanotaton

args<-commandArgs(T)
dir = args[1]
num = length(args)-1

name=args[2:num]
#plot
pdf("Peak_Annotation.pdf",width=10)
for (i in name){
  x<-read.table(paste0(i,"_gene.peaks.txt"))
  head(x)
  x$fz<-abs(x[,6])
  x<-x[order(x[,1],x[,4],x[,7]),]
  x<-x[!duplicated(x[,4]),]
  x<-x[,-7]
  x<-x[order(x[,1],x[,2],x[,6]),]
  genebody<-length(which((x[,5]==0) & (x[,6]>1500)))
  gene_downstream<-length(which((x[,5] != 0) & (x[,6]>0)))
  promoter <- length(which(((x[,5]!=0) & (x[,6]>= -4000) & (x[,6] < 0)) | ((x[,5]==0) &x[,6] <= 1500)))
  intergeric <- length(which((x[,5]!=0) & (x[,6]< -4000)))
  re<-t(data.frame(genebody, gene_downstream, promoter, intergeric))
  t<-as.matrix(re)
  names<-c("genebody", "intergenic_downstream", "promoter (-5 kb-1.5kb)", "intergeric_upstream")
  # 计算百分比
  piepercent = paste(round(100*t/sum(t)), "%")
  #绘图
  #png(paste0(i,"_gene_annotation.png"), width=1000,height=700)
  pie(t, label=piepercent , main= paste0(i,"_annotation"), col=c("#ED1C24","#22B14C","#FFC90E","#3f48CC"))
  legend("topright", names, cex=1.2, fill=cols)
  #dev.off()
}
dev.off()
