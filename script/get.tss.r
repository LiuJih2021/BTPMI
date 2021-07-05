positive<-read.table("positive.gtf",sep = "\t",quote = "")
TSS<-positive$V4
end<-positive$V5
promoter<-TSS-500
end2<-c(0,end)
for(i in 1:length(promoter)){if(promoter[i]<end2[i]){ promoter[i]=end2[i]}}
promoter<-cbind(promoter,TSS+100)
write.table(promoter,"tss.po.bed",sep = "\t", row.names =FALSE, col.names =FALSE,quote = FALSE)

negative<-read.table("negative.gtf",sep = "\t",quote = "")
 TSS<-negative$V4
 end<-negative$V5
 promoter<-end+500
 for(i in 1:(length(promoter)-1)){if(promoter[i]>TSS[i+1])promoter[i]=TSS[i+1]}
promoter<-cbind(end-100,promoter)
 write.table(promoter,"tss.ne.bed",sep = "\t", row.names =FALSE, col.names =FALSE,quote = FALSE)


