options(echo=FALSE);
param <- commandArgs(trailingOnly=T)
chr  = eval(paste(text=param[1]))
bim  = eval(paste(text=param[2]))
output_path = eval(paste(text=param[3]))

print(paste(bim,chr))

mydata    = read.table(bim,h=F)
dataUK10K = read.table(paste("/n/groups/price/ldsc/reference_files/UK10K/plink_files/ALSPAC.TWINSUK.QC.",chr,".bim",sep=""),h=F) 

start10kb   = pmax(mydata[,4]-(10*1000),min(mydata$V4))
end10kb     = pmin(mydata[,4]+(10*1000),max(mydata$V4))

diversity10kb=NULL

cpt=1; data10Mb=subset(dataUK10K, (dataUK10K$V4 > (cpt-1)*5000000+1) & (dataUK10K$V4<(cpt+1)*5000000) ) 

for (i in 1:nrow(mydata)){
	while ( (mydata$V4[i]+1000*1000) > ((cpt+1)*5000000) ) {
		cpt=cpt+1; data10Mb=subset(dataUK10K, (dataUK10K$V4 > (cpt-1)*5000000+1) & (dataUK10K$V4<(cpt+1)*5000000) ) 
	}
	datatemp=subset(data10Mb,data10Mb$V4>(mydata$V4[i]-(10*1000)) & data10Mb$V4<(mydata$V4[i]+(10*1000)))
	diversity10kb=c(diversity10kb,nrow(datatemp))
}

write(diversity10kb  /((end10kb  -start10kb  )/1000),
      file=paste(output_path, "/diversity.10kb.",chr,".annot",sep=""),ncolumn=1)
