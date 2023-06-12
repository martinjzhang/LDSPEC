# Step1: create a bim file with -10kb, and a bim file with +10kb, and an empty ped file
R
for (chr in 1:22) {
    data_path='/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/baseline_annot/bim'
    bim_file=paste('ukb_imp_chr',chr,'_v3',sep='')
    out_path='/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/baseline_annot/recomb_rate'
    # Code
	data=read.table(paste(data_path,'/',bim_file,".bim",sep=""),h=F)
	for (window in c(10)) {
		temp=cbind(data[,1:2],0,pmax(data[,4]-(window*1000),min(data$V4)),data[,5:6])
		write.table(temp,file=paste(out_path,'/',bim_file,".inf.",window,"kb.bim",sep=""),quote=F,col.names=F,row.names=F)
		temp=cbind(data[,1:2],0,pmin(data[,4]+(window*1000),max(data$V4)),data[,5:6])
		write.table(temp,file=paste(out_path,'/',bim_file,".sup.",window,"kb.bim",sep=""),quote=F,col.names=F,row.names=F)
	}
    ped = c("ID1","ID1",0,0,0,0,rep(1,2*nrow(temp)))
    write.table(t(ped),file=paste(out_path,'/',bim_file,".ped",sep=""),quote=F,col.names=F,row.names=F)
}
q()
n

# Step2: add cM to the 2 bim files
for CHR in {1..22}; do
    out_path=/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/baseline_annot/recomb_rate
    bim_file=${out_path}/ukb_imp_chr${CHR}_v3
    # Code
    for BORN in inf sup; do
        /n/groups/price/steven/soft/plink2_v1.90b3w/plink \
        --map ${bim_file}.$BORN.10kb.bim \
        --ped ${bim_file}.ped \
        --cm-map /n/groups/price/HGb37_REFERENCE/genetic_map_chr${CHR}_combined_b37.txt $CHR \
        --zero-cms \
        --make-bed \
        --out ${bim_file}.oxford.$BORN.10kb.cM
        rm ${bim_file}.oxford.$BORN.10kb.cM.fam ${bim_file}.oxford.$BORN.10kb.cM.bed ${bim_file}.oxford.$BORN.10kb.cM.log
    done
done
rm ${out_path}/*ped ${out_path}/*nosex

# Step3: compute recombination rate
R
for (chr in 1:22) {
    out_path='/n/groups/price/martin/data_GDREG/UKBimp_337K_MAF001/baseline_annot/recomb_rate'
    bim_file=paste('ukb_imp_chr',chr,'_v3',sep='')
    # Code
    data1    = read.table(paste(out_path,'/',bim_file,".oxford.inf.10kb.cM.bim",sep=""),h=F)
    data2    = read.table(paste(out_path,'/',bim_file,".oxford.sup.10kb.cM.bim",sep=""),h=F)
    rec_rate = (data2$V3-data1$V3)/((data2$V4-data1$V4)/1000000) #cM/Mb
    write(rec_rate,file=paste(out_path,'/',bim_file,".oxford.10kb.recrate.annot",sep=""),ncolumns=1)
}
q()
n

rm ${out_path}/*nosex
rm ${out_path}/*bim
gzip *.annot