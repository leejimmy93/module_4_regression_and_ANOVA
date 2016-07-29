

Check_File_Exists<-function(FileName){
	
	if(!file.exists(FileName)){
		Msg<-sprintf("File %s does not exist\n",FileName)
		stop(Msg)
	}

}

#
#	File.Annovar : .variant_function file of ANNOVAR
#	File.BIM : BIM file
#	File.SetID : filename of set ID file will be created
#	Include : types of variants included in sets
#
GetSetID<-function(File.Annovar, File.BIM, File.SetID, Include=c("exonic", "splicing")){

	#File.Annovar<-"./ExomeChr22.avinput.variant_function";File.BIM<-"./ExomeChr22.bim";File.SetID<-"./ExomeChr22.SetID";Include=c("exonic", "splicing")
	
	Check_File_Exists(File.Annovar)
	Check_File_Exists(File.BIM)
	
	Anno<-read.table(File.Annovar, header=FALSE, stringsAsFactors=FALSE)
	BIM<-read.table(File.BIM, header=FALSE, stringsAsFactors=FALSE)

	ID.include<-NULL
	for(i in 1:length(Include)){
		for(j in 1:2){
			temp1<-NULL
			if(j==1){
				temp<-sprintf(";%s",Include[i])
			} else {
				temp<-sprintf("^%s",Include[i])
			}
			temp1<-grep(temp, Anno$V1)
			ID.include<-union(ID.include, temp1)		
		}
	
	}
	
	Anno1<-Anno[sort(ID.include),c(1:7)]
	colnames(Anno1)<-c("Function", "Gene", "Chr", "Pos", "PosEnd", "A1", "A2")
	colnames(BIM)<-c("Chr", "ID", "gp", "Pos", "BIM.A1", "BIM.A2")
	Anno1$PosID<-paste(Anno1$Chr, Anno1$Pos, sep="-")
	BIM$PosID<-paste(BIM$Chr, BIM$Pos, sep="-")	
	
	# LOG
	msg1=sprintf("%d variants in ANNOVAR file [%s] belong to groups in Include.\n", length(Anno1$PosID), File.Annovar )
	msg2=sprintf("BIM file [%s] has %d variants.\n", File.BIM, length(BIM$PosID))
	cat(msg1)
	cat(msg2)	
	
	Data.A<-merge(Anno1, BIM, by.x="PosID", by.y="PosID", all.x=TRUE, all.y=FALSE)
	
	# Get entries with only 1 gene 
	ID1<-grep("[;,(]", Data.A$Gene)
	if(length(ID1) == 0){
		
		SetID<-cbind(SetID=Data.A$Gene, SNPID=Data.A$ID)
		
	} else {
		ID2<-(1:length(Data.A$Gene))[-ID1]
		SetID<-cbind(SetID=Data.A$Gene[ID2], SNPID=Data.A$ID[ID2])
		SNPID<-Data.A$ID[ID1]
			
		temp.s<-strsplit(Data.A$Gene[ID1], "[,(;]+")
		for(i in 1:length(ID1)){
			id1<-grep("N*_", temp.s[[i]])
			if(length(id1)> 0){
				temp.s1<-unique(temp.s[[i]][-id1])
			} else {
				temp.s1<-unique(temp.s[[i]])
			}
			
			SetID<-rbind(SetID, cbind(temp.s1, SNPID[i]))
		}
		
	}
	
	# sort by gene name 
	SetID<-SetID[order(SetID[,1]),]
	
	write.table(SetID, File.SetID, col.names=FALSE, row.names=FALSE, quote=FALSE)

	# LOG
	msg1=sprintf("%d genes (SNP sets) and total %d variants are saved in SetID file.\n", length(unique(SetID[,1])), length(SetID[,1]))
	cat(msg1)
	
}


