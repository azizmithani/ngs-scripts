library(gdata)

filename="Brachypodium_distachyon.gff"
out.filename="Brachypodium_distachyon_igv.gff"
data=trim(as.matrix(read.table(filename, sep="\t")))

rowsGenes = which(data[,3]=="gene")

for (i in 1:length(rowsGenes)) {
	row=rowsGenes[i]
	attributes = unlist(strsplit(data[row,9],";"))
	pos = grep("ID=",attributes)
	id=unlist(strsplit(attributes[pos],"="))[2]
	#data[row,9] = paste("gene ",id,";",data[row,9],sep="")
	data[row,9] = paste(data[row,9],";gene ",id,sep="")
}

write.table(data,out.filename,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
