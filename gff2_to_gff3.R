library(gdata)

filename <- "Brachypodium_distachyon.gff2"
outfile <- "Brachypodium_distachyon.gff3"

data <- as.matrix(read.table(filename, sep="\t"))
# convert the spaces to "=" in the attributes column
data[,9]=gsub(" ", "=", trim(data[,9]))

parents=matrix("",nrow=0,ncol=2)
colnames("element","value")
for i in 1:nrow(data) {
	# get the attributes
	attributes=unlist(strsplit(data[i,9],';'))
	# find the one that matches the type column (we need the id)
	pos = grep(data[i,3],attributes)
	
	# get the value 
	elementId = unlist(strsplit(attributes[pos],"="))[2]
	elementPosition=which(parents[,1]==data[i,3])
	if (length(elementPosition)>0) {
		parents[elementPosition] = elementId
	} else {
		parents = rbind(parents,c(data[i,3],elementId))
	}
	# set the id tag
	data[i,9]=gsub(data[i,3],"ID",data[i,9])
	
	# get the attributes
	attributes=unlist(strsplit(data[i,9],';'))
	# find the one that matches the type column (we need the id)
	pos = grep(data[i,3],attributes)
	
	grep("ID", data[i,9])
	if (data[i,3]) == "gene" {
		attributes=strsplit(data[i,9],c(';',' ','='))
		pos = grep("gene",data[i,9])
		
	}
}