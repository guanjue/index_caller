####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
count_file = args[1]
output_filename = args[2]
col_num = as.numeric(args[3])

data = read.table(count_file, header=F)

print(head(data))

counts = data[,col_num]

png(output_filename)
plot(c(1:length(counts)), counts, type='l')
dev.off()