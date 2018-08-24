cdSc <- readRDS('atac_matrix.binary.qc_filtered.rds')
dim(cdSc)

index <- sample(1:dim(cdSc)[2], 3000)

cd <- as.matrix(cdSc[,index])

cell_metadata <- read.table('cell_metadata.txt',sep = '\t', header = TRUE)
cell_metadata_sample <- cell_metadata[index,]

write.table(cd, 'mouse_atlast.csv', sep='\t')
write.table(cell_metadata_sample, 'cell_metadata_sample.csv', sep='\t')
