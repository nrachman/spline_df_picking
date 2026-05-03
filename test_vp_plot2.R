library(variancePartition)
df <- data.frame(Subject = c(0.1, 0.2), CellFreqs = c(0.5, 0.6))
p <- plotVarPart(df)
print(p$data)
