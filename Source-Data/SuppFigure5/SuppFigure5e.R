library(corrplot)
x = read.csv("multivariate_noise_correlation.csv")
row.names(x) = x[,1]
x= x[, 2:8]
colnames(x) = row.names(x)
corrplot(as.matrix(x), method = "circle", col= colorRampPalette(c('#313695',"white", '#A50026'))(10) )
