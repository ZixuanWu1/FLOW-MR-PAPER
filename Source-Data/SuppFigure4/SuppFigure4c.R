library(corrplot)
x = read.csv("k=4_noise_correlation.csv")
row.names(x) = x[,1]
x = x[, 2:5]
corrplot(as.matrix(x), method = "circle", col= colorRampPalette(c('#313695',"white", '#A50026'))(10) )
