library(corrplot)
x = read.csv("/Users/zixuanwu/Downloads/research/mediation_source/SuppFigure5/multivariate_noirse_correlation.csv")
row.names(x) = colnames(x)
corrplot(as.matrix(x), method = "circle", col= colorRampPalette(c('#313695',"white", '#A50026'))(10) )
