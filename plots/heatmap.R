
rm(list=ls())

source(paste(getwd(),"/environment.R",sep=""))
source(paste(getwd(),"/functions/MC_simulation_functions.R",sep=""))

#options(tz="CA")

set.seed(184)

#model parameters

N <- 20
p <- 100
p1 <- 40
p2 <- 25
X1_shift <- c(-1,1)
X2_shift <- c(-2,2)
sigma2_0 <- 1
sigma2_1 <- 1.5^2
beta<-2

train_sim <- gene_simulation(N,p,p1,p2,X1_shift,X2_shift,beta,sigma2_0,sigma2_1)

X <- train_sim$X

filename <- paste(figures_dir,"/heatmap_png.png",sep="")
png(filename,width=130,height=100,units="mm",res=300)
#tikz(file=filename,width = 4, height = 3)

my_colors <- colorRampPalette(c("green", "red"))
plot <- heatmap(t(X),scale="none",Rowv=NA,Colv=NA,col=my_colors(100),labRow = NA,labCol = NA,revC = TRUE)

text(0.35,-0.15,labels="DNA samples",adj=0.5,xpd=TRUE)
text(0.9,0.6,labels="Genes",adj=0.5,xpd=TRUE)

legend("bottomright", legend=c("min","max"), fill=my_colors(2))

print(plot)

dev.off()

       