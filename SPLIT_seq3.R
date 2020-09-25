require(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

#Data Rename
myData <- read.table("/home/samba/storage0/Shiomi/20200919MiSeq012/count/P7PCR03_S19_L001.counts.tsv.gz", header = TRUE)
myData1 <- dplyr::mutate_all(myData, ~gsub(.,pattern="mm10.*",replacement = "mouse"))
myData1 <- dplyr::mutate_all(myData1, ~gsub(.,pattern="GRCh38.*",replacement = "human"))
myData1 <- dplyr::mutate_all(myData1, ~gsub(.,pattern="ENSS.*",replacement = "pig"))

myData$species <- myData1$gene

#Data summarize
grouped_data <- group_by(myData, cell,species)


Res <- summarize(grouped_data, count = n())
Res <- tidyr::spread(Res, key = species, value = count)

Res$represent <- Res$cell
Res$human <- replace(Res$human, which(is.na(Res$human)), 0)
Res$mouse <- replace(Res$mouse, which(is.na(Res$mouse)), 0)
Res$pig <- replace(Res$pig, which(is.na(Res$pig)), 0)

Res$total <- Res$human+Res$mouse+Res$pig
encoded <- whitelist.umi_tools.list(Res,barcode$V1)
index <- data.frame(which(encoded$lv_total<=1))
colnames(index) <- c("matched")
matched <- encoded[index$matched,]

rt <- data.frame()
for (icnt in 1:16){
  rt[icnt,1] <- colSums(data.frame(matched[which(matched$third_index==icnt),]$count))
}
barplot(rt$V1)

#plot
gp<-ggplot(Res,aes(x=pig,y=human))
gp=gp+geom_point()
print(gp)

