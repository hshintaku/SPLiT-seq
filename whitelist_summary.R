f_whitelist <- "/home/samba/storage0/shintaku/20200919MiSeq012/P10PCR01_S1_L001.whitelist.txt"
whitelist <- read.table(f_whitelist, sep="\t")
colnames(whitelist) <- c('represent','variants', 'total','counts')
whitelist_sorted <- whitelist[order(whitelist$total,decreasing = T),]
whitelist_sorted$cumsum <- cumsum(whitelist_sorted$total)/sum(whitelist_sorted$total)
plot(whitelist_sorted$cumsum, xlim=c(0,20000))


barcode <- read.table("/home/samba/storage0/shintaku/SPLiT-seq/SPLiTbarcode2.txt")

encoded <- whitelist.umi_tools.list(whitelist,barcode$V1)


encoded2nd <- whitelist.umi_tools.encode(whitelist$represent,barcode$V1)
encoded2nd$count <- whitelist$total

counts <-data.frame()
for (icnt in 1:25){
  counts[icnt,1] <- colSums(data.frame(encoded[which(encoded$value==icnt-1),]$count))
}

plot(counts$V1)

hist(encoded[which(encoded$second_index==11),]$second)

index <- data.frame(which(encoded$lv_total<=8))
colnames(index) <- c("matched")
matched <- encoded[index$matched,]
hist(matched$count,30)


cells_1000counts <- matched[which(matched$count>1000),]
colSums(cells_1000counts,2)
colSums(matched,2)

rt <- data.frame()
for (icnt in 1:16){
  rt[icnt,1] <- colSums(data.frame(matched[which(matched$first_index==icnt),]$count))
}
barplot(rt$V1)

