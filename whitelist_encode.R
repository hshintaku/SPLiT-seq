whitelist.umi_tools.encode <- function (extracted,barcode){
  library(stringdist)
  distancematrix <- stringdistmatrix(extracted,barcode,method="lv")
  romin <- data.frame(apply(distancematrix,1,min))
  index <- data.frame(apply(distancematrix,1,which.min))
  romin <- cbind(romin,index)
  colnames(romin)<-c("value","index")
  return(romin)
}
whitelist.umi_tools.list <- function(whitelist,barcode){
  extracted1 <- data.frame(substr(whitelist$represent, 1,8))
  extracted2 <- data.frame(substr(whitelist$represent, 9,16))
  extracted3 <- data.frame(substr(whitelist$represent, 17,24))
  extracted <- data.frame(cbind(extracted1,extracted2,extracted3))
  colnames(extracted) <- c("first","second","third")
  extracted$count <- whitelist$total
  
  romin <- whitelist.umi_tools.encode(extracted$first,barcode)
  encoded1 <- romin
  romin <- whitelist.umi_tools.encode(extracted$second,barcode)
  encoded2 <- romin
  romin <- whitelist.umi_tools.encode(extracted$third,barcode)
  encoded3 <- romin
  total <- encoded1$value+encoded2$value+encoded3$value
  encoded <- data.frame(cbind(extracted1,encoded1,extracted2,encoded2,extracted3,encoded3,total,extracted$count))
  colnames(encoded)<- c("first_barcode","first","first_index","second_barcode","second","second_index","third_barcode","third","third_index","lv_total","count")
  return(encoded)

}