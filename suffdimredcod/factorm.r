library("MASS")
#library("expm")

factorm <- function(x){
svdx <- eigen(x%*%t(x))

hatF <- svdx$vector[, 1:k] * sqrt(n)
hatB  <- n^(-1)*t(x) %*% hatF
sB <- sign(hatB[1, ])
hatF<- t(apply(hatF, 1, '*', sB ))
cF <- cov(hatF )
svdcF <- svd(cF)
sovleF <- svdcF$u %*% diag((svdcF$d)^(-1/2))%*% t(svdcF$v)
hatF <- (hatF- matrix(apply(hatF, 2, mean), n, k, byrow = T)) %*% sovleF
return(hatF)


}
write.table(resmat, file = "datasrd.tab", row.names = FALSE, col.names = FALSE)
write.table(resmat3, file = "datasrd3.tab", row.names = FALSE, col.names = FALSE)
write.table(resmat4, file = "datasrd4.tab", row.names = FALSE,
col.names = FALSE)
write.table(resmat6, file = "datasrd6.tab", row.names = FALSE,
col.names = FALSE)
write.table(resmat7, file = "datasrd7.tab", row.names = FALSE,
col.names = FALSE)
write.table(beta, file = "beta.tab", row.names = FALSE, col.names = FALSE)


