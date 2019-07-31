

args = commandArgs(trailingOnly=TRUE)

a = read.table("../data/acceptable_single/sites/1.stats.txt", T)

mid = 32
dif = 15

begin = mid - dif
finish = mid + dif

pdf("../results/substitution.pdf")

plot(a$A[begin:finish], type='l', col='blue', lwd=5, xlab="Time", ylab="High fitness probability")
lines(a$B[begin:finish], col='red', lwd=5)

dev.off()