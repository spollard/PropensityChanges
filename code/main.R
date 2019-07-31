args = commandArgs(trailingOnly=TRUE)

#sites_dir = args[0]
#sites = args[1]

a = read.table("../data/acceptable/model", T)
b = apply(a, 2, mean)
for (d in names(b)) {cat(paste0(">",d, "\n",b[d],"\n"))}