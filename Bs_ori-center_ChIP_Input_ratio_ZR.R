#strain1
input1 = read.csv ("WGS_BWX5560_Xyl_clean.csv", header=T)
input1_reorder = c(input1$Value[2000001:4033458], input1$Value[1:1999542])
input1_reorder_last = c(input1$Value[1999543:2000000])
input1_reorder_norm = input1_reorder/  ((sum(input1_reorder, input1_reorder_last)/42)  / 1e6)
input1_reorder_last_norm = input1_reorder_last/  ((sum(input1_reorder, input1_reorder_last)/42)  / 1e6)
input1_bin = colSums(matrix(input1_reorder_norm,nrow=1000))/1000
input1_bin_last = sum(input1_reorder_last_norm)/458
input1_bin_whole = c(input1_bin, input1_bin_last)

########
strain1 = read.csv ("anti-mCherry_BWX5560_Xyl_clean.csv", header=T)
strain1_reorder = c(strain1$Value[2000001:4033458], strain1$Value[1:1999542])
strain1_reorder_last = c(strain1$Value[1999543:2000000])
strain1_reorder_norm = strain1_reorder/  ((sum(strain1_reorder, strain1_reorder_last)/42)  / 1e6)
strain1_reorder_last_norm = strain1_reorder_last/  ((sum(strain1_reorder, strain1_reorder_last)/42)  / 1e6)
strain1_bin = colSums(matrix(strain1_reorder_norm,nrow=1000))/1000
strain1_bin_last = sum(strain1_reorder_last_norm)/458
strain1_bin_whole = c(strain1_bin, strain1_bin_last)


#plotting
plot(seq(0, 4033, 1), strain1_bin_whole/input1_bin_whole, ylim=c(0, 10), type="h", lwd=1.5, axes=FALSE, xlab="", ylab="", col="red")

axis(1, at=c(0, 1000, 2000, 2033, 3033, 4033), labels=c(2000, 3000, 4000, 0, 1000, 2000))
axis(2)

#indicate ParS-1
#abline (v=2023589/1000, lty=1, lwd=2,col="green")
#indicate ParS-59 abline (v=1377806/1000, lty=1, lwd=2,col="green")

box()

title(ylab="ChIP/Input", xlab="genome position")

dev.copy2pdf(file="ChIP2input_Anti_mCherry_BWX5560_Xyl90min_ylim10.pdf", width = 10, height = 5)








