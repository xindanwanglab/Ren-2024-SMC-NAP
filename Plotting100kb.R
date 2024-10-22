


#strain1
input1 = read.csv ("WGS_BWX523_clean.csv", header=T)
input1_reorder = c(input1$Value[2000001:4033458], input1$Value[1:1999942])
input1_reorder_last = c(input1$Value[1999943:2000000])
input1_reorder_norm = input1_reorder/  ((sum(input1_reorder, input1_reorder_last)/42)  / 1e6)
input1_reorder_last_norm = input1_reorder_last/  ((sum(input1_reorder, input1_reorder_last)/42)  / 1e6)
input1_bin = colSums(matrix(input1_reorder_norm,nrow=100))/100
input1_bin_last = sum(input1_reorder_last_norm)/58
input1_bin_whole = c(input1_bin, input1_bin_last)

########
strain1 = read.csv ("Anti_HbsU_BWX523_clean.csv", header=T)
strain1_reorder = c(strain1$Value[2000001:4033458], strain1$Value[1:1999942])
strain1_reorder_last = c(strain1$Value[1999943:2000000])
strain1_reorder_norm = strain1_reorder/  ((sum(strain1_reorder, strain1_reorder_last)/42)  / 1e6)
strain1_reorder_last_norm = strain1_reorder_last/  ((sum(strain1_reorder, strain1_reorder_last)/42)  / 1e6)
strain1_bin = colSums(matrix(strain1_reorder_norm,nrow=100))/100
strain1_bin_last = sum(strain1_reorder_last_norm)/58
strain1_bin_whole = c(strain1_bin, strain1_bin_last)

#plotting
plot(seq(0, 40334, 1), strain1_bin_whole/input1_bin_whole,xlim=c(10510,11500), ylim=c(0, 10), type="l", lwd=1.5, axes=FALSE, xlab="", ylab="", col="red")

axis(1, at=c(10500, 10700, 10900, 11100, 11300, 11500), labels=c(3050, 3070, 3090, 3110, 3130, 3150))
axis(2)

#indicate ParS-1
abline (h=1, lty=2, lwd=2,col="black")
#indicate ParS-59 abline (v=1377806/1000, lty=1, lwd=2,col="green")

box()

title(ylab="ChIP/Input", xlab="genome position")

dev.copy2pdf(file="ChIP2input_Anti_HbsU_BWX523_100kb3050_ylim10.pdf", width = 10, height = 5)

