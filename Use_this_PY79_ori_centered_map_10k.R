
args = commandArgs(trailingOnly=TRUE)
outname = gsub('overlap.hm_corrected','',args[1])
r=read.table(args[1])


#k is the cutoff value
#
for (k in c(0.01,0.02,0.03)){

#for (k in c(0.01)){

#Define lower and upper cut-off
zlower = 0
zupper = k

#Function to trim the matrix down to min and max cut-off points
image_mod = function(m, zlower, zupper){
norm_m = m
for(i in 1:nrow(m) ){
  for (j in 1:ncol(m) ){
	if(m[i,j] >= zupper){
		norm_m[i,j] = zupper
	}
	if(m[i,j] <= zlower){
		norm_m[i,j] = zlower
	}
	if(m[i,j] < zupper & m[i,j] > zlower){
		norm_m[i,j] = m[i,j]
	}
  }
}
return(norm_m)
}
#######################################
#White blue purple black linear scheme
colvec= rgb(
rbind(
c(255, 255, 255)/255 ,
c(162, 192, 222)/255 ,
c(140, 137, 187)/255 ,
c(140, 87, 167)/255 ,
c(140, 45, 143)/255 ,
c(120, 20, 120)/255 ,
c(90, 15, 90)/255 ,
c(60, 10, 60)/255 ,
c(30, 5, 30)/255 ,
c(0, 0, 0)/255
)
)

#r=read.table("BWX4310_42C45m_30C15m_upIPTG15m-10k_overlap.hm_corrected")
#r=read.table("BWX4310_42C45m_30C15m_upIPTG20m-10k_overlap.hm_corrected")
#r=read.table("BWX4310_42C45m_30C15m_upIPTG25m-10k_overlap.hm_corrected")
#r=read.table("BWX4310_42C45m_30C15m_upIPTG30m-10k_overlap.hm_corrected")

m = as.matrix(r)
#Fill in the diagonal and +1/-1 subdiagonal with max cut-off point
diag(m) = zupper
diag(m[-nrow(m), -1]) = zupper
diag(m[-1, -nrow(m)]) = zupper
norm_m = m

#Generate the matrix with appropriate cutoff value
n= image_mod(norm_m, zlower, zupper)
# if 10kb bin, then the dimension of n is 404x404

#generate the 2by2 then 4by4 matrix
g = cbind(rbind(n,n), rbind(n,n))

#now need to cut this matrix to put ori in the middle of the map
#bin1: 0-10, bin2:10-20, bin3:20-30
starting_bin = 201
h = g[(starting_bin): (starting_bin + nrow(m)-1), (starting_bin): (starting_bin + nrow(m)-1)]

#Actual plotting
#pdf(file= paste ("~/Desktop/output/BWX4310_42C45m_30C15m_upIPTG15m-10k_", k, ".pdf", sep=""), useDingbats=FALSE)
#pdf(file= paste ("~/Desktop/output/BWX4310_42C45m_30C15m_upIPTG20m-10k_", k, ".pdf", sep=""), useDingbats=FALSE)
#pdf(file= paste ("~/Desktop/output/BWX4310_42C45m_30C15m_upIPTG25m-10k_", k, ".pdf", sep=""), useDingbats=FALSE)
#
pdf(file= paste (outname, k, ".pdf", sep=""), useDingbats=FALSE)


par(pty = "s")
image(seq(0, 4033459/1000, 10), seq(0, 4033459/1000,10),h, col=colorRampPalette(colvec)(256), zlim=c(zlower, zupper),xlab="genome position (kb)", ylab="genome position (kb)", useRaster=TRUE,axes=FALSE, xaxs="i")
axis(1, at=c(0, 1000, 2000,2040, 3000, 4000), labels=c(2000, 3000, 4000,0, 967, 1967))
axis(2, at=c(0, 1000, 2000,2040, 3000, 4000), labels=c(2000, 3000, 4000,0, 967, 1967))
axis(3, at=c(0, 1000, 2000,2040, 3000, 4000), labels=c(2000, 3000, 4000,0, 967, 1967))
axis(4, at=c(0, 1000, 2000,2040, 3000, 4000), labels=c(2000, 3000, 4000,0, 967, 1967))
box()
dev.off()
}




