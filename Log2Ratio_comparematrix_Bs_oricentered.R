#read matrices
m1 = as.matrix(read.table("BWX5567_xyl90min_DpnII-10k.hm_corrected"));
m2 = as.matrix(read.table("BWX3370_xyl90min_DpnII-10k.hm_corrected"));

#manual normalization of the two matrices
#for corrected matrices no need to normalize
#reads_m1 = 7492395;
#reads_m2 = 3149857;
reads_m1 = 1;
reads_m2 = 1;
nm1 = m1/reads_m1;
nm2 = m2/reads_m2;

#calculate log2ratio of two matrices
log2mat = log2(nm1/nm2);
m = log2mat
#replace NA with 0
m[is.na(m)] = 0
#Function to trim the matrix down to min and max cut-off points
zlower = -1
zupper = 1
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

#Generate the matrix with appropriate cutoff value
n= image_mod(m, zlower, zupper)
# if 10kb bin, then the dimension of n is 404x404

#generate the 2by2 then 4by4 matrix
g = cbind(rbind(n,n), rbind(n,n))

#now need to cut this matrix to put ori in the middle of the map
#bin1: 0-10, bin2:10-20, bin3:20-30
starting_bin = 201
h = g[(starting_bin-1): (starting_bin + nrow(m)-2), (starting_bin-1): (starting_bin + nrow(m)-2)]


#plot
pdf(file="./log2BWX5567xyl90mintoBWX3370xyl90min_-1_1_10k.pdf", useDingbats=FALSE)
par(pty = "s")
image(seq(0, 4036854/1000, 10), seq(0, 4036854/1000,10), h, col=colorRampPalette(c("blue","white","red"))(256), zlim=c(zlower, zupper),xlab="Genome position (kb)", ylab="Genome position (kb)", useRaster=TRUE,axes=FALSE, xaxs="i")
axis(1, at=c(0, 1000, 2000,2033, 3000, 4000), labels=c(2000, 3000, 4000,0, 967, 1967))
axis(2, at=c(0, 1000, 2000,2033, 3000, 4000), labels=c(2000, 3000, 4000,0, 967, 1967))
axis(3, at=c(0, 1000, 2000,2033, 3000, 4000), labels=c(2000, 3000, 4000,0, 967, 1967))
axis(4, at=c(0, 1000, 2000,2033, 3000, 4000), labels=c(2000, 3000, 4000,0, 967, 1967))
box()

dev.off()



