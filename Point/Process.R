setwd('~/Working/PTM_3D/Uniform')
source('~/Working/FlexEFT1D/Rscripts/Interpolate_WOA.R')
library(plot3D)
library(Rcpp)
sourceCpp("utils.cpp")
source('persp.R')
source('flat.R')
source('nndist.R')
oshima_persp()
oshima_persp(dir='~/Working/PTM_3D/Uniform')
oshima_persp(dir='~/Working/PTM_3D/Point2')
oshima_flat()
oshima_flat(dir='~/Working/PTM_3D/Uniform')
oshima_flat(dir='~/Working/PTM_3D/Point2')

system('ffmpeg -framerate 4 -i Output/Zoou%04d.png Output/Zoo_point2_20181010.avi')

system('ffmpeg -framerate 4 -i Output_persp/Zoo%04d.png Output_persp/Zoobath_point2_20181010.avi')

plot_nndist()
plot_nndist(dir='~/Working/PTM_3D/Uniform')
plot_nndist(dir='~/Working/PTM_3D/Point2')

trajfile  = paste0('Trajectory',DATE,'.pdf')
pdf(trajfile, width=8, height=8)
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(3.5,4,2,2),
             mgp     = c(2.3,1,0),
             cex.lab = 1.2,
             lwd     = 1.5,
             pch     = 16,
             mfrow   =c(1,1),
             cex.axis=1) 

#Calculate the individual trajectory
plot(Tout, parX[1,1,]/1E3, type='l',
     ylim=c(-50,10),
     xlab='Elapsed time (Days)', ylab='Horizontal position (km)')

legend('topleft', paste0('Group',1:4), col=1:4, lty=1)
for (i in 1:NG){
    for (j in 1:10){
        points(Tout, parX[i,j,]/1e3, type='l', col=i, lty=j)
    }
}
dev.off()

#Calculate DVM amplitude of each particle
final_day= (NT-ND):NT
dffz_min = apply(parZ[,,final_day], c(1,2), min)
dffz_max = apply(parZ[,,final_day], c(1,2), max)
dffz     = dffz_max - dffz_min
dffz     = apply(dffz, 1, mean)
dffz     = round(dffz,0)

#Calculate the horizontal moved distance of each particle 
dffx   = abs(parX[,,NT] - parX[,,1])/1E3
dff_m  = apply(dffx, 1, mean)
dff_sd = apply(dffx, 1, sd)

#Test whether statistically significant

t.test(dffx[1,], dffx[2,], paired = T)

#speed calculated from x positions
dtout = days_simulation*86400/(NT-1)
ucal  = array(NA, dim=c(NG, N_PAR, NT-1))
for (i in 2:NT){
    ucal[,,i-1] = (parX[,,i]-parX[,,(i-1)])/dtout
}

mean_ucal = apply(ucal, c(1,2), mean)

dispfile  = paste0('dispersion',DATE,'.pdf')
pdf(dispfile,width=5, height=9)
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(3.5,4,2,2),
             mgp     = c(2.3,1,0),
             cex.lab = 1.2,
             lwd     = 1.5,
             pch     = 16,
             mfrow   =c(2,1),
             cex.axis=1) 

boxplot(t(dffx), horizontal=T, 
        outline=F, yaxt='n', 
        xlab = 'Horizontally displaced distance (km)',
        ylab = 'Group'
        #ylab = 'Vertical migration amplitude (m)'
        )
axis(2, at = 1:length(dffz), labels=1:4)

boxplot(t(mean_ucal), horizontal=T, 
        outline=F, yaxt='n', 
        xlab = 'Average u experienced by each particle (m/s)',
        #ylab = 'Vertical migration amplitude (m)'
        ylab = 'Group'
        )
#axis(2, at = 1:length(dffz), labels=dffz)
axis(2, at = 1:length(dffz), labels=1:4)

dev.off()

#Find out the particle with the greatest displacement
imax = sapply(1:NG, function(i)which.max(dffx[i,]))
imin = sapply(1:NG, function(i)which.min(dffx[i,]))

op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(3.5,4,2,2),
             mgp     = c(2.3,1,0),
             cex.lab = 1.2,
             las     = 1,
             pch     = 16,
             mfrow   =c(1,1),
             cex.axis=1) 


plot(medianu[1,], dffx[1,], xlim=c(-0.05, 0.05), ylim=range(dffx), 
     xlab='Median u (m/s)', ylab='Horizontally displacement (km)') 
for (i in 2:NG){
    points(medianu[i,], dffx[i,], col=i)
}

#The distribution of u during one day

hist(parU[1, 1, (NT - ND):NT])


#Plot out the frequency distribution of u
pdf('histog_u.pdf',
                width=8, height=8,paper='a4')
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(3.5,4,2,2),
             mgp     = c(2.3,1,0),
             cex.lab = 1.2,
             las     = 1,
             pch     = 16,
             mfrow   =c(4,2),
             cex.axis=1) 

for (i in 1:NG){
    hist(parU[i,imax[i],], freq=F)
    hist(parU[i,imin[i],], freq=F)
}
dev.off()
#Check whether HDD relates to the average depth of each particle
#Calculate average depth of each particle
z_mean = apply(parZ, c(1,2), mean)
plot(z_mean[1,], dffx[1,], pch=16,
     xlab='Mean vertical position (m)',
     ylab='Horizontally displaced distance (km)')
for (i in 2:NG){
    points(z_mean[i,], dffx[i,], pch=16,col=i)
}

#Check the vertical position and u
wpar <- ncread(parfile,'w')
upar <- ncread(parfile,'u')
plot(parZ[1,1,], upar[1,1,], pch=16,
     xlab='Vertical position (m)',
     ylab='Horizontal velocity(m/s)')

#Obtain bathymetry for each particle
#Initial X position

#Read bathymetry file
bath = read.table('Depth.dat', header = T)

#Initial Bottom depth for the particles
parB = matrix(NA, nr = NG, nc = N_PAR)
for (i in 1:NG){
    parB[i,] = approx(x=bath[,1], y = bath[,2], xout = parX[i,,1])$y
}
parB = abs(parB)
pdf('Bathymetry_HDD.pdf',
                     width=6, height=6,paper='a4')
op <- par( font.lab  = 1,
             family  ="serif",
             mar     = c(3.5,4,2,2),
             mgp     = c(2.3,1,0),
             cex.lab = 1.2,
             lwd     = 1.5,
             pch     = 16,
             mfrow   =c(1,1),
             cex.axis=1) 


plot(as.numeric(parB), as.numeric(dffx), 
    type = 'n', xlab = 'Initial bottom depth (m)', ylab = 'Horizontal displaced distance (km)')
for (i in 1:NG){
    points(parB[i,], dffx[i,], pch=16, col = i)
}
txt = paste0('DVM=', dffz,'m')
legend('topleft', txt, pch=16, col=1:3 )
dev.off()


#Select particles with initial bottom depth shallower than 120 m

wB = which(parB < 120)
NW = length(wB)
hist(parU[1,wB[1],])
#setwd('~/Working/SUNTANS_2D_ridge/images')
#for (i in 24:71){
#    file1 = paste0('middle_diurnal_5mm_00',i,'.png')
#    k     = i-23
#    file2 = paste0('middle_diurnal_5mm_',sprintf('%02d',k),'.png')
#    file.copy(from = file1, to = file2)
#}

#system('ffmpeg -framerate 4 -i middle_diurnal_5mm_%04d.png out.avi')
