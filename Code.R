library(gdata)
library(raster)
library(ggplot2)
library(FactoMineR)
library(corrplot)
library(factoextra)
library(RColorBrewer)
library(geomorph)
library(sp)
library(ncdf4)

#### 1. Load and prepare data ####

# Load climate data
load(".../Data.RData")

# Vectors for periods and variables name
P <- c("00k", "06k", "09k", "40k", "60k", "66k", "115k")
var <- c("drysoil_frac", "humtot", "hyd_stress", "precip", "qsurf", "sols",
         "tsol", "tz_ampl_day", "w10m")

# Extract values on EH2
for (p in 1:length(P)){
  for (v in 1:length(var)){
    rast <- get(apropos(P[p])[v])
    rast <- rast[46:47,70:71,] # extract values on EH2
    mv(from="rast", to=paste(P[p], sep="_", var[v]))
  }
}

# Standardize data between periods
for (v in 1:length(var)){
  # apropos(var[v])
  var_data <- NULL
  for (i in 1:length(P)){
    var_data <- rbind(var_data, get(apropos(paste("_",sep="",var[v]))[i]))
  }
  var_data_norm <- (var_data-mean(var_data))/sd(var_data)
  for (i in 1:length(P)){
    k=i*4-3;j=i*4
    x <- var_data_norm[k:j,]
    mv(from="x", to=apropos(paste("_",sep="",var[v]))[i])
  }
}

# Mean and sd per period
clim <- NULL
for (p in 1:length(P)){
  apropos(P[p])
  r <- NULL
  for (v in 1:length(var)){
    m <- mean(get(apropos(P[p])[v]))
    sd <- sd(get(apropos(P[p])[v]))
    r <- c(r, m, sd)
  }
  clim <- rbind(clim,r)
}
rownames(clim) <- P
colnames(clim) <- c("drysoil_frac_m", "drysoil_frac_sd",
                    "humtot_m", "humtot_sd", "hyd_stress_m", "hyd_stress_sd",
                    "precip_m","precip_sd", "qsurf_m", "qsurf_sd",
                    "sols_m", "sols_sd", "tsol_m", "tsol_sd",
                    "tz_ampl_day_m", "tz_ampl_day_sd", "w10m_m", "w10m_sd")
clim

# D1 and D2
D1 <- clim[c(1,2,3,4,4,5,7),]
# L0=00k, L1c=06k, L2c=09k, L3esr=40k, L4esr=40k, L5esr=60k, L8esr=115k
D2 <- clim[c(1,2,3,5,6,7,7,7,7),]
# L0=00k, L1c=06k, L2c=09k, L3osl=60k, L4osl=66k, L5osl=115k, L6osl=115k, L7osl=115k, L8osl=115k

#________________________________________________________________
#________________________________________________________________

#
#### 2. Covariance analysis ####

#________________________________
# PCA

# pca
pca <- PCA(D1, ncp=20) # change for D2
summary(pca)
# Plot
plot.PCA(pca, axes=c(1, 2), choix="ind")
# principal components contributions
fviz_contrib(pca, choice = "var", axes = 1, top = 20)
fviz_contrib(pca, choice = "var", axes = 2, top = 20)

# Biplot
fviz_pca_biplot(pca, repel = TRUE, col.var = "#2E9FDF", col.ind = "#696969")

#________________________________
# Covariance analyses between climate data and paleoenvironmental indicators

# /!\ D1 has no L6/L7
# /!\ dO18/dC13 has no L0,L1
# /!\ THI has no L0

# ISOTOPES

# Layers for Isotopes
L <- c("L2", "L3", "L4a", "L5", "L6", "L7", "L8")
# dO18 (Jeffrey 2016)
dO18_m <- c(-2.6, -3.7, -3.9, -4.3, -3.7, -5.5, -5.2); names(dO18_m) <- L
dO18_sd <- c(0.7, 0.8, 0.7, 0.7, 0.6, 0.8, 0.8); names(dO18_sd) <- L
# dC13 (Jeffrey 2016)
dC13_m <- c(-8.8, -8.7, -9.0, -9.4, -8.9, -9.5, -9.5); names(dC13_m) <- L
dC13_sd <- c(0.7, 0.4, 0.4, 0.6, 0.5, 0.8, 0.7); names(dC13_sd) <- L
# MAP (Jeffrey 2016)
MAP_min <- c(268, 360, 393, 442, 381, 612, 550); names(MAP_min) <- L
MAP_max <- c(405, 576, 594, 668, 543, 980, 897); names(MAP_min) <- L
MAP <- rowMeans(cbind(MAP_max, MAP_min))
# Concatenate
isotopes = cbind(dO18_m=dO18_m, dO18_max=dO18_m+dO18_sd, dO18_min=dO18_m-dO18_sd,
                 dC13_m=dC13_m, dC13_max=dC13_m+dC13_sd, dC13_min=dC13_m-dC13_sd,
                 MAP=MAP, MAP_max=MAP_max, MAP_min=MAP_min)
# Tests
# D1
I <- isotopes[-c(5,6),]
C <- pca$ind$coord[-c(1,2),2]; names(C) <- rownames(I) # PC1 or PC2 (change pc number)
C <- pca$ind$coord[-c(1,2),]; rownames(C) <- rownames(I) # All climate variables
# D2
I <- isotopes
C <- pca$ind$coord[-c(1,2),1]; names(C) <- rownames(I) # PC1 or PC2 (change pc number)
C <- pca$ind$coord[-c(1,2),]; rownames(C) <- rownames(I) # All climate variables
# Pairwise correlations
M <- cor(I, C) # R2
p.mat <- cor.mtest2(I, C) # p-values
cbind(M, p.mat) # alpha = 0.05
# 2B-pls
pls <- two.b.pls(I, C)
pls
plot(pls)


# THI

# Layers for THI
L <- c("L1", "L2", "L3", "L4a", "L5", "L6", "L7", "L8")
# THI % (Stoetzel 2011, Jeffrey 2016) / Oasis is not there because it is constant
forest <- c(22, 15, 17, 17, 15, 17, 13, 20); names(forest) <- L
bush <- c(24, 25, 25, 25, 22, 25, 22, 23); names(bush) <- L
Steppe <- c(31, 40, 33, 33, 38, 35, 40, 32); names(Steppe) <- L
wetland <- c(17, 10, 18, 18, 15, 16, 12, 15); names(wetland) <- L
rocky <- c(4, 8, 5, 5, 8, 5, 11, 8); names(rocky) <- L
# Concatenate
THI = cbind(forest=forest, bush=bush, Steppe=Steppe, wetland=wetland, rocky=rocky)
# Tests
# D1
I <- THI[-c(6,7),]
C <- pca$ind$coord[-1,1]; names(C) <- rownames(I) # PC1 or PC2 (change pc number)
C <- pca$ind$coord[-1,]; rownames(C) <- rownames(I) # All climate variables
# D2
I <- THI
C <- pca$ind$coord[-1,1]; names(C) <- rownames(I) # PC1 or PC2 (change pc number)
C <- pca$ind$coord[-1,]; rownames(C) <- rownames(I) # All climate variables
# Pairwise correlations
M <- cor(I, C) # R2
p.mat <- cor.mtest2(I, C) # p-values
cbind(M, p.mat) # alpha = 0.05
# 2B-pls
pls <- two.b.pls(I, C)
pls
plot(pls)


#________________________________________________________________
#________________________________________________________________

#
#### 2. Climate visualizations ####

#________________________________
# Climate variation over EH2 sequence

# D1
col <- colorRampPalette(c("white", "burlywood1", "sandybrown", "hotpink4"))(20) # D1
corrplot(scale(D1), is.corr = FALSE, tl.col = "black", tl.srt = 45, col=col, method="color")
# D2
col <- colorRampPalette(c("white", "paleturquoise", "lightseagreen", "paleturquoise4"))(20) # D2
corrplot(scale(D2), is.corr = FALSE, tl.col = "black", tl.srt = 45, col=col, method="color")
# => data scaled for visualization purpose (so that all m and sd are on the same scale)

#________________________________
# Dissimilarity maps

# Standardize crop of the area
cells

# Load crop data for Ctrl
path1
setwd(path1)
list <- list.files(pattern = "ymm.nc"); list
for (i in 1:length(list)){
  rast <- stack(list[i]) # data for var i
  rast <- rast[cells]; names(rast) <- cells # extract values for the crop
  mv(from="rast", to=paste("Ctrl", sep="_", var[i]))
}
# Load data for all periods
for (p in 1:length(P)){
  path <- get(paste("path",sep="",p))
  setwd(path)
  list <- list.files(pattern = "ymm.nc")
  for (i in 1:length(list)){
    rast <- stack(list[i]) # data for var i
    rast <- rast[46:47,70:71,] # extract values on EH2
    #rast <- rast[46,70,] # one cell
    mv(from="rast", to=paste(P[p], sep="_", var[i]))
  }
}

# Normalize each variable
for (v in 1:length(var)){
  var_data <- NULL
  for (i in 1:(length(P)+1)){
    var_data <- rbind(var_data, get(apropos(paste("_",sep="",var[v]))[i]))
  }
  var_data_norm <- (var_data-mean(na.omit(var_data)))/sd(na.omit(var_data))
  for (i in 1:length(P)){
    k=i*4-3;j=i*4
    x <- var_data_norm[k:j,]
    mv(from="x", to=apropos(paste("_",sep="",var[v]))[i])
  }
  x <- var_data_norm[29:402,]
  mv(from="x", to=apropos(paste("Ctrl_",sep="",var[v])))
}

# Mean and sd (past periods)
tab.P <- NULL
for (p in 1:length(P)){
  apropos(P[p])
  r <- NULL
  for (v in 1:length(var)){
    m <- mean(get(apropos(P[p])[v]))
    sd <- sd(get(apropos(P[p])[v]))
    r <- c(r, m, sd)
  }
  tab.P <- rbind(tab.P,r)
}
rownames(tab.P) <- P
colnames(tab.P) <- c("drysoil_frac_m", "drysoil_frac_sd",
                     "humtot_m", "humtot_sd", "hyd_stress_m", "hyd_stress_sd",
                     "precip_m","precip_sd", "qsurf_m", "qsurf_sd",
                     "sols_m", "sols_sd", "tsol_m", "tsol_sd",
                     "tz_ampl_day_m", "tz_ampl_day_sd", "w10m_m", "w10m_sd")
tab.P <- tab.P[-1,] # remove 00k

# Mean and sd (each cells of Ctrl)
tab.Ctrl <- NULL
for (c in 1:length(cells)){ # for each cell
  r <- NULL
  for (v in 1:length(var)){ # for each variable
    m <- mean(get(apropos("Ctrl")[v])[c,])
    sd <- sd(get(apropos("Ctrl")[v])[c,])
    r <- c(r, m, sd)
  }
  tab.Ctrl <- rbind(tab.Ctrl, r)
  print(paste(c,"/",length(cells)))
}
rownames(tab.Ctrl) <- paste("cell_",sep="",cells)
colnames(tab.Ctrl) <- c("drysoil_frac_m", "drysoil_frac_sd",
                        "humtot_m", "humtot_sd", "hyd_stress_m", "hyd_stress_sd",
                        "precip_m","precip_sd", "qsurf_m", "qsurf_sd",
                        "sols_m", "sols_sd", "tsol_m", "tsol_sd",
                        "tz_ampl_day_m", "tz_ampl_day_sd", "w10m_m", "w10m_sd")

# Distance vector
# for L1
tab.P; L <- 3 # Choose layer
rownames(tab.P)[L]
head(tab.Ctrl)
distance <- NULL
for (i in 1:nrow(tab.Ctrl)) {
  test_table <- rbind(tab.P[L,], tab.Ctrl[i,])
  dist_matrix <- stats::dist(test_table)
  distance <- c(distance, dist_matrix)
  print(paste(i,"/",nrow(tab.Ctrl)))
}
head(distance)
Sim <- (1-(distance/max(distance)))*100
# Similarity raster
rast <- stack(list[1])
rast <- rasterFromCells(rast, cells = cells) # raster of the size of the area of "cells"
values(rast) <- Sim # replacement of values by similarity index for the chosen layer
# Outline map
bmap=map_data("world", xlim=c(extent(rast)[1], extent(rast)[2]), ylim=c(extent(rast)[3], extent(rast)[4]))
# Plot map
x_df <- as.data.frame(rast, xy = TRUE)
# Plot
rast_dist <- ggplot() +
  geom_polygon(data=bmap, aes(x=long,y=lat,group=group),
               inherit.aes=F, colour='black', fill="black", lwd=0.35) +
  geom_raster(data = x_df, aes(x = x, y = y, fill = layer)) +
  coord_quickmap()+
  scale_fill_gradientn(colours = c(rep("transparent",7), c(rep("turquoise3",3))),
                       na.value = "transparent", limits=c(0,100))+
  theme_void()+
  theme(plot.background = element_rect(fill='white', color=NA))
rast_dist

#________________________________________________________________
#________________________________________________________________

