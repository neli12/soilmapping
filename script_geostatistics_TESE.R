#Carregar pacotes necesarios##
require(raster)
require(dplyr)
require(ggplot2)
library(gridExtra)
library(sp)
library(gstat)
library(rgdal)
require(hydroGOF)
require(dlookr)
require(corrplot)
library(e1071) 
library(rcompanion)
library(car)
library(farver)

list.files()
data.SYSIs <- read.csv("dados_p_semivariogram_color.csv", h=T,
                  sep=",")

data.singles <- read.csv("dados_preditos_cor_p_semivariogram.csv", h=T,
                       sep=",")
summary(data.SYSIs)
summary(data.singles )
coordinates(data.SYSIs) <- ~Long+Lat


coordinates(data.singles) <- ~Long+Lat


sand <- data.SYSIs$chroma
sandSYSI <- data.SYSIs$chromaSYSI
sandSEN <- data.SYSIs$chromaSEN
sandComb <- data.SYSIs$chromaCom
sandL8OLI <- data.singles$chromaL8.OLI.pred.data
sandS2MSI <- data.singles$chromaS2.MSI.pred.data

gSand <- gstat(id="sand", formula = chroma~1, data = data.SYSIs)
gSandSYSI <- gstat(id="sandSYSI", formula = chromaSYSI ~ 1, data = data.SYSIs)
gSandSEN <- gstat(id="sandSEN", formula = chromaSEN ~ 1, data = data.SYSIs)
gSandComb <- gstat(id="sandComb", formula = chromaCom ~ 1, data = data.SYSIs)
gSandL8OLI <- gstat(id="sandL8OLI", formula = chromaL8.OLI.pred.data ~ 1, data = data.singles)
gSandS2MSI <- gstat(id="sandS2MSI", formula = chromaS2.MSI.pred.data ~ 1, data = data.singles)

print(max(dist(data.SYSIs@coords/2)))
print(min(dist(data.SYSIs@coords)))

sand_exp <- gstat::variogram(gSand, width = 2800, cutoff = 30000)
sand_exp1 <- gstat::variogram(gSandComb, width = 2200, cutoff = 30000)
sand_exp2 <- gstat::variogram(gSandSYSI, width = 2200, cutoff = 30000)
sand_exp3 <- gstat::variogram(gSandSEN,width = 2200, cutoff = 30000)
sand_exp4 <- gstat::variogram(gSandL8OLI, width = 2200, cutoff = 30000)
sand_exp5 <- gstat::variogram(gSandS2MSI, width = 2200, cutoff = 30000)
plot(sand_exp)
plot(sand_exp1)
plot(sand_exp2)
plot(sand_exp3)
plot(sand_exp4)
plot(sand_exp5)

Sand.vgm <- fit.variogram(sand_exp, vgm(0.6, "Sph", 5000, 0.2))
SandComb.vgm <- fit.variogram(sand_exp1, vgm(1000, "", 5000, 2000))
SanSYSI.vgm <- fit.variogram(sand_exp2, vgm(1000, "Sph", 5000, 2000))
SandSEN.vgm <- fit.variogram(sand_exp3, vgm(1000, "Sph", 5000, 2000))
SandL8OLI.vgm <- fit.variogram(sand_exp4, vgm(0.8, "Exp", 5000, 0.2))
SandS2MSI.vgm <- fit.variogram(sand_exp5, vgm(0.8, "Exp", 5000, 0.2))


plot(sand_exp, Sand.vgm)
plot(sand_exp1, SandComb.vgm)
plot(sand_exp2, SanSYSI.vgm)
plot(sand_exp3, SandSEN.vgm)
plot(sand_exp4, SandL8OLI.vgm)
plot(sand_exp5, SandS2MSI.vgm)
SandL8OLI.vgm
SandS2MSI.vgm

{#Validação cruzada
##Experimental#
xvalid.exp <- krige.cv(sandSEN ~ 1, locations = data.SYSIs, model = vgm1)
plot(xvalid.exp$var1.pred ~ data.SYSIs$SandSYSI, cex = 1.2, lwd = 2, ylim=c(0,950))
abline(0, 1, col = "red", lwd = 2)
lm_exp <- lm(xvalid.exp$var1.pred ~ data.SYSIs$SandSYSI)
abline(lm_exp, col = "green", lwd = 2)
r2_exp <- summary(lm_exp)$r.squared
rmse_exp <- hydroGOF::rmse(xvalid.exp$var1.pred,  data.SYSIs$SandSYSI)
r2_exp
rmse_exp

###Spherical##
xvalid.sph <- krige.cv(sandSEN ~ 1, locations = data.SYSIs, model = vgm1.1)
plot(xvalid.sph$var1.pred ~ data.SYSIs$SandSYSI, cex = 1.2, lwd = 2)
abline(0, 1, col = "red", lwd = 2)
lm_sph <- lm(xvalid.sph$var1.pred ~  data.SYSIs$SandSYSI)
abline(lm_sph, col = "green", lwd = 2)
r2_sph <- summary(lm_sph)$r.squared
rmse_sph <- hydroGOF::rmse(xvalid.sph$var1.pred, data.SYSIs$SandSYSI)
r2_sph
rmse_sph

###Gaussian##
xvalid.gau <- krige.cv(SandSYSI ~ 1, locations = data.SYSIs, model = vgm1.2)
plot(xvalid.gau$var1.pred ~ data.SYSIs$SandSYSI, cex = 1.2, lwd = 2)
abline(0, 1, col = "red", lwd = 2)
lm_gau <- lm(xvalid.gau$var1.pred ~ data.SYSIs$SandSYSI)
abline(lm_gau, col = "green", lwd = 2)
r2_gau <- summary(lm_gau)$r.squared
rmse_gau <- hydroGOF::rmse(xvalid.gau$var1.pred, data.SYSIs$SandSYSI)
r2_gau
rmse_gau

}

vgLine <- rbind(#cbind(variogramLine(Sand.vgm, maxdist = max(sand_exp$dist)), id = "Measured"),
                #cbind(variogramLine(SandComb.vgm, maxdist = max(sand_exp3$dist)), id = "Predicted by SYSI Combined"),
                #cbind(variogramLine(SanSYSI.vgm, maxdist = max(sand_exp1$dist)), id = "Predicted by SYSI L8-OLI"),
                #cbind(variogramLine(SandSEN.vgm, maxdist = max(sand_exp2$dist)), id = "Predicted by SYSI S2-MSI"), 
                cbind(variogramLine(SandL8OLI.vgm, maxdist = max(sand_exp3$dist)), id = "Predicted by L8-OLI"),
                cbind(variogramLine(SandS2MSI.vgm, maxdist = max(sand_exp3$dist)), id = "Predicted by S2-MSI"))


p6 <- ggplot(data = sand_exp) + geom_point(aes(x=dist, y=gamma), col = "black", size=3) + 
  geom_point(data = sand_exp1, aes(x=dist, y=gamma), col = "blue", size=3) + 
  geom_point(data=sand_exp2, aes(x=dist, y=gamma), col = "green", size=3) + 
  geom_point(data=sand_exp3, aes(x=dist, y=gamma), col = "violet", size=3) + 
  geom_point(data=sand_exp4, aes(x=dist, y=gamma), col = "red", size=3) + 
  geom_point(data=sand_exp5, aes(x=dist, y=gamma), col = "orange", size=3) + 
  theme_bw() + #geom_line(data = vgLine, aes(x=dist, y=gamma, colour = id), size = 1.2) + 
  scale_color_manual(values = c("Measured" = "black", "Predicted by SYSI Combined" = "blue", "Predicted by SYSI L8-OLI" = "green",
                                "Predicted by SYSI S2-MSI" = "violet", "Predicted by L8-OLI" = "red", "Predicted by S2-MSI" = "orange")) + 
  ylab("") + xlab("") + xlim(0,30000) + theme(axis.text=element_text(size=22)) + 
  theme(axis.title=element_text(size=24)) + theme(legend.text = element_text(size = 22)) + theme(legend.position = c(0.7, 0.12)) + theme(legend.title = element_blank())
p6
dev.off()

require(ggpubr)
tiff("Semivariogram_all_DEF.tif", width = 4500, height = 6500, res = 300)
ggarrange(p1, p4, p2, p5, p3,p6, nrow=3, ncol=2, common.legend = TRUE, legend="bottom")

dev.off()

write.csv(train.color, "train_color.csv")
##Som##
tiff("Semivariogram_Som.tif", width = 2600, height = 2200, res = 300)

p2 <- ggplot(data = var_exp) + geom_point(aes(x=dist, y=gamma), col = "red") + 
  geom_point(data = var_exp1, aes(x=dist, y=gamma), col = "green") + 
  geom_point(data=var_exp2, aes(x=dist, y=gamma), col = "violet") + 
  geom_point(data=var_exp3, aes(x=dist, y=gamma), col = "blue") + theme_bw() + 
  geom_line(data = vgLine, aes(x=dist, y=gamma, colour = id), size = 1.2) + 
  scale_color_manual(values = c("Measured" = "red", "Predicted by SYSI L8-OLI" = "green", "Predicted by SYSI Fused" = "blue",
                                "Predicted by SYSI S2-MSI" = "violet")) + ylab("") + xlab("") + xlim(0,30000) + theme(axis.text=element_text(size=22)) + 
  theme(axis.title=element_text(size=24)) + theme(legend.text = element_text(size = 22)) + theme(legend.position = c(0.7, 0.12)) + theme(legend.title = element_blank())
p2
dev.off()

##Som##
tiff("Semivariogram_Som.tif", width = 2600, height = 2200, res = 300)

p3 <- ggplot(data = var_exp) + geom_point(aes(x=dist, y=gamma), col = "red") + 
  geom_point(data = var_exp1, aes(x=dist, y=gamma), col = "green") + 
  geom_point(data=var_exp2, aes(x=dist, y=gamma), col = "violet") + 
  geom_point(data=var_exp3, aes(x=dist, y=gamma), col = "blue") + theme_bw() + 
  geom_line(data = vgLine, aes(x=dist, y=gamma, colour = id), size = 1.2) + 
  scale_color_manual(values = c("Measured" = "red", "Predicted by SYSI L8-OLI" = "green", "Predicted by SYSI Fused" = "blue",
                                "Predicted by SYSI S2-MSI" = "violet")) + ylab("") + xlab("Distance (m)") + xlim(0,30000) + theme(axis.text=element_text(size=22)) + 
  theme(axis.title=element_text(size=24)) + theme(legend.text = element_text(size = 22)) + theme(legend.position = c(0.7, 0.5)) + theme(legend.title = element_blank())
p3
dev.off()

tiff("Semivariogram_hue.tif", width = 2600, height = 2200, res = 300)

p4 <- ggplot(data = var_exp) + geom_point(aes(x=dist, y=gamma), col = "red") + 
  geom_point(data = var_exp1, aes(x=dist, y=gamma), col = "green") + 
  geom_point(data=var_exp2, aes(x=dist, y=gamma), col = "violet") + 
  geom_point(data=var_exp3, aes(x=dist, y=gamma), col = "blue") + theme_bw() + 
  geom_line(data = vgLine, aes(x=dist, y=gamma, colour = id), size = 1.2) + ylim(0,7) + 
  scale_color_manual(values = c("Measured" = "red", "Predicted by SYSI L8-OLI" = "green", "Predicted by SYSI Fused" = "blue",
                                "Predicted by SYSI S2-MSI" = "violet")) + ylab("") + xlab("") + xlim(0,30000) + theme(axis.text=element_text(size=22)) + 
  theme(axis.title=element_text(size=24)) + theme(legend.text = element_text(size = 22)) + theme(legend.position = c(0.7, 0.12)) + theme(legend.title = element_blank())
p4
dev.off()

tiff("Semivariogram_value.tif", width = 2600, height = 2200, res = 300)

p5 <- ggplot(data = var_exp) + geom_point(aes(x=dist, y=gamma), col = "red") + 
  geom_point(data = var_exp1, aes(x=dist, y=gamma, colour=id), col = "green") + 
  geom_point(data=var_exp2, aes(x=dist, y=gamma), col = "violet") + 
  geom_point(data=var_exp3, aes(x=dist, y=gamma), col = "blue") + theme_bw()  + ylim(0,0.9) + 
  scale_color_manual(values = c("Measured" = "red", "Predicted by SYSI L8-OLI" = "green", "Predicted by SYSI Fused" = "blue",
                                "Predicted by SYSI S2-MSI" = "violet")) + ylab("") + xlab("") + xlim(0,30000) + theme(axis.text=element_text(size=22)) + 
  theme(axis.title=element_text(size=24)) + theme(legend.text = element_text(size = 22)) + theme(legend.position = c(0.7, 0.12)) + theme(legend.title = element_blank())
p5
dev.off()


tiff("Semivariogram_croma.tif", width = 2600, height = 2200, res = 300)

p6 <- ggplot(data = var_exp) + geom_point(aes(x=dist, y=gamma), col = "red") + 
  geom_point(data = var_exp1, aes(x=dist, y=gamma), col = "green") + 
  geom_point(data=var_exp2, aes(x=dist, y=gamma), col = "violet") + 
  geom_point(data=var_exp3, aes(x=dist, y=gamma), col = "blue") + theme_bw()  + ylim(0,0.8) + 
  scale_color_manual(values = c("Measured" = "red", "Predicted by SYSI L8-OLI" = "green", "Predicted by SYSI Fused" = "blue",
                                "Predicted by SYSI S2-MSI" = "violet")) + ylab("") + xlab("Distance (m)") + xlim(0,30000) + theme(axis.text=element_text(size=22)) + 
  theme(axis.title=element_text(size=24)) + theme(legend.text = element_text(size = 22)) + theme(legend.position = c(0.7, 0.12)) + theme(legend.title = element_blank())
p6
dev.off()
