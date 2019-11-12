# #################################################################### ####
# Title: Comparing the recovery of richness, structure and biomass in  ####
#        naturally regrowing and planted reforestation                 #### 
# Author: Timothy L Staples                                            ####
# Details: Minimal code to reproduce results and figures               ####  
# #################################################################### ####

rm(list=ls())
setwd("LOCATION OF DOWNLOADED REPOSITORY")

# Functions ####

sapply(list.files(path="./Functions", pattern=".R", full.names=TRUE),
       source)

# Colours for plotting
colours<-read.csv("./Data/plot.colours.csv")

# Libraries ####

# loads packages if installed, installs them if not
package.loader(c("lme4", "mgcv", "FD", "multcomp", "dendextend", "grImport",
                 "merTools", "plotrix"))

# PLOT COMPARISONS ####

#                               read in data ####

site <- read.csv("./Data/plot_attribute_data.csv", header=TRUE)
site$plot.type <- relevel(site$plot.type, "Remnant")


response.vars<-c("log.biomass.area",
                 "log.density",
                 "log.rare.rich")

#                               models ####

plot.type.models<-lapply(response.vars, function(x){
  
  temp.data<-site[!is.na(site[,x]),]
  
  lmer(as.formula(paste0(x, "~ plot.type + (1|preRE/planting)")), 
       data=temp.data)
})

plot.type.glht<-lapply(plot.type.models, function(x){
  glht(x, linfct = mcp(plot.type="Tukey"))
})

plot.glht.summ<-lapply(plot.type.glht, function(x){
  do.call("cbind", summary(x)$test[-(c(1,2,7))])
})
names(plot.glht.summ)<-response.vars

# which models show a difference between two groups?
plot.glht.sub<-plot.glht.summ[which(sapply(plot.glht.summ, function(x){
  ifelse(sum(x[,4] <=0.05) > 0, TRUE, FALSE)
}))]

glht.output<-do.call("rbind", lapply(1:length(plot.glht.sub), function(x){
  temp<-as.data.frame(plot.glht.sub[[x]])
  temp<-round(temp, 3)
  temp$var<-rep(names(plot.glht.sub)[x], dim(temp)[1])
  return(temp)
}))

write.csv(glht.output, "./Outputs/plot comparison glht.csv")

#                               Figure 1 ####

pdf("./Plots/Figure 1.pdf", height=4.5, width=1.9, useDingbats=FALSE)
par(mfrow=c(3,1), oma=c(2.5,3,1,0.5), mar=c(0,0,0,0), 
    las=1, tck=-0.02, ps=8, mgp=c(3,0,0))

names(plot.type.models)<-response.vars
temp.names<-cbind(response.vars,
                  c("Plot aboveground biomass ln(kg/ha",
                    "Plant density ln(plants/ha",
                    "Rarefied species richness (n = 12)"))

lapply(1:3, function(x){
  
  iless<-update(plot.type.models[[x]], .~. -1)
  raw.data<-plot.type.models[[x]]@frame
  head(raw.data)
  
  temp.coefs<-summary(iless)$coefficients
  levels<-substr(rownames(temp.coefs),
                 nchar("plot.type")+1, nchar(rownames(temp.coefs)))
  
  ylims<-rbind(log(c(1000,600000)),
               log(c(100,9000)),
               log(c(1,9)))[x,]
  
  plot(x=NULL, y=NULL, xlim=c(0.5,4.5), ylim=ylims,
       axes=FALSE, xlab="", ylab="")
  box()
  
  # Y-axes
  if(x==1){
    
    axis(side=2, at=log(c(1000,10000,100000,500000)), labels=c(1,10,100,500), mgp=c(3,0.5,0))
    axis(side=2, at=log(c(seq(1000,10000,1000),seq(10000,100000,10000), 
                          seq(100000,1000000,100000))), tck=-0.01,
         labels=NA, mgp=c(3,0.5,0))
  }
  
  if(x==2){
    
    axis(side=2, at=log(c(100,1000,5000)), labels=c(100,1000,5000), mgp=c(3,0.5,0))
    axis(side=2, at=log(c(seq(100,1000,100),seq(1000,10000,1000))), tck=-0.01,
         labels=NA, mgp=c(3,0.5,0))
  }
  
  if(x==3){
    
    axis(side=2, at=log(c(1,2,5)), labels=c(1,2,5), mgp=c(3,0.5,0))
    axis(side=2, at=log(seq(1,10,1)), tck=-0.01,
         labels=NA, mgp=c(3,0.5,0))
  }
  
  mtext(side=2, text= c("Plot aboveground biomass (Mg/ha)",
                        "Plant density (plants/ha)",
                        "Rarefied species richness (n = 12)")[x], line=c(1.75,2,1.5)[x], 
        las=0, cex=0.7)
  
  if(x==3){
    axis(side=1, at=1:4, 
         labels=c("Remnant", NA, NA, NA), cex=0.8)
    
    mtext(side=1, at=2, line=0, text="Planting", cex=0.6)
    
    mtext(side=1, at=3:4, line=0.6, 
          text=c("Young\nregrowth", "Old\nregrowth"), cex=0.6)
  } else {
    axis(side=1, at=1:4, 
         labels=NA, cex=0.8)
  }
  
  # raw points
  sapply(1:4, function(y){
    
    box.data<-plot.type.models[[x]]@frame[plot.type.models[[x]]@frame[,2]==levels[y],]
    
    #boxplot(box.data[,1], at=y, add=TRUE)
    points(y=box.data[,1], x=jitter(rep(y, dim(box.data)[1]), amount=0.2), 
           col=rgb(colours[y,5],
                   colours[y,6],
                   colours[y,7], 0.4),
           pch=16)
    
  })
  
  points(x=1:4, y=temp.coefs[,1], pch=16,
         col="black")
  
  segments(x0=1:4, x1=1:4,
           y0=temp.coefs[,1] + 1.96*temp.coefs[,2],
           y1=temp.coefs[,1] - 1.96*temp.coefs[,2], lwd=1.5,
           col="black")  
  
  # significance labels
  sig.labs <- rbind(c("A", "C", "B", "AB"),
                    c("A", "B", "A", "AB"),
                    c("A", "B", "B", "B"))[x,]
  
  text.adj <- c(0.35,0.2,0.1)[x]
  
  text(x=1:4, y=temp.coefs[,1] + 1.96*temp.coefs[,2] + text.adj,
       labels=sig.labs, font=2)
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.95, "y"), 
       labels=paste0("(",letters[1:3],")")[x], font=2,
       adj=0, cex=1.1)
  
})
dev.off()

# BIOMASS ACCRUAL MODEL ####
#                               data prep ####

model.data <- site

model.data$plot.type<-as.character(model.data$plot.type)
model.data$plot.type[model.data$plot.type %in% 
                        c("Young regrowth", "Old regrowth")] =  "Natural regrowth"
model.data$plot.type<-as.factor(model.data$plot.type)
model.data$plot.type<-relevel(model.data$plot.type, "Remnant")

model.data<-droplevels(model.data[model.data$plot.type !="Remnant", ])

#                               model ####

null.model<-lmer(log(biomass.area) ~ plot.type + (1|preRE/planting),
                 data=model.data, REML=FALSE)
summary(null.model)

age.model<-lmer(log(biomass.area) ~ log.age * plot.type + (1|preRE/planting),
                data=model.data,
                REML=FALSE)
summary(age.model)
anova(null.model, age.model)

age.quad.model<-lmer(log(biomass.area) ~ (log.age + I(log.age^2)) * plot.type + (1|preRE/planting),
                     data=model.data,
                     REML=FALSE,
                     control=lmerControl(optimizer="bobyqa"))
summary(age.quad.model)
anova(age.model, age.quad.model)

#                               Figure 2 ####

predict.points<-lapply(c("Planting","Natural regrowth"), function(x){
  
  temp.data<-model.data[model.data$plot.type==x,]
  factor.level<-unique(as.numeric(model.data$plot.type[model.data$plot.type == x]))
  
  new.data<-data.frame(log.age=seq(min(temp.data$log.age),
                                   max(temp.data$log.age), 
                                   length=200),
                       plot.type=factor(rep(x, 200),
                                        levels=levels(temp.data$plot.type)),
                       planting = "a",
                       preRE = "a")
  
  points <- predictInterval(age.quad.model, 
                            newdata = new.data,
                            n.sims = 3999,
                            which="fixed",
                            level=0.95,
                            include.resid.var=FALSE,
                            type="linear.prediction")
  
  return(cbind(new.data, points))
  
})

pdf("./Plots/Figure 2.pdf", height=2.2, width=2.5, useDingbats=FALSE)
par(mar=c(1.5,2,1,1), ps=6, tck=-0.01, mgp=c(3,0.2,0), las=1)
plot(x=NULL, y=NULL, 
     xlim=log(c(4, 150)), 
     ylim=c(log(1200), log(500000)), axes=FALSE)
box()
abline(v=log(85), lty="dashed")
axis(side=2, at=log(c(1000,5000,10000,50000,100000,500000)),
     labels=c(1,5,10,50,100,500))
axis(side=2, at=log(c(seq(1000,10000,1000),
                      seq(10000,100000,10000),
                      seq(100000,700000,100000))), tck=-0.005, labels=NA)
mtext(side=2, line=1, text="Plot aboveground biomass (Mg/ha)", las=0)

axis(side=1, at=log(c(1, 2, 5, 10, 20, 50)),
     labels=c(1, 2, 5, 10, 20, 50), mgp=c(3,-0.2,0))
axis(side=1, at=log(c(seq(1,20,1),
                      seq(20,80,10))), tck=-0.005, labels=NA)
axis(side=1, at=log(120), label="Remnant", mgp=c(3,-0.2,0))
mtext(side=1, line=0.5, text=expression("Regrowing forest age (years)"))

# Raw points
with(model.data,
     points(log(biomass.area) ~ jitter(log.age, amount=0.025),
            col=rgb(colours[c(4,2),5],
                    colours[c(4,2),6],
                    colours[c(4,2),7], 0.6)[plot.type], 
            cex=0.6, pch=c(17,15)[plot.type], lwd=0.5))

# separate regrowth polygon into 3 segments, based on where we have data and 
# where we don't. Not running at the moment
if(FALSE){
  reg.1<-max(which(predict.points[[2]]$log.age <= log(5)))
  reg.2<-max(which(predict.points[[2]]$log.age < log(10)))
  
  with(predict.points[[2]][1:reg.1,],
       polygon(x=c(log.age, rev(log.age)),
               y=c(fit + 1.96* se.fit,
                   rev(fit - 1.96* se.fit)),
               border=NA, col=rgb(0,0.5,0,0.5)))
  
  with(predict.points[[2]][reg.1:reg.2,],
       polygon(x=c(log.age, rev(log.age)),
               y=c(fit + 1.96* se.fit,
                   rev(fit - 1.96* se.fit)),
               border=NA, col=rgb(0,0.5,0,0.2)))
  
  with(predict.points[[2]][reg.2:200,],
       polygon(x=c(log.age, rev(log.age)),
               y=c(fit + 1.96* se.fit,
                   rev(fit - 1.96* se.fit)),
               border=NA, col=rgb(0,0.5,0,0.5)))
  
  with(predict.points[[2]][1:reg.1,],
       points(fit ~ log.age, type='l', col="darkgreen", lwd=2))
  
  with(predict.points[[2]][reg.1:reg.2,],
       points(fit ~ log.age, type='l', col="darkgreen", lwd=2, lty="22"))
  
  with(predict.points[[2]][reg.2:200,],
       points(fit ~ log.age, type='l', col="darkgreen", lwd=2))
}

# Remnant points
rem.points<-site[site$plot.type=="Remnant",]
with(rem.points, 
     points(log(biomass.area) ~ jitter(rep(log(120), dim(rem.points)[1]),
                                       amount=0.2),
            col=rgb(colours[1,5], 
                    colours[1,6], 
                    colours[1,7], 0.4), 
            cex=0.6, pch=16, lwd=0.5))

base.coefs<-summary(plot.type.models[[1]])$coefficients

points(y=base.coefs[1,1], x=log(120), pch=16, col=rgb(colours[1,2],
                                                      colours[1,3],
                                                      colours[1,4]), cex=0.8)
segments(y0=base.coefs[1,1] + 1.96* base.coefs[1,2],
         y1=base.coefs[1,1] - 1.96* base.coefs[1,2],
         x0=log(120), x1=log(120),
         col=rgb(colours[1,2],
                 colours[1,3],
                 colours[1,4]), lwd=1.5)

# CI Polygons
polygon(x=c(predict.points[[2]]$log.age, rev(predict.points[[2]]$log.age)),
        y=c(predict.points[[2]]$upr,
            rev(predict.points[[2]]$lwr)),
        border=NA, col=rgb(colours[4,2],
                           colours[4,3],
                           colours[4,4], 0.5))

polygon(x=c(predict.points[[1]]$log.age, rev(predict.points[[1]]$log.age)),
        y=c(predict.points[[1]]$upr,
            rev(predict.points[[1]]$lwr)),
        border=NA, col=rgb(colours[2,2],
                           colours[2,3],
                           colours[2,4], 0.5))

# Slopes
points(predict.points[[1]]$fit ~ predict.points[[1]]$log.age,
       type='l', lwd=1.5, col=rgb(colours[2,2],
                                  colours[2,3],
                                  colours[2,4]))

points(predict.points[[2]]$fit ~ predict.points[[2]]$log.age,
       type='l', lwd=1.5, col=rgb(colours[4,2],
                                  colours[4,3],
                                  colours[4,4]))

rect(xleft=log(25), 
     xright=log(85), 
     ybottom=par("usr")[3], 
     ytop=log(3500))

legend(x=log(85),
       y=relative.axis.point(-0.05, "y"),
       pch=c(15,17,16), legend=c("Planting", "Regrowth", "Remnant"),
       col=rgb(colours[c(2,4,1),2],
               colours[c(2,4,1),3],
               colours[c(2,4,1),4]),
       y.intersp=0.5, x.intersp=0.5, xjust=1, yjust=0,
       pt.cex=0.6, bty="n")

box()
dev.off()

summary(age.quad.model)
write.csv(round(summary(age.quad.model)$coefficients, 3), 
          "./Outputs/age quad model coefs.csv")
# SIZE CLASS DISTRIBUTION MODEL ####

#                               data prep ####

# cut by arbitrary number of points
cut.data <- read.csv("./Data/size_class_data.csv")

cut.data$Var1<-as.numeric(as.character(cut.data$Var1))

cut.data.bin<-cut.data
cut.data.bin$Freq <- ifelse(cut.data.bin$Freq > 0 , 1, 0)
cut.data.nozeros<-cut.data[cut.data$Freq > 0,]

bin.counts<-with(cut.data.bin[cut.data.bin$Freq==1,], as.data.frame(table(Var1, plot.type)))
bin.counts<-bin.counts[bin.counts$Freq==0,]
# model by mass and height as thin plate spline

# model will only converge if we truncate splines in each category so there are no
# 0 plot.type x mass bins
cut.data.bin<-cut.data.bin[!paste0(cut.data.bin$Var1, cut.data.bin$plot.type) %in%
                             paste0(bin.counts[,1], bin.counts[,2]),]

#                               models ####

cut.model.occur<-gamm(Freq ~ plot.type + s(Var1, bs="cr", k=6, by=plot.type),
                      random=list(preRE = ~1,
                                  planting=~1,
                                  pl.p=~1),
                      family=binomial,
                      data=cut.data.bin)

write.csv(cbind(summary(cut.model.occur$gam)$p.table,
                summary(cut.model.occur$gam)$s.table), "./Outputs/size class occurrence coefs.csv")

cut.model.count<-gamm(log(Freq) ~ plot.type + s(Var1, bs="cr", k=6, by=plot.type),
                      random=list(preRE = ~1,
                                  planting=~1,
                                  pl.p=~1),
                      family=gaussian,
                      data=cut.data.nozeros)

write.csv(cbind(summary(cut.model.count$gam)$p.table,
                summary(cut.model.count$gam)$s.table), "./Outputs/size class density coefs.csv")

#                               Figure 3 ####

count.gam.predict<-lapply(levels(cut.data.nozeros$plot.type), function(x){
  
  temp.data<-cut.data.nozeros[cut.data.nozeros$plot.type==x,]
  
  pred.data=data.frame(Var1=rep(seq(min(temp.data$Var1),
                                    max(temp.data$Var1),
                                    length=200)),
                       plot.type=rep(x, 200))
  
  temp.predict<-do.call("cbind", predict(cut.model.count$gam, 
                                         newdata=pred.data,
                                         se.fit=TRUE, type="response"))
  
  return(cbind(pred.data, temp.predict))
  
})
names(count.gam.predict)<-levels(cut.data$plot.type)

occur.gam.predict<-lapply(levels(cut.data$plot.type), function(x){
  
  temp.data<-cut.data.bin[cut.data.bin$plot.type==x,]
  
  pred.data=data.frame(Var1=rep(seq(min(temp.data$Var1),
                                    max(temp.data$Var1),
                                    length=200)),
                       plot.type=rep(x, 200))
  
  temp.predict<-do.call("cbind", predict(cut.model.occur$gam, 
                                         newdata=pred.data,
                                         se.fit=TRUE, type="link"))
  
  return(cbind(pred.data, temp.predict))
  
})
names(occur.gam.predict)<-levels(cut.data$plot.type)

pdf("./Plots/Figure 3.pdf", height=3.75, width=4, useDingbats=FALSE)
par(mfrow=c(2,2), mar=c(0,0,0,0), oma=c(2.5,3,1.5,1), las=1, 
    tck=-0.025, mgp=c(3,0.5,0), ps=8)

#                                     Occurrence subplots ####

ylims<-c(0,1)
xlims<-summary(cut.data$Var1)[c(1,6)]+c(-0.2,0.25)

plot(y=NULL, x=NULL, type="n",
     ylim=ylims, xlim=xlims, yaxs="i", xlab="", ylab="", axes=FALSE)

axis(side=1, at=log(c(0.01,0.1,1,10,100,1000)), 
     labels=NA, mgp=c(3,0.2,0))
box()

axis(side=1, at=log(c(seq(0.1,1,0.1),
                      seq(1,10,1), seq(10,100,10), seq(100,1000,100),
                      seq(1000,10000,1000))), labels=NA, tck=-0.01)

axis(side=2, at=seq(0,1,0.2), labels=seq(0,1,0.2))

mtext(side=2, line=1.5, text=expression("Probability of occurrence"), las=0, cex=0.8)

temp.predict<-occur.gam.predict[["Young regrowth"]]

temp.data<-droplevels(cut.data.bin[cut.data.bin$plot.type=="Young regrowth" &
                                     cut.data.bin$Var2 %in% temp.predict$Var2,])

points(jitter(temp.data$Freq, amount=0.05) ~ jitter(temp.data$Var1, amount=0.1),
       col=rgb(colours[3,2],
               colours[3,3],
               colours[3,4], 0.7),
       pch=c(1, 16)[as.factor(temp.data$Var2)])

# Young regrowth
with(temp.predict[1:200,],
     polygon(x=c(Var1, rev(Var1)),
             y=c(plogis(fit + 1.96*se.fit),
                 rev(plogis(fit - 1.96*se.fit))),
             col=rgb(colours[3,2],
                     colours[3,3],
                     colours[3,4], 0.5), border=NA))

with(temp.predict[1:200,],
     points(plogis(fit) ~ Var1, type="l",
            col=rgb(colours[3,2],
                    colours[3,3],
                    colours[3,4]), lwd=1.5))


# Plantings
temp.predict<-occur.gam.predict[["Planting"]]
temp.data<-droplevels(cut.data.bin[cut.data.bin$plot.type=="Planting" &
                                     cut.data.bin$Var2 %in% temp.predict$Var2,])

points(jitter(temp.data$Freq, amount=0.05) ~ jitter(temp.data$Var1, amount=0.1),
       col=rgb(colours[2,2],
               colours[2,3],
               colours[2,4], 0.7),
       pch=c(1, 16)[as.factor(temp.data$Var2)])

with(temp.predict[1:200,],
     polygon(x=c(Var1, rev(Var1)),
             y=c(plogis(fit + 1.96*se.fit),
                 rev(plogis(fit - 1.96*se.fit))),
             col=rgb(colours[2,2],
                     colours[2,3],
                     colours[2,4], 0.5), border=NA))

with(temp.predict[1:200,],
     points(plogis(fit) ~ Var1, type="l",
            col=rgb(colours[2,2],
                    colours[2,3],
                    colours[2,4]), lwd=1.5))


text(x=relative.axis.point(0.03, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(a)", font=2, cex=1, adj=0)

mtext(side=3, line=0.2, font=2, text="Young stands")

text(x=relative.axis.point(0.5, "x"), y=relative.axis.point(0.5, "y"),
     labels=c("Young regrowth", "Planting"), col=rgb(colours[3:2,2],
                                                     colours[3:2,3],
                                                     colours[3:2,4],1))


plot(y=NULL, x=NULL, type="n",
     ylim=ylims, xlim=xlims, yaxs="i", xlab="", ylab="", axes=FALSE)

axis(side=1, at=log(c(0.01,0.1,1,10,100,1000)), 
     labels=NA, mgp=c(3,0.2,0))
box()

axis(side=1, at=log(c(seq(0.1,1,0.1),
                      seq(1,10,1), seq(10,100,10), seq(100,1000,100),
                      seq(1000,10000,1000))), labels=NA, tck=-0.01)

axis(side=2, at=seq(0,1,0.2), labels=NA)

temp.predict<-occur.gam.predict[["Old regrowth"]]

temp.data<-droplevels(cut.data.bin[cut.data.bin$plot.type=="Young regrowth" &
                                     cut.data.bin$Var2 %in% temp.predict$Var2,])

points(jitter(temp.data$Freq, amount=0.05) ~ jitter(temp.data$Var1, amount=0.1),
       col=rgb(colours[4,2],
               colours[4,3],
               colours[4,4], 0.7),
       pch=c(1, 16)[as.factor(temp.data$Var2)])

# Planting
with(temp.predict[1:200,],
     polygon(x=c(Var1, rev(Var1)),
             y=c(plogis(fit + 1.96*se.fit),
                 rev(plogis(fit - 1.96*se.fit))),
             col=rgb(colours[4,2],
                     colours[4,3],
                     colours[4,4], 0.5), border=NA))

with(temp.predict[1:200,],
     points(plogis(fit) ~ Var1, type="l",
            col=rgb(colours[4,2],
                    colours[4,3],
                    colours[4,4]), lwd=1.5))

temp.predict<-occur.gam.predict[["Remnant"]]

# Young regrowth
with(temp.predict[1:200,],
     polygon(x=c(Var1, rev(Var1)),
             y=c(plogis(fit + 1.96*se.fit),
                 rev(plogis(fit - 1.96*se.fit))),
             col=rgb(colours[1,2],
                     colours[1,3],
                     colours[1,4], 0.5), border=NA))

with(temp.predict[1:200,],
     points(plogis(fit) ~ Var1, type="l",
            col=rgb(colours[1,2],
                    colours[1,3],
                    colours[1,4]), lwd=1.5))

text(x=relative.axis.point(0.03, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(b)", font=2, cex=1, adj=0)

mtext(side=3, line=0.2, font=2, text="Old stands")

text(x=relative.axis.point(0.5, "x"), y=relative.axis.point(0.5, "y"),
     labels=c("Old regrowth", "Remnant"), col=rgb(colours[c(4,1),2],
                                                  colours[c(4,1),3],
                                                  colours[c(4,1),4],1))

#                                     Abundance subplots ####

ylims<-log(c(4,400))
xlims<-summary(cut.data$Var1)[c(1,6)]+c(-0.2,0.25)

plot(y=NULL, x=NULL, type="n",
     ylim=ylims, xlim=xlims, yaxs="i", xlab="", ylab="", axes=FALSE)

axis(side=1, at=log(c(0.01,0.1,1,10,100,1000)), 
     labels=c(0.01,0.1,1,10,100,1000), mgp=c(3,0.2,0))
box()

axis(side=1, at=log(c(seq(0.1,1,0.1),
                      seq(1,10,1), seq(10,100,10), seq(100,1000,100),
                      seq(1000,10000,1000))), labels=NA, tck=-0.01)
mtext(side=1, line=1.2, at=par("usr")[2], text="Plant aboveground biomass (kg)", cex=0.8)

axis(side=2, at=log(c(0,1,10,100,800)+1), labels=c(0,1,10,100,800))
axis(side=2, at=log(c(seq(1,10,1), seq(10,100,10), seq(100,800,100))+1), labels=NA, tck=-0.01)

mtext(side=2, line=1.5, text="Plant density (plants/ha)", las=0, cex=0.8)

temp.predict<-count.gam.predict[["Planting"]]

with(temp.predict[1:200,],
     polygon(x=c(Var1, rev(Var1)),
             y=c(fit + 1.96*se.fit,
                 rev(fit - 1.96*se.fit)),
             col=rgb(colours[2,2],
                     colours[2,3],
                     colours[2,4], 0.5), border=NA))

with(temp.predict[1:200,],
     points(fit ~ Var1, type="l",
            col=rgb(colours[2,2],
                    colours[2,3],
                    colours[2,4]), lwd=1.5))

with(temp.predict[201:400,],
     polygon(x=c(Var1, rev(Var1)),
             y=c(fit + 1.96*se.fit,
                 rev(fit - 1.96*se.fit)),
             border=rgb(colours[2,2],
                        colours[2,3],
                        colours[2,4], 1), col=rgb(1,1,1,0)), lty="dashed")

with(temp.predict[201:400,],
     points(fit ~ Var1, type="l",
            col=rgb(colours[2,2],
                    colours[2,3],
                    colours[2,4]), lwd=1.5, lty="dashed"))

temp.predict<-count.gam.predict[["Young regrowth"]]

with(temp.predict[1:200,],
     polygon(x=c(Var1, rev(Var1)),
             y=c(fit + 1.96*se.fit,
                 rev(fit - 1.96*se.fit)),
             col=rgb(colours[3,2],
                     colours[3,3],
                     colours[3,4], 0.5), border=NA))

with(temp.predict[1:200,],
     points(fit ~ Var1, type="l",
            col=rgb(colours[3,2],
                    colours[3,3],
                    colours[3,4]), lwd=1.5))

text(x=relative.axis.point(0.03, "x"),
     y=relative.axis.point(0.92, "y"),
     labels="(c)", font=2, cex=1, adj=0)

text(x=relative.axis.point(0.5, "x"), y=relative.axis.point(0.5, "y"),
     labels=c("Young regrowth", "Planting"), col=rgb(colours[3:2,2],
                                                     colours[3:2,3],
                                                     colours[3:2,4],1))

plot(y=NULL, x=NULL, type="n",
     ylim=ylims, xlim=xlims, yaxs="i", xlab="", ylab="", axes=FALSE)

axis(side=1, at=log(c(0.01,0.1,1,10,100,1000)), 
     labels=c(0.01,0.1,1,10,100,1000), mgp=c(3,0.2,0))
box()

axis(side=1, at=log(c(seq(0.1,1,0.1),
                      seq(1,10,1), seq(10,100,10), seq(100,1000,100),
                      seq(1000,10000,1000))), labels=NA, tck=-0.01)

temp.predict<-count.gam.predict[["Old regrowth"]]

with(temp.predict[1:200,],
     polygon(x=c(Var1, rev(Var1)),
             y=c(fit + 1.96*se.fit,
                 rev(fit - 1.96*se.fit)),
             col=rgb(colours[4,2],
                     colours[4,3],
                     colours[4,4], 0.5), border=NA))

with(temp.predict[1:200,],
     points(fit ~ Var1, type="l",
            col=rgb(colours[4,2],
                    colours[4,3],
                    colours[4,4]), lwd=1.5))

temp.predict<-count.gam.predict[["Remnant"]]

with(temp.predict[1:200,],
     polygon(x=c(Var1, rev(Var1)),
             y=c(fit + 1.96*se.fit,
                 rev(fit - 1.96*se.fit)),
             col=rgb(colours[1,2],
                     colours[1,3],
                     colours[1,4], 0.5), border=NA))

with(temp.predict[1:200,],
     points(fit ~ Var1, type="l",
            col=rgb(colours[1,2],
                    colours[1,3],
                    colours[1,4]), lwd=1.5))

text(x=relative.axis.point(0.03, "x"),
     y=relative.axis.point(0.92, "y"),
     labels="(d)", font=2, cex=1, adj=0)

text(x=relative.axis.point(0.5, "x"), y=relative.axis.point(0.5, "y"),
     labels=c("Old regrowth", "Remnant"), col=rgb(colours[c(4,1),2],
                                                  colours[c(4,1),3],
                                                  colours[c(4,1),4],1))

dev.off()

# SMALL PLANTS IN PLANTINGS PLOT ####

small.groups <- readRDS("./Data/small_planting_groups_mat.RDS")
small.reg.groups <- readRDS("./Data/small_regrowth_groups_mat.RDS")

pdf("./Plots/Figure 4.pdf", height=3, width=6)

par(mfrow=c(1,2), mar=c(0,0,0,0), oma=c(3,3.5,1,1), ps=8, tck=-0.025, mgp=c(3,0.5,0),
    las=1)

plot(x=NULL, y=NULL, xlim=c(0.5,8.5), ylim=c(0,1000), axes=FALSE,
     xlab="", ylab="")
axis(side=1, at=1:8, mgp=c(3,0.2,0))
axis(side=2)
mtext(side=1, at=par("usr")[2], line=1.25, text="Growth functional group")
mtext(side=2, line=2, text="Count of small plants (< 40 kg)", las=0)
box()

sapply(colnames(small.groups), function(x){
  
  temp.sort<-sort(small.groups[,x])
  temp.sort<-temp.sort[temp.sort > 0]
  
  temp.mat<-c(0, temp.sort)
  
  for(i in 1:length(temp.mat)){
    
    rect(ybottom=ifelse(i>1, sum(temp.mat[1:(i-1)]), 0),
         ytop=sum(temp.mat[1:i]),
         xleft=as.numeric(x)-0.3,
         xright=as.numeric(x)+0.3,
         col=c("grey80", "grey50")[(i %% 2) + 1])
    
  }
  
  top.sort<-temp.sort>=20
  
  if(sum(top.sort)>0){
    text(labels=names(temp.sort[top.sort]),
         x=as.numeric(x),
         y=sapply(names(temp.sort[top.sort]), function(y){
           sum(temp.sort[1:which(names(temp.sort)==y)])}) - 0.5*temp.sort[top.sort],
         pos=2, font=3, cex=0.8)
    
    segments(x0=as.numeric(x), x1=as.numeric(x)-0.3,
             y0=sapply(names(temp.sort[top.sort]), function(y){
               sum(temp.sort[1:which(names(temp.sort)==y)])}) - 0.5*temp.sort[top.sort],
             y1=sapply(names(temp.sort[top.sort]), function(y){
               sum(temp.sort[1:which(names(temp.sort)==y)])}) - 0.5*temp.sort[top.sort])
  }
  
})

segments(x0=5-0.3, x1=5+0.3, y0=0, y1=0)
text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(a) Plantings", font=2, adj=0)

# ANd compared to young regrowth....

plot(x=NULL, y=NULL, xlim=c(0.5,8.5), ylim=c(0,1000), axes=FALSE)
axis(side=1, at=1:8, mgp=c(3,0.2,0))
axis(side=2, labels=NA)
box()

sapply(colnames(small.reg.groups), function(x){
  
  print(x)
  
  temp.sort<-sort(small.reg.groups[,x])
  temp.sort<-temp.sort[temp.sort > 0]
  
  temp.mat<-c(0, temp.sort)
  
  for(i in 1:length(temp.mat)){
    
    rect(ybottom=ifelse(i>1, sum(temp.mat[1:(i-1)]), 0),
         ytop=sum(temp.mat[1:i]),
         xleft=as.numeric(x)-0.3,
         xright=as.numeric(x)+0.3,
         col=c("grey80", "grey50")[(i %% 2) + 1])
    
  }
  
  top.sort<-temp.sort>=20
  
  if(sum(top.sort)>0){
    text(labels=names(temp.sort[top.sort]),
         x=as.numeric(x),
         y=sapply(names(temp.sort[top.sort]), function(y){
           sum(temp.sort[1:which(names(temp.sort)==y)])}) - 0.5*temp.sort[top.sort],
         pos=2, font=3, cex=0.8)
    
    segments(x0=as.numeric(x), x1=as.numeric(x)-0.3,
             y0=sapply(names(temp.sort[top.sort]), function(y){
               sum(temp.sort[1:which(names(temp.sort)==y)])}) - 0.5*temp.sort[top.sort],
             y1=sapply(names(temp.sort[top.sort]), function(y){
               sum(temp.sort[1:which(names(temp.sort)==y)])}) - 0.5*temp.sort[top.sort])
  }
  
})

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(b) Young regrowth", font=2, adj=0)

dev.off()
# FUNCTIONAL GROUP RICHNESS MODELS ####

#                               seed functional groups ####

seed.count <- read.csv("./Data/seed_fungroup_count.csv")
seed.mm <- readRDS("./Data/seed_model_matrix.RDS")
seed.dend <- readRDS("./Data/seed_dendrogram.RDS")
seed.SM <- readRDS("./Data/seed_mass_group_means.RDS")
seed.class <- readRDS("./Data/seed_disp_class.RDS")
seed.removal <- readRDS("./Data/seed_missing_groups.RDS")
seed.split <- readRDS("./Data/seed_group_membership.RDS")

seed.model<-glmmPQL(count ~ 0 + seed.mm, random=list(~1|preRE,
                                                     ~1|planting,
                                                     ~1|pl.p),
                    family=negative.binomial(theta = 1),
                    data=seed.count)

write.csv(summary(seed.model)$tTable, "./Outputs/seed group model coefs.csv")

#                               growth functional groups ####

growth.count <- read.csv("./Data/growth_fungroup_count.csv")
growth.mm <- readRDS("./Data/growth_model_matrix.RDS")
growth.dend <- readRDS("./Data/growth_dendrogram.RDS")
growth.SLA <- readRDS("./Data/SLA_group_means.RDS")
growth.WD<- readRDS("./Data/wood_dens_group_means.RDS")
growth.MH <- readRDS("./Data/max_height_group_means.RDS")
growth.split <- readRDS("./Data/growth_group_membership.RDS")
n.groups <- readRDS("./Data/n_fixer.RDS")

growth.model<-glmmPQL(count ~ -1 + growth.mm, random=list(~1|preRE,
                                                          ~1|planting,
                                                          ~1|pl.p),
                      family=negative.binomial(theta = 1),
                      data=growth.count)

write.csv(summary(growth.model)$tTable, "./Outputs/growth group model coefs.csv")
#                               Figure 5 ####
#                                     seed data prep ####

# rerun seed models to get reference level coefs and SEs
seed.coef.df<-as.data.frame(summary(seed.model)$tTable)

seed.coef.df$plot.type<-substr(rownames(seed.coef.df),
                               regexpr("R|P|O|Y", rownames(seed.coef.df)),
                               nchar(rownames(seed.coef.df)))
seed.coef.df$plot.type[!seed.coef.df$plot.type %in% c("Planting", 
                                                      "Young regrowth", 
                                                      "Old regrowth")] = "Remnant"
seed.coef.df$group<-substr(rownames(seed.coef.df),
                           regexpr("seed\\.group", rownames(seed.coef.df))+10,
                           regexpr("seed\\.group", rownames(seed.coef.df))+10)
seed.coef.df$group[is.na(as.numeric(seed.coef.df$group))]=1

# calculate and write comparison to remnants 
seed.cut.mat<-model.matrix(seed.model, data=seed.count)
colnames(seed.cut.mat)<-gsub("seed.mm", "", colnames(seed.cut.mat))
colnames(seed.cut.mat)
seed.cut.mat[,"plot.area"]=0
seed.cut.mat<-unique.matrix(seed.cut.mat, MARGIN=1)
seed.marg.mat<-seed.cut.mat # save this to predict actual coefficients

# create contrast matrix to compare group estimates to remnants
seed.cont.mat<-seed.cut.mat
seed.cont.mat[,1:8]<-0
seed.cont.mat<-seed.cont.mat[rowSums(seed.cont.mat)>0,]

# get group names and plot types for each row of marginal matrix
ind<-which(seed.cont.mat !=0, arr.ind=TRUE)
ind<-ind[order(ind[,1]),]

colnames(seed.cont.mat)
plot.type.ind<-ind[ind[,2] %in% 9:11,]
group.ind<-ind[!ind[,2] %in% 9:11,]

colnames(seed.cont.mat)
plot.types<-gsub("plot.type", "", colnames(seed.cont.mat)[plot.type.ind[,2]])

groups<-rep(1, length(plot.types))
groups[group.ind[,1]]<-substr(colnames(seed.cont.mat)[group.ind[,2]],
                              nchar("seed.group")+1,
                              regexpr(":", colnames(seed.cont.mat)[group.ind[,2]])-1)

seed.marg.groups<-cbind(groups, plot.types)
seed.keeps<-!duplicated(paste0(seed.marg.groups[,1],
                               seed.marg.groups[,2]))

seed.marg.groups<-seed.marg.groups[seed.keeps,]
seed.cont.mat<-seed.cont.mat[seed.keeps,]

seed.contrast<-summary(glht(seed.model, linfct=seed.cont.mat))

seed.contrast.df<-data.frame(group=seed.marg.groups[,1],
                             plot.type=seed.marg.groups[,2],
                             estimate=seed.contrast$test$coefficients,
                             se=seed.contrast$test$sigma,
                             t=seed.contrast$test$tstat,
                             p=seed.contrast$test$pvalues)

seed.contrast.df[,-(1:2)] = round(seed.contrast.df[,-(1:2)], 3)
write.csv(seed.contrast.df, "./Outputs/seed contrast test.csv")

# get plot coefficients and error

seed.marg.df<-data.frame(seed.mm=I(seed.marg.mat))
gsub("seed\\.mm\\.", "seed\\.mm", colnames(seed.marg.df))

seed.rem.df<-data.frame(group=substr(seed.removal, 1,1),
                        plot.type=substr(seed.removal, 12, nchar(seed.removal)),
                        stringsAsFactors = FALSE)

seed.rem.rows<-apply(seed.rem.df, 1, function(x){
  
  #group column
  if(x[1]==1){
    gr.col<-which(rowSums(seed.marg.mat[,2:8])==0)
  } else {gr.col<-which(seed.marg.mat[,2:8][,grepl(x[1], 
                                                   colnames(seed.marg.mat[,2:8]))] == 1)}
  
  #plot type column
  pl.col<-which(seed.marg.mat[,9:11][,grepl(x[2], colnames(seed.marg.mat[,9:11]))] == 1)
  
  return(gr.col[gr.col %in% pl.col])
  
})

seed.marg.mat<-seed.marg.mat[-seed.rem.rows,]

seed.preds<-data.frame(estimate=seed.marg.mat%*%fixef(seed.model),
                       se=sqrt(diag(seed.marg.mat %*% tcrossprod(vcov(seed.model),seed.marg.mat))))
seed.preds$plot.type= "Remnant" 

plot.type.ind<-which(seed.marg.mat[,9:11] == 1, arr.ind=TRUE)
seed.preds$plot.type[plot.type.ind[,1]] =
  c("Planting", "Young regrowth", "Old regrowth")[plot.type.ind[,2]]

seed.preds$group<-1
colnames(seed.marg.mat)
group.ind<-which(seed.marg.mat[,2:8] == 1, arr.ind=TRUE)
seed.preds$group[group.ind[,1]] = (2:8)[group.ind[,2]]

#                                     growth data prep ####
growth.coef.df<-as.data.frame(summary(growth.model)$tTable)

growth.coef.df$plot.type<-substr(rownames(growth.coef.df),
                                 regexpr("R|P|O|Y", rownames(growth.coef.df)),
                                 nchar(rownames(growth.coef.df)))
growth.coef.df$plot.type[!growth.coef.df$plot.type %in% c("Planting", 
                                                          "Young regrowth", 
                                                          "Old regrowth")] = "Remnant"

growth.coef.df$group<-substr(rownames(growth.coef.df),
                             regexpr("growth\\.group", rownames(growth.coef.df))+12,
                             regexpr("growth\\.group", rownames(growth.coef.df))+12)
growth.coef.df$group[is.na(as.numeric(growth.coef.df$group))]=1

# calculate and write comparison to remnants 
growth.cut.mat<-model.matrix(growth.model, data=growth.count)
colnames(growth.cut.mat)<-gsub("growth.mm", "", colnames(growth.cut.mat))
growth.cut.mat[,"plot.area"]=0
growth.cut.mat<-unique.matrix(growth.cut.mat, MARGIN=1)
growth.marg.mat<-growth.cut.mat # save this to predict actual coefficients

# create contrast matrix to compare group estimates to remnants
# rows
growth.cont.mat<-growth.cut.mat
growth.cont.mat[,1:8]<-0
growth.cont.mat<-growth.cont.mat[rowSums(growth.cont.mat)>0,]

# get group names and plot types for each row of marginal matrix
ind<-which(growth.cont.mat !=0, arr.ind=TRUE)
ind<-ind[order(ind[,1]),]
plot.type.ind<-ind[ind[,2] %in% 9:11,]
group.ind<-ind[!ind[,2] %in% 9:11,]

plot.types<-gsub("plot.type", "", colnames(growth.cont.mat)[plot.type.ind[,2]])

groups<-rep(1, length(plot.types))
groups[group.ind[,1]]<-substr(colnames(growth.cont.mat)[group.ind[,2]],
                              nchar("growth.group")+1,
                              regexpr(":", colnames(growth.cont.mat)[group.ind[,2]])-1)

growth.marg.groups<-cbind(groups, plot.types)
growth.keeps<-!duplicated(paste0(growth.marg.groups[,1],
                                 growth.marg.groups[,2]))

growth.marg.groups<-growth.marg.groups[growth.keeps,]
growth.cont.mat<-growth.cont.mat[growth.keeps,]

growth.contrast<-summary(glht(growth.model, linfct=growth.cont.mat))

growth.contrast.df<-data.frame(group=growth.marg.groups[,1],
                               plot.type=growth.marg.groups[,2],
                               estimate=growth.contrast$test$coefficients,
                               se=growth.contrast$test$sigma,
                               t=growth.contrast$test$tstat,
                               p=growth.contrast$test$pvalues)

growth.contrast.df[,-(1:2)] = round(growth.contrast.df[,-(1:2)], 3)
write.csv(growth.contrast.df, "./Outputs/growth contrast test.csv")

# get plot coefficients and error

growth.rem.df<-data.frame(group=substr(growth.removal, 1,1),
                          plot.type=substr(growth.removal, 12, nchar(growth.removal)),
                          stringsAsFactors = FALSE)

growth.rem.rows<-apply(growth.rem.df, 1, function(x){
  
  #group column
  if(x[1]==1){
    gr.col<-which(rowSums(growth.marg.mat[,2:8])==0)
  } else {gr.col<-which(growth.marg.mat[,2:8][,grepl(x[1], 
                                                     colnames(growth.marg.mat[,2:8]))] == 1)}
  
  #plot type column
  pl.col<-which(growth.marg.mat[,9:11][,grepl(x[2], colnames(growth.marg.mat[,9:11]))] == 1)
  
  return(gr.col[gr.col %in% pl.col])
  
})

growth.marg.mat<-growth.marg.mat[-growth.rem.rows,]

growth.preds<-data.frame(estimate=growth.marg.mat%*%fixef(growth.model),
                         se=sqrt(diag(growth.marg.mat %*% tcrossprod(vcov(growth.model), 
                                                                     growth.marg.mat))))

growth.preds$plot.type= "Remnant" 
colnames(growth.marg.mat)
plot.type.ind<-which(growth.marg.mat[,9:11] == 1, arr.ind=TRUE)
growth.preds$plot.type[plot.type.ind[,1]] =
  c("Planting", "Young regrowth", "Old regrowth")[plot.type.ind[,2]]

growth.preds$group<-1
group.ind<-which(growth.marg.mat[,2:8] == 1, arr.ind=TRUE)
growth.preds$group[group.ind[,1]] = (2:8)[group.ind[,2]]

#                                     plot setup ####
pdf("./Plots/Figure 5.pdf", height=9, width=3.74, useDingbats=FALSE)
#

split.screen(rbind(c(0.05,0.65,0.65,0.98),
                   c(0.5,0.99,0.65,0.98),
                   
                   c(0.05,0.65, 0.16, 0.57),
                   c(0.5,0.99,0.16,0.57),
                   
                   c(0.05,0.99,0.65,0.98),
                   c(0.05,0.99,0.16,0.57),
                   
                   c(0.05,0.99,0.01,0.13)))

#                                     seed dendro ####
screen(1)
par(mar=c(0,0,0,5), oma=c(0,0,0,0), ps=8, tck=-0.02, mgp=c(3,0.3,0))

plot(dendrapply(as.dendrogram(seed.dend), 
                function(x) { attr(x, "height") <- log(attr(x, "height")+1); x }), 
     horiz=TRUE, yaxs="i", leaflab="none", axes=FALSE, ylab="", xlab="",
     ylim=c(-1,68.50), xlim=c(3,0), lwd=0.5, hang=1)

dend.pars<-par("usr")
dend.nodes<-as.data.frame(cbind(get_nodes_xy(seed.dend, type = c("rectangle"), center = FALSE,
                                             horiz = TRUE),
                                get_nodes_attr(seed.dend, attribute="label")))

dend.nodes$V1<-as.numeric(as.character(dend.nodes$V1))
dend.nodes<-dend.nodes[!is.na(dend.nodes$V3),]
dend.nodes$group<-seed.split[match(dend.nodes$V3, names(seed.split))]
dend.groups<-do.call("rbind", lapply(split(dend.nodes, f=dend.nodes$group), 
                                     function(x){c(min(x$V1),max(x$V1))}))
group.centers<-tapply(dend.nodes$V1, dend.nodes$group, mean)

text(x=relative.axis.point(0.08, "x"), 
     y=sort(group.centers), labels=8:1, cex=1.2, font=2,
     col=c("grey80", "grey50"))

seed.group.decode<-cbind(names(group.centers[order(group.centers)]),
                         8:1)
write.csv(seed.group.decode, "./Outputs/seed.groups.decoded.csv")

par(xpd=NA)
rect(xleft=par("usr")[1], xright=relative.axis.point(2.76, "x"),
     ybottom=dend.groups[names(sort(group.centers))[c(2,4,6,8)],1]-0.5, 
     ytop=dend.groups[names(sort(group.centers))[c(2,4,6,8)],2]+0.5,
     col=rgb(0,0,0,0.1), border=NA)

rect(xleft=-0.2, xright=-0.3, ybottom=group.centers-0.5, ytop=group.centers+2.5, lty="11")
rect(xleft=-0.2, xright=-0.3, ybottom=group.centers-0.5, ytop=(group.centers-0.5)+seed.SM, 
     col="black")
segments(x0=-0.15, x1=-0.35, y0=group.centers-0.5, y1=group.centers-0.5)
text(x=-0.25, y=group.centers-1.5, labels="M", adj=0.5, cex=0.8)
text(x=-0.6, y=group.centers-1.5, labels=seed.class, adj=0.5, cex=0.8)

picture.ys<-0.555 + ((par("usr")[3] + group.centers)/par("usr")[4])*(0.98-0.555)

PostScriptTrace("./Data/wind.ps")
wind<-readPicture("wind.ps.xml")
PostScriptTrace("./Data/fleshy.ps")
fleshy<-readPicture("fleshy.ps.xml")
PostScriptTrace("./Data/elaisome.ps")
elaisome<-readPicture("elaisome.ps.xml")
PostScriptTrace("./Data/unassisted.ps")
unassisted<-readPicture("unassisted.ps.xml")

lapply(picture.ys[seed.class=="Wi"], function(x){
  grid.picture(wind, 
               x=0.435, 
               y=x+0.03, 
               width=0.03, height=0.03)
})

lapply(picture.ys[seed.class=="Fl"], function(x){
  grid.picture(fleshy, 
               x=0.435, 
               y=x+0.055, 
               width=0.04, height=0.04)
})

lapply(picture.ys[seed.class=="El"], function(x){
  grid.picture(elaisome, 
               x=0.435, 
               y=x+0.055, 
               width=0.04, height=0.04)
})

lapply(picture.ys[seed.class=="Un"], function(x){
  grid.picture(unassisted, 
               x=0.435, 
               y=x+0.08, 
               width=0.03, height=0.03)
})

par(xpd=FALSE)

axis(side=1, mgp=c(3,0,0))
dendro.pars<-par("usr")
mtext(side=1, line=0.75, text="ln(Ward's minimum variance)")
close.screen(1)

#                                     seed points ####
screen(2)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), ps=8, tck=-0.02, mgp=c(3,0.3,0))
plot(x=NULL, y=NULL, xlim=c(-5.6, 1.25), ylim=dend.pars[3:4], 
     yaxs="i", axes=FALSE, ylab="", xlab="")

axis(side=1, at=log(c(0.001,0.01,0.1,1,2,3)), labels=c(0.001,0.01,0.1,1,2,3),
     mgp=c(3,0,0))
axis(side=1, at=log(c(seq(0.001,0.01,0.001), seq(0.01,0.1,0.01), seq(0.1,1,0.1),
                      1,3,1)), labels=NA, tck=-0.01)
mtext(side=1, line=0.75, text="Mean species per plot")

plot.type.pos<-group.centers[seed.preds$group] + 
  c(-2.25, -0.75, 0.75, 2.25)[as.factor(seed.preds$plot.type)]

segments(x0=seed.preds$estimate + 1.96*seed.preds$se,
         x1=seed.preds$estimate - 1.96*seed.preds$se,
         y0=plot.type.pos, 
         y1=plot.type.pos, lwd=1.5,
         col=rgb(colours[match(seed.preds$plot.type, colours[,1]),2],
                 colours[match(seed.preds$plot.type, colours[,1]),3],
                 colours[match(seed.preds$plot.type, colours[,1]),4]))

points(x=seed.preds$estimate,
       y=plot.type.pos,
       pch=c(25,22,21,24)[as.factor(seed.preds$plot.type)],
       cex=0.8,
       bg=rgb(colours[match(seed.preds$plot.type, colours[,1]),2],
              colours[match(seed.preds$plot.type, colours[,1]),3],
              colours[match(seed.preds$plot.type, colours[,1]),4]))
close.screen(2)

#                                     growth dendro ####
screen(3)

par(mar=c(0,0,0,5), oma=c(0,0,0,0), ps=8, tck=-0.02, mgp=c(3,0.3,0))

growth.dend <- set(growth.dend, "labels_cex", 0.5)

plot(dendrapply(as.dendrogram(growth.dend), 
                function(x) { attr(x, "height") <- log(attr(x, "height")+1); x }), 
     horiz=TRUE, axes=FALSE, yaxs="i", leaflab="none", ylab="", xlab="",
     ylim=c(-1,68.50), xlim=c(2.5,0), lwd=0.5)

dend.pars<-par("usr")
dend.nodes<-as.data.frame(cbind(get_nodes_xy(growth.dend, type = c("rectangle"), center = FALSE,
                                             horiz = TRUE),
                                get_nodes_attr(growth.dend, attribute="label")))

dend.nodes$V1<-as.numeric(as.character(dend.nodes$V1))
dend.nodes<-dend.nodes[!is.na(dend.nodes$V3),]
dend.nodes$group<-growth.split[match(dend.nodes$V3, names(growth.split))]
dend.groups<-do.call("rbind", lapply(split(dend.nodes, f=dend.nodes$group), 
                                     function(x){c(min(x$V1),max(x$V1))}))
group.centers<-tapply(dend.nodes$V1, dend.nodes$group, mean)

growth.group.decode<-cbind(names(group.centers[order(group.centers)]),
                           8:1)
write.csv(growth.group.decode, "./Outputs/growth.groups.decoded.csv")

par(xpd=NA)
rect(xleft=par("usr")[1], xright=relative.axis.point(2.76, "x"),
     ybottom=dend.groups[names(sort(group.centers))[c(2,4,6,8)],1]-0.5, 
     ytop=dend.groups[names(sort(group.centers))[c(2,4,6,8)],2]+0.5,
     col=rgb(0,0,0,0.1), border=NA)

text(x=relative.axis.point(0.08, "x"),
     y=sort(group.centers), labels=8:1, cex=1.2, font=2,
     col=c("grey80", "grey50"))

# SLA
rect(xleft=-0.2, xright=-0.3, ybottom=group.centers-0.5, 
     ytop=group.centers+2.5,
     lty="11")
rect(xleft=-0.2, xright=-0.3, ybottom=group.centers-0.5, 
     ytop=(group.centers-0.5)+growth.SLA, 
     col="black")

# WD
rect(xleft=-0.35, xright=-0.45, ybottom=group.centers-0.5, ytop=group.centers+2.5, lty="11")
rect(xleft=-0.35, xright=-0.45, ybottom=group.centers-0.5, 
     ytop=(group.centers-0.5)+growth.WD, 
     col="black")

# MH
rect(xleft=-0.5, xright=-0.6, ybottom=group.centers-0.5, ytop=group.centers+2.5, lty="11")
rect(xleft=-0.5, xright=-0.6, ybottom=group.centers-0.5, 
     ytop=(group.centers-0.5)+growth.MH, 
     col="black")

draw.ellipse(x=rep(-0.775, sum(names(group.centers) %in% n.groups)), 
             y=group.centers[names(group.centers) %in% n.groups]+1,
             a=0.1, b=1, lwd=1)
text(x=-0.775, y=group.centers[names(group.centers) %in% n.groups]+1, labels="N", adj=0.5, cex=0.8)

segments(x0=-0.15, x1=-0.65, y0=group.centers-0.5, y1=group.centers-0.5)
text(x=-0.25, y=group.centers-1.5, labels="S", adj=0.5, cex=0.8)
text(x=-0.4, y=group.centers-1.5, labels="W", adj=0.5, cex=0.8)
text(x=-0.55, y=group.centers-1.5, labels="H", adj=0.5, cex=0.8)

par(xpd=FALSE)

axis(side=1, mgp=c(3,0,0))
dendro.pars<-par("usr")
mtext(side=1, line=0.75, text="ln(Ward's minimum variance)")
close.screen(3)

#                                     growth points ####
screen(4)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), ps=8, tck=-0.02, mgp=c(3,0.3,0))
plot(x=NULL, y=NULL, xlim=c(-5.6, 1.25), ylim=dend.pars[3:4], 
     yaxs="i", axes=FALSE, xlab="", ylab="")

axis(side=1, at=log(c(0.001,0.01,0.1,1,2,3)), labels=c(0.001,0.01,0.1,1,2,3),
     mgp=c(3,0,0))
axis(side=1, at=log(c(seq(0.001,0.01,0.001), seq(0.01,0.1,0.01), seq(0.1,1,0.1),
                      1,3,1)), labels=NA, tck=-0.01)
mtext(side=1, line=0.75, text="Mean species per plot")

plot.type.pos<-group.centers[growth.preds$group] + 
  c(-2.25, -0.75, 0.75, 2.25)[as.factor(growth.preds$plot.type)]
#plot.type.pos[growth$group==5]<-1.5

segments(x0=growth.preds$estimate + 1.96*growth.preds$se,
         x1=growth.preds$estimate - 1.96*growth.preds$se,
         y0=plot.type.pos, y1=plot.type.pos, lwd=1.5,
         col=rgb(colours[match(growth.preds$plot.type, colours[,1]),2],
                 colours[match(growth.preds$plot.type, colours[,1]),3],
                 colours[match(growth.preds$plot.type, colours[,1]),4]))

points(x=growth.preds$estimate,
       y=plot.type.pos,
       pch=c(25,22,21,24)[as.factor(growth.preds$plot.type)],
       cex=0.8,
       bg=rgb(colours[match(growth.preds$plot.type, colours[,1]),2],
              colours[match(growth.preds$plot.type, colours[,1]),3],
              colours[match(growth.preds$plot.type, colours[,1]),4]))
close.screen(4)

#                                     other screens ####
screen(5)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), ps=8, tck=-0.02, mgp=c(3,0.3,0))
plot.new()
#box()
par(xpd=NA)
text(x=relative.axis.point(0.01, "x"),
     y=relative.axis.point(1.01, "y"), labels="(a) Seed dispersal", font=2, adj=0)
par(xpd=FALSE)
close.screen(5)

screen(6)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), ps=8, tck=-0.02, mgp=c(3,0.3,0))
plot.new()
#box()
par(xpd=NA)
text(x=relative.axis.point(0.01, "x"),
     y=relative.axis.point(1.01, "y"), labels="(b) Growth and structure", 
     font=2, adj=0)
par(xpd=FALSE)
close.screen(6)

screen(7)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), ps=8, tck=-0.02, mgp=c(3,0.3,0))
plot(x=NULL, y=NULL, xlim=c(0,12), ylim=c(0,1), axes=FALSE)

gap<-seq(0,0.15, length=2)[2]

text(x=1.25, y=0.825, labels="Categorial traits", font=2)

text(x=0.5, y=0.8 - rep(gap, 5)*c(1:5),
     labels=c("Wind dispersed", "Elaisome", "Fleshy", "Unassisted",
              "Nitrogen fixer"), adj=0)

PostScriptTrace("./Data/wind.ps")
wind<-readPicture("wind.ps.xml")
grid.picture(wind, 
             x=0.1, 
             y=0.09, 
             width=0.03, height=0.03)

PostScriptTrace("./Data/elaisome.ps")
elaisome<-readPicture("elaisome.ps.xml")
grid.picture(elaisome, 
             x=0.1, 
             y=0.07, 
             width=0.03, height=0.03)

PostScriptTrace("./Data/fleshy.ps")
fleshy<-readPicture("fleshy.ps.xml")
grid.picture(fleshy, 
             x=0.1, 
             y=0.055, 
             width=0.03, height=0.03)

PostScriptTrace("./Data/unassisted.ps")
unassisted<-readPicture("unassisted.ps.xml")
grid.picture(unassisted, 
             x=0.1, 
             y=0.035, 
             width=0.03, height=0.03)

draw.ellipse(x=0.2, 
             y=0.8-gap*5,
             a=0.2, b=0.06, lwd=1)
text(x=0.2, y=0.8-gap*5, labels="N", adj=0.5, cex=0.8)

text(x=5.5, y=0.825, labels="Continuous traits", font=2)

text(x=4.5, y=0.8 - rep(gap, 4)*c(1:4),
     labels=c("Maximum height", "Seed mass", 
              "Specific leaf area", "Wood density"), adj=0)

text(x=4, y=0.8 - rep(gap, 4)*c(1:4),
     labels=c("H", "M", "S", "W"), font=2)

text(x=10, y=0.825, labels="Stand categories", font=2)

text(x=9.5, y=0.8 - rep(gap, 4)*c(1:4),
     labels=c("Remnant", "Planting",
              "Young regrowth", "Old regrowth"), adj=0)

points(x=rep(9, 4), y=0.8 - rep(gap, 4)*c(1:4),
       pch=c(21,22,24,25),
       bg=rgb(colours[1:4,5],
              colours[1:4,6],
              colours[1:4,7]),
       cex=1.2)
close.screen(7)

close.screen(all.screens=TRUE)
dev.off()
