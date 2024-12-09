## ----setup, include=FALSE, cache=FALSE--------------------
# set global chunk options
library('knitr')
opts_chunk$set(fig.path="wapat-",
               fig.align="center",
               fig.show="hold",
               echo=TRUE,
               results="markup",
               fig.width=10,
               fig.height=10, out.width='\\linewidth',
               out.height='\\linewidth',
               cache=FALSE,
               dev='png',
               concordance=TRUE,
               error=FALSE)
opts_knit$set(aliases = c(h = 'fig.height',
              w = 'fig.width',
              wo='out.width',
              ho='out.height'))
options(replace.assign=TRUE,width=60)
set.seed(9567)


## ----eval=FALSE-------------------------------------------
## ## If devtools is not yet installed, type
## install.packages("devtools")
## 
## ## Install the package badgertub
## devtools::install_github("ClementCalenge/QBloodWB", ref="main")


## ----load-QBloodWB----------------------------------------
library(QBloodWB)


## ----load-dataset-----------------------------------------
data(bloodWB)
str(bloodWB)


## ----extract-data-of-interest-----------------------------
a <- bloodWB[,1:15]
time <- bloodWB$TimeSinceColl


## ----first-pca, fig.width=12, fig.height=6, out.width='\\linewidth', out.height='0.5\\linewidth'----
library(ade4)
pc0 <- dudi.pca(a, scannf=FALSE)
par(mfrow=c(1,2))
barplot(pc0$eig, main="Eigenvalues")
s.label(pc0$li)


## ----second-pca, fig.width=3.93, fig.height=9.84, out.width='0.4\\linewidth', out.height='\\linewidth'----
## Remove WB #19
ab <- a[-19,]

## Carries out the PCA
pc <- dudi.pca(ab, scan=FALSE)

## adds WB 19 as supplementary individual, and insert it in the results:
pc$li <- rbind(pc$li[-19,],suprow(pc, a[19,])$lisup[1,1],
               pc$li[19,,drop=FALSE])

## Load the packages to allow the plot
library(ggplot2)
library(gridExtra)

## The segmented regression (threshold set at 2)
axis1 <- pc$li[,1]
re <- (time>2)*time
mod1 <- lm(axis1~re)

## First plot (screeplot)
df <- data.frame(noval=1:length(pc$eig),Eig=pc$eig)
gra <- ggplot(df)+geom_bar(aes(x=noval,y=Eig), stat="Identity")+
    xlab("Axes of the PCA")+ylab("Eigenvalues")+
    ggtitle("(A)")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    theme(plot.title = element_text(face="bold",hjust=1))
df2 <- pc$li  |> dplyr::mutate(Time=time)
df3 <- pc$co  |> tibble::rownames_to_column()  |>
dplyr::arrange(dplyr::desc(Comp1))
df3$rowname <- factor(df3$rowname, levels=df3$rowname)
grb <- ggplot(df3)+geom_point(aes(x=Comp1, y=rowname))+
    geom_vline(xintercept = 0)+xlab("Correlation with the first component of the PCA")+
    ylab("Biochemical Parameters")+xlim(-1,1)+ggtitle("(B)")+
    theme(plot.title = element_text(face="bold",hjust=1))
grc <- ggplot(df2)+geom_segment(aes(x = 0, y = coefficients(mod1)[1],
                                    xend = 2, yend = coefficients(mod1)[1]),
                                colour="#F8766D", size=1.5)+
    geom_segment(aes(x = 2, y = coefficients(mod1)[1],
                     xend = 6,
                     yend = coefficients(mod1)[1]+coefficients(mod1)[2]*6),
                 colour="#F8766D", linewidth=1.5, alpha=0.8)+
    geom_point(aes(Time, Axis1))+xlab("Time after death")+
    ylab("Score on the first axis of the PCA")+ggtitle("(C)")+
    theme(plot.title = element_text(face="bold",hjust=1))+
    NULL
grid.arrange(gra, grb, grc, ncol=1)


## ---------------------------------------------------------
(su <- round(100*pc$eig[1]/sum(pc$eig)))


## ----rss-different-thresholds, fig.width=5, fig.height=5, out.width='0.5\\linewidth', out.height='0.5\\linewidth'----
ssr <- rssSRCuts(time, axis1, cutlim=c(0,5))
plot(c(0:5), ssr, ty="l", ylim=range(ssr),
     xlab="Threshold",
     ylab="Residual sum of squares")


## ----two-models-segmented-unsegmented---------------------
## Segmented regression
segtime <- (time>2)*time
mod1 <- lm(axis1~segtime)
## Classical regression
mod2 <- lm(axis1~time)

summary(mod1)
summary(mod2)


## ----anova-segmented-unsegmented--------------------------
anova(mod2,mod1)


## ----bootstrap-threshold-choice---------------------------
set.seed(777) ## for reproducibility

bootSR(time, axis1, cutlim=c(1,4), nBoot=1000)

