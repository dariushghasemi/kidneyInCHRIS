#!/usr/bin/Rscript


# ----------------------------
# take the date and time
format(Sys.time(), "%a %b %e %H:%M:%S %Y")

# date
today.date <- format(Sys.Date(), "%d-%b-%y")
print(paste("Today is: ", today.date))

# ----------------------------
# User-defined variables
BASE.dir  <- "/home/dghasemisemeskandeh/projects/gwas"
OUT.dir   <- paste0(BASE.dir, "/03_gwas_qc/outputs/")
GWAS.file <- paste0(BASE.dir, "/Output/ReformedScheme/eGFRw.log.Res.txt.gz")
OUT.mh    <- paste0(OUT.dir, today.date, "_MH_plot_eGFRw.log.Res_filtered.png")
OUT.qq    <- paste0(OUT.dir, today.date, "_QQ_plot_eGFRw.log.Res_filtered_widen1.png")
separator <- "\t"
max.no.rows <- 2000600 #10000 #

# ----------------------------
# File upload
# ----------------------------
library(data.table)

U <- fread(GWAS.file, header = TRUE, sep = separator, stringsAsFactors = F, nrows=max.no.rows) #
names(U) = c("chr","BEG","END","MARKER_ID","NS", "AC","CALLRATE","GENOCNT","MAF","STAT","PVALUE","BETA","SEBETA","R2") 

# ----------------------------
head(U)
str(U)
# ----------------------------
# Data Preparation
# ----------------------------
#names(U) [names(U)=="#CHROM"] <- "chr1"
names(U) [names(U)=="BEG"] <- "position"
#names(U) [names(U)=="PVALUE"] <- "P.value"

# ----------------------------
# Genomic-Control (GC) Correction
# ----------------------------

# Function to estimate the GC factor lambda on the meta-analysis results
gc_lambda1 <- function(p)
{
   # p = p-value vector
   z <- qnorm(0.5 * p, mean=0.0, sd= 1.0, lower.tail=F)
   median(z * z, na.rm=TRUE) / qchisq(0.5, df=1, lower.tail=F)
}

# Lambda estimation
mylambda <- gc_lambda1(U$PVALUE)
GIF <- print(paste("Lambda =", round(mylambda, 3)))

# print time
format(Sys.time(), "%a %b %e %H:%M:%S %Y")

# Generate GC corrected P-values and standard errors
U$pval_gc <- pchisq(qchisq(U$PVALUE, df=1, lower.tail=F)/mylambda, df=1, lower.tail=F)
#U$se_gc   <- U$StdErr*(mylambda^0.5)


# ----------------------------------------------
# Graphical output
# ----------------------------------------------

#------------ QQ plot PNG -------------

QQplot_png <- function(p, png_file_name = "QQ_plot.png")
{
   obs <- -log(p,10)
   N <- length(obs) ## number of p-values
   
   ## create the null distribution (-log10 of the uniform)
   null <- -log(1:N/N,10)
   MAX <- max(c(obs,null), na.rm=TRUE)
   
   ## the jth order statistic from an Uniform(0,1) sample
   ## has a beta(j,n-j+1) distribution
   ## (Casella & Berger, 2002, 2nd edition, pg 230, Duxbury)
   c95 <- qbeta(0.95,1:N,N-(1:N)+1)
   c05 <- qbeta(0.05,1:N,N-(1:N)+1)
   
   lc95 <- -log(c95,10)
   lc05 <- -log(c05,10)
   
   # Axis labels
   x.lab <- expression(Expected~~-log[10](italic(p-value)), cex.axis=1.5)#
   y.lab <- expression(Observed~~-log[10](italic(p-value)), cex.axis=1.5)#
   
   png(png_file_name, width=1300, height=1300)

   ## plot confidence lines
   plot(null, lc95, ylim=c(0,MAX), xlim=c(0,max(null)), type="l", axes=FALSE, xlab="", ylab="")
   par(new=T, mar=c(7,6,5,1))
   plot(null, lc05, ylim=c(0,MAX), xlim=c(0,max(null)), type="l", axes=FALSE, xlab="", ylab="")
   polygon(c(null, rev(null)), c(lc05, rev(lc95)), col="grey80", border="white")
   lines(c(0,max(null)),c(0, max(null)), lwd=2)

   ## add diagonal
   par(new=T, mar=c(7,6,5,1))

   ## add qqplot
   qqplot(null, 
          obs, 
          ylim=c(0,MAX),
          xlim=c(0,max(null)), 
          main="", 
          pch=16, 
          bg="grey20", 
          xlab="", 
          ylab="", 
          cex=1,
          font=2,
          las=1,           #oriantation of axix lables
          #line=.2,        # Space beween axis title and axis lables
          #cex.main: Size of main title
          #cex.lab=2.5,    #axis labels size (text describing axis)
          #cex.axis=2,    #axis text size (values indicating axis tick labels)
          #mgp = c(4.5, 1, 0)  # Space between axis title and axis lables
          )
  
  title(xlab=x.lab, mgp=c(4.5,1,0), cex.lab=2.5, font=2)
  title(ylab=y.lab, mgp=c(3,1,0),   cex.lab=2.5, font=2)

   legend("topleft", inset=.05, legend = GIF, cex = 2.5)
   #text("topleft", paste("Lambda = "))

   dev.off()
}

#QQplot_png(U$PVALUE, OUT.qq) # original

QQplot_png_org <- function(p, png_file_name = "QQ_plot.png")
{
   obs <- -log(p,10)
   N <- length(obs) ## number of p-values
   
   ## create the null distribution (-log10 of the uniform)
   null <- -log(1:N/N,10)
   MAX <- max(c(obs,null), na.rm=TRUE)
   
   ## the jth order statistic from an Uniform(0,1) sample
   ## has a beta(j,n-j+1) distribution
   ## (Casella & Berger, 2002, 2nd edition, pg 230, Duxbury)
   c95 <- qbeta(0.95,1:N,N-(1:N)+1)
   c05 <- qbeta(0.05,1:N,N-(1:N)+1)
   
   lc95 <- -log(c95,10)
   lc05 <- -log(c05,10)
   
   # Axis labels
   x.lab <- expression(Expected~~-log[10](italic(p-value)),cex.axis=1.5)
   y.lab <- expression(Observed~~-log[10](italic(p-value)),cex.axis=1.5)
   
   png(png_file_name, width=1500, height=1500)
   
   ## plot confidence lines
   plot(null, lc95, ylim=c(0,MAX), xlim=c(0,max(null)), type="l", axes=FALSE, xlab="", ylab="")
   par(new=T)
   plot(null, lc05, ylim=c(0,MAX), xlim=c(0,max(null)), type="l", axes=FALSE, xlab="", ylab="")
   polygon(c(null,rev(null)), c(lc05, rev(lc95)), col="grey80", border="white")
   lines(c(0,max(null)),c(0, max(null)), lwd=2)
   
   ## add diagonal
   par(new=T)#, mar=c(5,6,2,2)+2
   
   ## add qqplot
   qqplot(null,obs, ylim=c(0,MAX),xlim=c(0,max(null)), main="", pch=16, bg="grey20", xlab="", ylab="", cex=1, las=1, cex.axis=1.5, cex.lab=4.5)
   legend("topleft", inset=.05, legend = GIF, cex = 2.5)
   title(xlab=x.lab, mgp=c(4.5,1,0), line=4, cex.lab=2.5, font=2)
   title(ylab=y.lab, mgp=c(4.5,5,0), line=6, cex.lab=2.5, font=2)

   dev.off()
}

QQplot_png_org(U$PVALUE, OUT.qq) # original

# print time
format(Sys.time(), "%a %b %e %H:%M:%S %Y")

quit()

#------------- Manhattan plot PNG -------------

MH_plot_PNG <- function(A, MH_file_name = "MH_plot.png", mytitle = "", color1 = "deeppink1", color2 = "deepskyblue1", pch.type = 16)
  {
    A <- A[order(A$chr, A$position),]
    
    A$log10p <- -log10(A$pval_gc)
    n.snps <- nrow(A)
    ymax <- max(A$log10p, na.rm=TRUE)+1
    A$abs_position <- 1:n.snps
    
    png(MH_file_name, width=2400, height=1400)
    
    par(mar=c(3,5,1,1))
    plot(A$abs_position, A$log10p, type="n", yaxt="n", xaxt="n", xlab="", ylab="", main=mytitle, xlim=c(0,n.snps), ylim=c(0,ymax))
    chr <-  1:22
    middle.points <- numeric(22)
    for (i in 1:22)
    {
      idx <- A$chr==i
      points(A$abs_position[idx], A$log10p[idx], col=ifelse(i%%2==0, color1, color2), pch=pch.type, cex=0.8)   
      middle.points[i] <- median(A$abs_position[idx], na.rm=TRUE)
    }
    axis(side=2, at=seq(from=0, to=ymax, by=1), labels=seq(from=0, to=ymax, by=1), tick=T, cex.axis=1, las=1)
    axis(side=1, at=middle.points, labels=as.character(c(1:22)), tick=F, cex.axis=1.5, las=1, pos=0)
    mtext(text="Chromosomes", side=1, cex=1.5, line=1.3)
    mtext(text=expression(-log[10]* (p-value)), side=2, line=3, cex=1.5)
    abline(h=round(-log10(.05/343),0), lty=2, col="black")
    #abline(h=-log10(5*10^-8),lty=2,col="black")
    
    dev.off()
  }
  
  
#MH_plot_PNG(U, OUT.mh)

### sbatch --wrap 'Rscript 03-2_gwas_qc.R'  -c 4 --mem-per-cpu=8GB -J "03-2.R"
