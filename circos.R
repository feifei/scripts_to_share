library(circlize)

##############
# Reading in WB files
##############
karyotype_file <- "data/karyotype.txt"
gc_file <- "data/gc.txt"
crmp1_file <- "data/crmp1.txt"
crmp2_file <- "data/crmp2.txt"
crp_file <- "data/crp.txt"
rRNA_file <- "data/5srna.txt" #5S purple
h4_file <- "data/h4.txt"
rt_file <- "data/rt.txt"
coding_density_file <- "data/coding_density.txt"
snps_file <- "data/snps.txt"

# uniq_file <- "data/uniq.txt"
link1_file <- "data/blastn.link1.txt"
link2_file <- "data/blastn.link2.txt"


karyotype <- read.table(file = karyotype_file, 
                           colClasses=c("character", "numeric", "numeric", "character", "character"))
gc <- read.table(file = gc_file, 
                    colClasses=c("character", "numeric", "numeric", "numeric"))
crmp1 <- read.table(file = crmp1_file,
                     colClasses=c("character", "numeric", "numeric"))
crmp2 <- read.table(file = crmp2_file,
                      colClasses=c("character", "numeric", "numeric"))
crp <- read.table(file = crp_file,
                      colClasses=c("character", "numeric", "numeric"))
rRNA <- read.table(file = rRNA_file,
                     colClasses=c("character", "numeric", "numeric"))
h4 <- read.table(file = h4_file,
                   colClasses=c("character", "numeric", "numeric"))
rt <- read.table(file = rt_file,
                 colClasses=c("character", "numeric", "numeric"))
# uniq <- read.table(file = uniq_file,
#                       colClasses=c("character", "numeric", "numeric", "character"))
# names(uniq) <- c("chromosome", "start", "end", "col")
# uniq <- uniq[order(uniq$col),] 

coding_density <- read.table(file = coding_density_file,
                                colClasses = c("character", "numeric", "numeric", "numeric"))
snps <- read.table(file = snps_file,
                      colClasses=c("character", "numeric", "numeric", "numeric"))


link1 <- read.table(file = link1_file,
                       colClasses=c("character", "numeric", "numeric"))
link2 <- read.table(file = link2_file,
                       colClasses=c("character", "numeric", "numeric"))


############
# Setting up colors
############
vvlgrey = rgb(240, 240, 240, maxColorValue=255)
vvlgreen = rgb(237, 248, 233, maxColorValue=255)

crmp1_col = rgb(230, 159, 0, maxColorValue=255) # orange
crmp2_col = rgb(175, 141, 195, maxColorValue=255) # purple
crp_col = 

rRNA_col = rgb(86, 180, 233, maxColorValue=255) # sky blue
h4_col = rgb(0, 114, 178, maxColorValue=255) # blue
rt_col = rgb(251, 154, 153, maxColorValue=255) # pink
# rgb(0, 158, 115, maxColorValue=255) # bluish green


link_col = rgb(180, 180, 180, maxColorValue=255) # grey
# link_col = rand_color(nrow(link1), transparency = 0.5)

############
# Setting up some cutoffs and values
############
snps_cutoff <- 10
snps$col <- sapply(snps$V4, function(i) {
  ifelse (i >= snps_cutoff, "red", "dark grey")
})
max_snps <- max(snps$V4)

coding_cutoff <- 0.5
coding_density$col <- sapply(coding_density$V4, function(i) {
  ifelse (i <= coding_cutoff, "green", "dark grey")
})

#############
# circos
#############
pdf(file="plots/circos.pdf", width=6.65, height=6.65, colormodel="cmyk") #, pointsize=9)

circos.clear()
circos.par("start.degree" = 75)
circos.par("gap.degree" = c(2,2,2,2,2,2,2,2,30))
circos.par("track.height" = 0.055)
circos.par("track.margin" = c(0.011, 0.011))
circos.par("points.overflow.warning" = FALSE)
circos.par("cell.padding" = c(0, 0, 0, 0))

# Initialize chr
circos.initializeWithIdeogram(karyotype, plotType = c("ideogram", "labels"), sort.chr = FALSE)
# Highlighting telomeres, since the default is covered by the border
telomeres = karyotype[karyotype$V5=='acen',]
for (row in 1:nrow(telomeres)) {
  chrid <- telomeres[row, "V1"]
  start <- telomeres[row, "V2"]
  end <- telomeres[row, "V3"]

  if(start == 1) {
    circos.segments(1000, 0, 1000, 1, col = "red", sector.index = chrid, lwd = 1)
  } else {
    circos.segments(end-1000, 0, end-1000, 1, col = "red", sector.index = chrid, lwd = 1)
  }
}
text(0, sum(get.cell.meta.data("yplot"))/2, "Chromosomes", cex = 0.85)

# GC
circos.genomicTrack(gc, numeric.column = 4, ylim = c(0,100),
                    panel.fun = function(region, value, ...){
                      circos.genomicLines(region, value, col = "dark grey", lwd = 0.5)
                    }, bg.border = NA, bg.col = vvlgrey)
circos.yaxis(side = "left", at = c(100, 0), labels.cex = 0.5, tick = TRUE, sector.index = "chr1", lwd = 0.8)
text(0, sum(get.cell.meta.data("yplot"))/2, "GC%", cex = 0.85)


# 
# 5S rRNA, crmp1, RT
circos.genomicTrack(rRNA, ylim = c(0,1),
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, ytop = 1, ybottom = 0,
                                         col = rRNA_col, border = rRNA_col)
                    }, bg.border = NA, bg.col = vvlgrey)

rRNA_track_idx <- get.cell.meta.data("track.index")



circos.genomicTrack(rt, ylim = c(0,1),
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, ytop = 1, ybottom = 0,
                                         col = rt_col, border = rt_col,
                                         track.index = rRNA_track_idx)
                    }, bg.border = NA, bg.col = NA,
                    track.index = rRNA_track_idx)

circos.genomicTrack(crmp1, ylim = c(0,1),
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, ytop = 1, ybottom = 0,
                                         col = crmp1_col, border = crmp1_col,
                                         track.index = rRNA_track_idx)
                    }, bg.border = NA, bg.col = NA,
                    track.index = rRNA_track_idx)



text(0, sum(get.cell.meta.data("yplot"))/2, expression("rRNA/RT/CRMP1"), cex = 0.70)

# 
# CRMP2, H4 in one track
circos.genomicTrack(crmp2, ylim = c(0,1),
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, ytop = 1, ybottom = 0,
                                         col = crmp2_col, border = crmp2_col)
                    }, bg.border = NA, bg.col = vvlgrey)

crmp_track_idx <- get.cell.meta.data("track.index")



# circos.genomicTrack(crp, ylim = c(0,1),
#                     panel.fun = function(region, value, ...) {
#                       circos.genomicRect(region, ytop = 1, ybottom = 0,
#                                          col = crp_col, border = crp_col,
#                                          track.index = crmp1_track_idx)
#                     }, bg.border = NA, bg.col = NA,
#                     track.index = crmp1_track_idx)
circos.genomicTrack(h4, ylim = c(0,1),
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, ytop = 1, ybottom = 0,
                                         col = h4_col, border = h4_col,
                                         track.index = crmp_track_idx)
                    }, bg.border = NA, bg.col = NA,
                    track.index = crmp_track_idx)

text(0, sum(get.cell.meta.data("yplot"))/2, expression("CRMP2/H4"), cex = 0.7)


# Coding density
circos.genomicTrack(coding_density, numeric.column = 4, ylim = c(0, 1),
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, col = value$col, cex = 0.2, pch = 19)
                    }, bg.border = NA, bg.col = vvlgreen, track.height = 0.1)

circos.yaxis(side = "left", at = c(0, 1), labels.cex = 0.5, tick = TRUE, sector.index = "chr1", lwd = 0.8)
text(0, sum(get.cell.meta.data("yplot"))/2 + 0.02, "Coding%/", cex = 0.6)
text(0, sum(get.cell.meta.data("yplot"))/2 - 0.03, "5kbp", cex = 0.6)


# SNPs
circos.genomicTrack(snps, numeric.column = 4, ylim = c(0,max(snps$V4)),
                    panel.fun = function(region, value, ...){
                      circos.genomicPoints(region, value, col=value$col, cex = 0.2, pch = 19)
                    }, bg.border = NA, bg.col = vvlgreen, track.height = 0.1)
circos.yaxis(side = "left", at = c(max(snps$V4), 0), labels.cex = 0.5, tick = TRUE, sector.index = "chr1", lwd = 0.8)
text(0, sum(get.cell.meta.data("yplot"))/2 + 0.02, "#SNPs/", cex = 0.6)
text(0, sum(get.cell.meta.data("yplot"))/2 - 0.03, "kbp", cex = 0.6)

# duplication links within genome
circos.genomicLink(link1, link2, col=link_col, border = NA)

circos.clear()

legend("topright", inset = c(0.03, 0.04),
       legend = c("5S rRNA", "RT", "CRMP1", "CRMP2", "Histone H4"),
       col = c(rRNA_col, rt_col, crmp1_col, crmp2_col, h4_col),
       lty = c(1, 1, 1, 1, 1), lwd = 1.5, y.intersp = 0.75,
       cex = 0.8, bty = "n", pt.cex = 0.3, seg.len = 1.5)

# 
# legend("bottomright", inset = c(0.08, 0.05), 
#        legend = c("RT", expression(paste(psi, "RT", sep = "")), 
#                   "18S ", expression(paste(psi, "28S", sep = ""))), 
#        col = c(transpose_col, p_transpose_col, "red", "black"),
#        lty = c(1, 1, 1, 1), 
#        seg.len = c(1.5, 1.5, 1, 0.5),
#        lwd = 1.5, y.intersp = 0.75, 
#        cex = 0.8, bty = "n", pt.cex = 0.3)
# 
legend("bottomright", inset = c(0.06, 0.04),
       legend = c(paste("Coding%<", coding_cutoff), paste("# SNPs>=", snps_cutoff)),
       col = c( "green", "red"),
       lty = c(NA, NA),
       pch = c(19, 19),
       lwd = 1.5, y.intersp = 0.75,
       cex = 0.8, bty = "n", pt.cex = 0.3)

dev.off()