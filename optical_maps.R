library(XML)
library(seqinr)
library(grid)
seq_file = "data/spiro.scf.fasta"
scaffolds <- read.fasta(file = seq_file)
max_chr_size <- length(scaffolds[['chr1']])
# Find the chromosomes non-gap sequnces
in_silico_boxed_regions <- list()
num_chr <- 9
for (num in 1:num_chr) {
  scfid <- paste("chr", num, sep = "")
  index <- which(is.na(match(scaffolds[[scfid]], 'n')))
  intervals <- split(index, cumsum(c(1, diff(index) != 1)))
  starts <- sapply(intervals, min)
  ends <- sapply(intervals, max)
  in_silico_boxed_regions[[num]] <- list(starts = starts, ends = ends)
}

# In silico optical maps
library(purrr)
in_silico_file <- 'data/spiro.scf.fasta.NheI'
in_silico_maps <- as.list(readLines(in_silico_file))
in_silico_maps <- sapply(map(map(in_silico_maps, strsplit, split=' '), unlist), function(x) as.numeric(x))

in_silico_chr_tags <- c()
for (n in 1:num_chr) {
  in_silico_chr_tags <- append(in_silico_chr_tags, c(paste(n, "[NheI] (in silico)")))
}


# Physical optical maps
data <- xmlParse("data/chromosomes.maps")
xml_data <- xmlToList(data)

physical_maps <- list()
j <- 1
physical_maps_order <- c(1, 7, 6, 3, 2, 4, 5, 8, 9)
for (i in 1:length(xml_data)){
  if (names(xml_data[i]) == "RESTRICTION_MAP") {
    info <- xml_data[i]$RESTRICTION_MAP$.attrs
    if (startsWith(info['ID'], "tig")) {
      next
    }
    map <- list()
    print(info['ID'])
    # Get the fragment sizes
    sizes <- sapply(xml_data[i]$RESTRICTION_MAP$FRAGMENTS, `[`, 'S')
    sizes <- sizes[!is.na(sizes)]
    d <- data.frame(fragment_size = sizes, sites = cumsum(sizes))
    map <- append(map, info)
    map <- append(map, d)
    chr_id <- match(j, physical_maps_order)
    physical_maps[[chr_id]] <- map
    j <- j + 1
  }
}

physical_boxed_regions <- list()
for (j in 1:num_chr) {
    physical_boxed_regions[[j]] <- list(starts = 1, ends = max(physical_maps[[j]]$sites))
}


physical_chr_tags <- c()
for (n in 1:num_chr) {
  if (n==2 || n == 3) {
    physical_chr_tags <- append(physical_chr_tags, c(paste(n, "[NheI] (incomplete map, right end low coverage)")))
  } else {
  physical_chr_tags <- append(physical_chr_tags, c(paste(n, "[NheI]")))
  }
}



mapGrob <- function(sites, ...) {
  enzymeGrob <- gList()
  for (i in 1:length(sites)){
    site <- sites[i]
    enzymeGrob[[i]] <- segmentsGrob(site, 0, site, 1, default.units = "native", 
                                    gp = gpar(col="grey40", lwd = 0.3), 
                                    name = paste("enzyme", i, sep="."))
  }
  return(enzymeGrob)
}






pdf(file = "plots/optical_maps.pdf", height = 8, width = 6.75)
grid.newpage()
vp_x = unit(0, "npc") + unit(0.8, "line")
vp_w = unit(1, "npc") - unit(1.4, "line")

in_silico_shift_to_right <- c(-20000, -20000, -5000, 15000, -22000, -25000, -7000, -8000, -5000)

# Draw physical maps first
cell_height <- 1/8*3/5*1/2
text_height <- 1/8*2/5*1/2
section_height <- 1/8 # Draw 9th chromosome together with 8th
for (j in 1:length(physical_maps)) {
  if (j == 9) {
    text_vp_x = vp_x + unit(1, 'mm') + unit(3/5, 'npc')
    text_vp_y = unit(1 - (j - 1 - 1) * section_height, "npc")
    map_vp_x = vp_x + unit(3/5, 'npc')
    map_vp_y = unit(1 - (j - 1 - 1) * section_height - text_height , "npc")
  } else {
    text_vp_x = vp_x + unit(1, 'mm')
    text_vp_y = unit(1 - (j - 1) * section_height, "npc")
    map_vp_x = vp_x 
    map_vp_y = unit(1 - (j - 1) * section_height - text_height , "npc")
  }
  # draw text
  pushViewport(viewport(x = text_vp_x, y = text_vp_y,
                        height=unit(text_height, "npc"), width = vp_w,
                        xscale = c(0, max_chr_size), yscale = c(0, 1),
                        just = c("left", "top"), name = paste("chr", j, sep = ".")))
  grid.text(physical_chr_tags[j], x = 0, y = 0.5, just = c("left", "center"), default.units = "native", gp = gpar(cex=0.75))
  upViewport()
  
  # draw physical map
  pushViewport(viewport(x = map_vp_x, y = map_vp_y , 
                        height=unit(cell_height, "npc"), width = vp_w,
                        xscale = c(0, max_chr_size), yscale = c(0, 1),
                        just = c("left", "top"), name = paste("map", j, sep = ".")))
  # draw backbone
  grid.lines(c(0, max(physical_boxed_regions[[j]]$ends)), 
             c(0.5,0.5), default.units = "native", name = "dna", gp = gpar(lwd = 1.2))

  # draw non-gap boxes
  for (i in 1: length(physical_boxed_regions[[j]]$starts)){
    start <- physical_boxed_regions[[j]]$starts[i]
    end <- physical_boxed_regions[[j]]$ends[i]
    grid.polygon(x = c(start, start, end, end), y = c(0, 1, 1, 0), gp = gpar(fill = "lightgrey", lwd = 0.3, col='grey40'),
                 name = paste("box", j, i, sep="."), default.units= "native")
  }
  # Draw the restriction enzyme highlights
  map <- physical_maps[[j]]
  grid.draw(mapGrob((map$sites)))
  upViewport()
}
  
# Draw in silico chromosomes    
for (j in 1:length(in_silico_maps)) {
  if (j == 9) {
    text_vp_x = vp_x + unit(1, 'mm') + unit(3/5, 'npc')
    text_vp_y = unit(1 - (j -1 - 1) * section_height - text_height - 2 * cell_height, "npc")
    map_vp_x = vp_x + unit(3/5, 'npc')
    map_vp_y = unit(1 - (j - 1 - 1) * section_height - text_height - cell_height, "npc") 
  } else {
    text_vp_x = vp_x + unit(1, 'mm')
    text_vp_y = unit(1 - (j - 1) * section_height - text_height - 2 * cell_height, "npc")
    map_vp_x = vp_x 
    map_vp_y = unit(1 - (j - 1) * section_height - text_height - cell_height, "npc") 
  }
  # draw in silico map
  pushViewport(viewport(x = map_vp_x, y = map_vp_y,
                        height=unit(cell_height, "npc"), width = vp_w,
                        xscale = c(0, max_chr_size), yscale = c(0, 1),
                        just = c("left", "top"), name = paste("map", j, sep = ".")))
  # draw backbone line
  grid.lines(c(in_silico_shift_to_right[j], max(in_silico_boxed_regions[[j]]$ends) + in_silico_shift_to_right[j]),
            c(0.5,0.5), default.units = "native", name = "dna", gp = gpar(lwd = 1.2))
  # draw non-gap boxes
  for (i in 1: length(in_silico_boxed_regions[[j]]$starts)){
    start <- in_silico_boxed_regions[[j]]$starts[i] + in_silico_shift_to_right[j]
    end <- in_silico_boxed_regions[[j]]$ends[i] + in_silico_shift_to_right[j]
    grid.polygon(x = c(start, start, end, end), y = c(0, 1, 1, 0), gp = gpar(fill = "lightgrey", lwd = 0.3, col='grey40'),
                 name = paste("box", j, i, sep="."), default.units= "native")
  }
  
  # Draw the restriction enzyme highlights
  map <- in_silico_maps[[j]]
  grid.draw(mapGrob((map + in_silico_shift_to_right[j])))
  upViewport()
  # draw text
  pushViewport(viewport(x = text_vp_x, y = text_vp_y,
                        height=unit(text_height, "npc"), width = vp_w,
                        xscale = c(0, max_chr_size), yscale = c(0, 1),
                        just = c("left", "top"), name = paste("chr", j, sep = ".")))
  grid.text(in_silico_chr_tags[j], x = 0, y = 0.5, just = c("left", "center"), default.units = "native", gp = gpar(cex=0.75))
  upViewport()
}

# Draw scale
pushViewport(viewport(x = vp_x , y = unit(0.2, "npc") , 
                      height=unit(0.1, "npc"), width = vp_w, 
                      xscale = c(0, max_chr_size), just = c("left", "top"), name = "dna_scale_vp"))

grid.lines(c(max_chr_size-200000, max_chr_size), c(0.6,0.6), default.units = "native", name = "dna_scale", gp = gpar(lwd = 1.2))
grid.text("0.2 Mbp", x = max_chr_size - 100000, y = 0.7, just = "bottom", default.units = "native", gp = gpar(cex=0.75))
upViewport()

pushViewport(viewport(x = unit(0, "npc"), y = unit(1, "npc"), height=unit(1, "npc"),
                      xscale = c(0, max_chr_size), just = c("left", "top"), name = "a"))
upViewport()


dev.off()
