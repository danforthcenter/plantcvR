# BEGIN: Loading Libraries
# ########################
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(scales)
library(gtools)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(car)
library(reshape2)
library(grid)
library(stringr)
library(gridExtra)
# END: Loading Libraries
# ########################




########################################################################################
# Analyze VIS data
########################################################################################
############################################
# Zoom correction
############################################
zoom.lm = lm(zoom.camera ~ zoom, data=data.frame(zoom=c(1,6000), zoom.camera=c(1,6)))

# Download data for a reference object imaged at different zoom levels
if (!file.exists('zoom_calibration_data.txt')) {
  download.file('http://files.figshare.com/2084101/zoom_calibration_data.txt',
                'zoom_calibration_data.txt')
}

# Read zoom calibrartion data
z.data = read.table(file="zoom_calibration_data.txt", sep="\t", header=TRUE)

# Calculate px per cm
z.data$px_cm = z.data$length_px / z.data$length_cm

# Convert LemnaTec zoom units to camera zoom units
z.data$zoom.camera = predict(object = zoom.lm, newdata=z.data)

# Zoom correction for area
area.coef = coef(nls(log(rel_area) ~ log(a * exp(b * zoom.camera)),
                     z.data, start = c(a = 1, b = 0.01)))
area.coef = data.frame(a=area.coef[1], b=area.coef[2])
area.nls = nls(rel_area ~ a * exp(b * zoom.camera),
               data = z.data, start=c(a=area.coef$a, b=area.coef$b))

# Zoom correction for length
len.poly = lm(px_cm ~ zoom.camera + I(zoom.camera^2),
              data=z.data[z.data$camera == 'VIS SV',])

############################################
# Read data and format for analysis
############################################
# Planting date
lt1a_planting_date = as.POSIXct("2016-05-05")
lt1b_planting_date = as.POSIXct("2016-05-16")

# Read VIS data
vis.data = read.table(file = "plantcv_results.csv", sep = ',', header = TRUE)
# Read metadata
metadata = read.table(file = "lt1.barcodes.csv", sep = ",", header = TRUE)

# Use barcodes to assign genotype and treatment labels
vis.data = merge(vis.data, metadata, by.x = "plantbarcode.x", by.y = "Barcodes")

# Convert timestamp from text to date-time
vis.data$timestamp = ymd_hms(vis.data$timestamp)

# Days after planting
vis.data$dap = NA
vis.data[vis.data$measurementlabel.x == "TM015_F_051616",]$dap = 
  as.numeric(vis.data[vis.data$measurementlabel.x == "TM015_F_051616",]$timestamp - 
               lt1a_planting_date)
vis.data[vis.data$measurementlabel.x == "TM016_F_052716",]$dap = 
  as.numeric(vis.data[vis.data$measurementlabel.x == "TM016_F_052716",]$timestamp - 
               lt1b_planting_date)

# Integer day
vis.data$day = as.integer(vis.data$dap)

# Convert LemnaTec zoom units to camera zoom units
vis.data$zoom = 1
vis.data$zoom.camera = predict(object = zoom.lm, newdata = vis.data)
vis.data$rel_area = predict(object = area.nls, newdata = vis.data)
vis.data$px_cm = predict(object = len.poly, newdata = vis.data)

############################################
# Build traits table
############################################
traits = data.frame(plantbarcode = vis.data$plantbarcode.x, timestamp = vis.data$timestamp,
                    genotype = vis.data$Genotype, dap = vis.data$dap, day = vis.data$day,
                    treatment = vis.data$Treatment)

# Zoom correct TV and SV area
traits$area = (vis.data$sv0_area / vis.data$rel_area) +
  (vis.data$sv90_area / vis.data$rel_area)

# Zoom correct height
traits$height = ((vis.data$sv0_height_above_bound / vis.data$px_cm) +
                   (vis.data$sv90_height_above_bound / vis.data$px_cm)) / 2

traits$group = paste(traits$genotype, traits$treatment, sep = "-")

# Remove day 17 (lots of missing data)
traits = traits[traits$day != 17,]

############################################
# Leaf area analysis
############################################
genotypes = levels(factor(traits$genotype))
for (genotype in genotypes) {
  leaf.plot = ggplot(traits[traits$genotype == genotype,], aes(x = dap, y = area, color = group)) +
    geom_line(aes(group = plantbarcode)) +
    geom_smooth(method = "loess") +
    scale_x_continuous("Days after planting") +
    scale_y_continuous("Projected leaf area (px)") +
    theme_bw()
  ggsave(filename = paste(genotype, ".pdf", sep = ""), plot = leaf.plot)
}

days = as.integer(levels(factor(as.integer(traits$dap))))
#days = days[days<47]
groups = levels(factor(traits$group))

area.means = c()
for (day in days) {
  for (group in groups) {
    area.mean = as.numeric(mean(traits[traits$group == group & 
                                         traits$day == day,]$area))
    if (is.na(area.mean) | is.nan(area.mean) | is.infinite(area.mean)) {
      area.mean = 0
    }
    area.means = c(area.means,area.mean)
  }
}

area.means = rescale(area.means, c(0,1))
area.mat = matrix(area.means, nrow = length(groups), ncol = length(days))
rownames(area.mat) = groups
colnames(area.mat) = days

color.palette = colorRampPalette(c("white", "#006699"), space = "rgb")
pdf(file = "leafarea_heatmap_blue.pdf", width = 8, height = 8, pointsize = 8, useDingbats = FALSE)
heatmap.2(area.mat,
          Rowv = TRUE,
          Colv = FALSE,
          dendrogram = 'row',
          scale = 'none',
          col = color.palette(256),
          trace = 'none',
          symbreaks = FALSE,
          cex.main = 0.75,
          keysize = 1.4)
dev.off()

# Remove misbehaved lines
traits = traits[traits$genotype != 'BAz9504',]
traits = traits[traits$genotype != 'PI330858',]
traits = traits[traits$genotype != 'PI330181',]
traits = traits[traits$genotype != 'PI329440',]
traits = traits[traits$genotype != 'PI329333',]
traits = traits[traits$genotype != 'PI505722',]
traits = traits[traits$genotype != 'PI329338',]
traits$genotype = droplevels(traits$genotype)

genotypes = levels(factor(traits$genotype))
treat.100 = c()
for (day in days) {
  for (genotype in genotypes) {
    area.mean = as.numeric(mean(traits[traits$genotype == genotype & 
                                         traits$day == day & 
                                         traits$treatment == 100,]$area))
    if (is.na(area.mean) | is.nan(area.mean) | is.infinite(area.mean)) {
      area.mean = 0
    }
    treat.100 = c(treat.100, area.mean)
  }
}
treat.100.mat = matrix(treat.100, nrow = length(genotypes), ncol = length(days))
rownames(treat.100.mat) = genotypes
colnames(treat.100.mat) = days

treat.30 = c()
for (day in days) {
  for (genotype in genotypes) {
    area.mean = as.numeric(mean(traits[traits$genotype == genotype & 
                                         traits$day == day & 
                                         traits$treatment == 30,]$area))
    if (is.na(area.mean) | is.nan(area.mean) | is.infinite(area.mean)) {
      area.mean = 0
    }
    treat.30 = c(treat.30, area.mean)
  }
}
treat.30.mat = matrix(treat.30, nrow = length(genotypes), ncol = length(days))
rownames(treat.30.mat) = genotypes
colnames(treat.30.mat) = days

treat.100.mat[treat.100.mat == 0] = 1
treat.30.mat[treat.30.mat == 0] = 1

color.palette = colorRampPalette(c("darkorange", "white", "#006699"), space = "rgb")
pdf(file = "treat100_vs_treat75.pdf", width = 8, height = 8, pointsize = 8, useDingbats = FALSE)
heatmap.2(log(treat.100.mat/treat.30.mat, 2),
          Rowv = TRUE,
          Colv = FALSE,
          dendrogram = 'row',
          scale = "none",
          col = color.palette(256),
          trace = 'none',
          symbreaks = TRUE)
dev.off()
