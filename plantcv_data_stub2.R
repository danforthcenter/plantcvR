source("~/GitHub/plantcvR/plantcv_functions.R")

############################################
# Read data and format for analysis
############################################
# Filenames

## snapshot.file <- "SnapshotInfo.csv"
## barcode.file <- "TM015_barcodes.csv"
## plantcv.file <- "TM015_results.csv"

# Read VIS data
snapshot.data <- read.table(file=snapshot.file, sep=",", header=TRUE)
barcode.data <- read.table(file=barcode.file, sep =",", header=TRUE)
plantcv.data <- read.table(file=plantcv.file, sep = ",", header=TRUE)

## Fixing plantcv.data's plantbarcode header
## All the datasets require a column named genotype, treatment, and plantbarcode
#names(plantcv.data)[1] <- "plantbarcode"
names(barcode.data) <- c("plantbarcode", "genotype", "treatment", "replicate")
names(snapshot.data)[3] <- "plantbarcode"

## Cleaning snapshot.data.
## Testing massive removal.

snapshot.data <- snapshot.data[, names(snapshot.data) %in% c("plantbarcode", "timestamp", "weight.before", "weight.after", "water.amount", "measurement.label")]

snapshot.data <- snapshot.data[snapshot.data$weight.before != -1, ]

## Merging snapshot info and genotype treatment data from barcode.data, test to figure out which rows do not merge
snapshot.gt.data <- merge(barcode.data, snapshot.data, all = TRUE)
excluded <- snapshot.gt.data[rowSums(is.na(snapshot.gt.data)) > 0,]

## Fully merging once non-merged rows are inspected
snapshot.gt.data <- merge(barcode.data, snapshot.data, all = FALSE)

## Testing merging with plantcv.data and snapshot.gt.data
plantcv.gt.data <- merge(barcode.data, plantcv.data, all = TRUE, by = "plantbarcode")

excluded <- plantcv.gt.data[rowSums(is.na(plantcv.gt.data)) > 0, ]

## Fully merging once non-merged rows are inspected
plantcv.gt.data <- merge(barcode.data, plantcv.data, all = FALSE)

vis.data <- plantcv.gt.data

## Zoom-calibration
## Zoom correction
############################################
zoom.lm <- lm(zoom.camera ~ zoom, data=data.frame(zoom=c(1,6000), zoom.camera=c(1,6)))

# Download data for a reference object imaged at different zoom levels
if (!file.exists('zoom_calibration_data.txt')) {
  download.file('http://files.figshare.com/2084101/zoom_calibration_data.txt',
                'zoom_calibration_data.txt')
}

# Read zoom calibrartion data
z.data <- read.table(file="zoom_calibration_data.txt", sep="\t", header=TRUE)

# Calculate px per cm
z.data$px_cm <- z.data$length_px / z.data$length_cm

# Calculate area for each row
z.data$area_cm <- ifelse(z.data$reference == z.data$reference[[1]], (13.2*13.2), (13.2*3.7))

# Calculate px**2 per cm**2
z.data$px2_cm2 <- z.data$area_px / z.data$area_cm

# Convert LemnaTec zoom units to camera zoom units
z.data$zoom.camera <- predict(object = zoom.lm, newdata=z.data)

# Zoom correction for area
area.coef <- coef(nls(log(rel_area) ~ log(a * exp(b * zoom.camera)),
                     z.data, start = c(a = 1, b = 0.01)))
area.coef <- data.frame(a=area.coef[1], b=area.coef[2])
area.nls <- nls(rel_area ~ a * exp(b * zoom.camera),
               data = z.data, start=c(a=area.coef$a, b=area.coef$b))

# Zoom correction for length
len.poly <- lm(px_cm ~ zoom.camera + I(zoom.camera^2),
  data=z.data[z.data$camera == 'VIS SV',])

# Experimental zoom correction for area
area.poly <- lm(px2_cm2 ~ zoom.camera + I(zoom.camera^2),
  data = z.data)

# Convert LemnaTec zoom units to camera zoom units
vis.data$zoom = 1
#vis.data$zoom <- vis.data$tv0_zoom
#vis.data$zoom <- as.integer(gsub('z', '', vis.data$zoom))
vis.data$zoom.camera <- predict(object = zoom.lm, newdata = vis.data)
vis.data$rel_area <- predict(object = area.nls, newdata = vis.data)
vis.data$px_cm <- predict(object = len.poly, newdata = vis.data)
vis.data$px2_cm2 <- predict(object = area.poly, newdata = vis.data)

# Planting date
## planting_date <- as.POSIXct("2016-05-05")

# Convert timestamp from text to date-time
vis.data$timestamp = ymd_hms(vis.data$timestamp)

# Days after planting
vis.data$dap = NA
vis.data$dap = as.numeric(vis.data$timestamp - planting_date)


#vis.data$dap <- as.numeric(vis.data$date - planting_date)
vis.data$day <- as.integer(vis.data$dap)

# Adjusted day for odd and even image days
vis.data$imageday <- 0
vis.data$imageday <- vis.data$day

# Hack
vis.data = vis.data[vis.data$day != 17,]

############################################
# Build traits table
############################################
traits <- data.frame(plantbarcode = vis.data$plantbarcode, timestamp = vis.data$timestamp,
  genotype = vis.data$genotype, treatment = vis.data$treatment, replicate = vis.data$replicate,
  dap = vis.data$dap, day = vis.data$day, imageday = vis.data$imageday)

## # Zoom correct TV and SV area
## traits$tv_area <- vis.data$tv0_area / vis.data$rel_area
## traits$sv_area <- (vis.data$sv0_area / vis.data$rel_area) +
##                  (vis.data$sv90_area / vis.data$rel_area)
## traits$area <- traits$tv_area + traits$sv_area

# Experimental zoom correct TV and SV area
traits$tv_area <- vis.data$tv0_area / vis.data$px2_cm2
traits$sv_area <- (vis.data$sv0_area / vis.data$px2_cm2) +
  (vis.data$sv90_area / vis.data$px2_cm2)
traits$area <- traits$tv_area + traits$sv_area

# Adding in hull area
traits$tv_hull.area <- vis.data$tv0_hull.area / vis.data$px2_cm2
traits$sv_hull.area <- (vis.data$sv0_hull.area / vis.data$px2_cm2) +
  (vis.data$sv90_hull.area / vis.data$px2_cm2)
traits$hull.area <- traits$tv_hull.area + traits$sv_hull.area

# Zoom correct height
traits$height <- ((vis.data$sv0_height_above_bound / vis.data$px_cm) +
                    (vis.data$sv90_height_above_bound / vis.data$px_cm)) / 2

traits$group <- paste(traits$genotype,'-',traits$treatment,sep='')

barcodes <- unique(sort(traits$plantbarcode))

traits <- traits[with(traits, order(day)), ]

refcols <- c("plantbarcode", "timestamp", "genotype", "treatment", "replicate", "group", "dap", "day", "imageday")
traits <- traits[, c(refcols, setdiff(names(traits), refcols))]

## Quick and easy removal of duplicated extra day
## traits <- traits[!duplicated(traits[, c("genotype", "treatment", "replicate", "group", "day")]), ]

######################################## BEGIN: Outlier Detection And Removal #####################################

if(outlier.removal) {
  ## Removing outliers.
  ## get the image dataset without outliers
  ## traits <- c("solidity", "sv_area", "tv_area", "extent_x", "extent_y", "height_above_bound", "wue", "volume")
  trait.vec <- c("sv_area", "tv_area", "area", "height")
  traits.out <- traits

  for(trait in trait.vec) {
    traits.out[, trait][which(traits.out[, trait]==0)] <- NA
  }

  system.time(
    dataset <- outlier.detection(traits.out, col.day="imageday", fill=TRUE)
  )

  no.outliers <- sum(traits.out[, trait.vec]!=dataset[, trait.vec], na.rm=TRUE)

  traits <- dataset
}
######################################## END: Outlier Detection And Removal #######################################

if(filter.zeros) {
  ## Removing all traits with all traits <= 0:
  traits <- filter_(traits, paste0(paste(names(traits)[c(9:15)], collapse = " > 0 & "), " > 0"))
}


traits = traits[traits$genotype != 'BAz9504',]
traits = traits[traits$genotype != 'PI330858',]
traits = traits[traits$genotype != 'PI330181',]
traits = traits[traits$genotype != 'PI329440',]
traits = traits[traits$genotype != 'PI329333',]
traits = traits[traits$genotype != 'PI505722',]
traits = traits[traits$genotype != 'PI329338',]
traits$genotype = droplevels(traits$genotype)


#################### Making an average of the replicates ####################
traits.avg <- aggregate(cbind(tv_area, sv_area, area, height, hull.area)~genotype+treatment+day+imageday, traits, mean)
numeric.traits.avg <- grep("area|height", names(traits.avg))
non.num.traits.avg <- seq_along(traits.avg)[seq_along(traits.avg) %ni% numeric.traits.avg]

## Cool statement that can be used as internal function to pass and
## aggregate within function instead of creating many aggregate data frames
## dfr <- aggregate(as.formula(paste(x, " ~ genotype + treatment + day +imageday", sep = "")), dfr, mean)

######################################## BEGIN: Making data frames ################################################

## Variables used throughout the script
treatments <- unique(sort(as.character(traits$treatment)))
breaks <- treatments
labels <- treatments
genotypes <- unique(sort(as.character(traits$genotype)))


######################################## Water Use Efficiency (WUE!) ###############################################
water.data <- snapshot.gt.data

## Add a genotype x treatment group column.
water.data$group <- paste(water.data$genotype,'-',water.data$treatment,sep='')

## Add a column for days after planting since the planting date.
water.data$date <- ymd_hms(water.data$timestamp)
water.data$dap = NA
#water.data$date <- as.POSIXct(water.data$timestamp, origin = "1970-01-01")
water.data$dap = as.numeric(water.data$date - planting_date)

#water.data$dap = as.numeric(water.data$date - planting_date)
water.data$day <- as.integer(water.data$dap)
water.data$imageday <- water.data$day

# Hack
water.data = water.data[water.data$day != 17,]
water.data = water.data[water.data$genotype != 'BAz9504',]
water.data = water.data[water.data$genotype != 'PI330858',]
water.data = water.data[water.data$genotype != 'PI330181',]
water.data = water.data[water.data$genotype != 'PI329440',]
water.data = water.data[water.data$genotype != 'PI329333',]
water.data = water.data[water.data$genotype != 'PI505722',]
water.data = water.data[water.data$genotype != 'PI329338',]
water.data$genotype = droplevels(water.data$genotype)
## water.data[water.data$day %% 2 == 0,]$imageday <- water.data[water.data$day %% 2 == 0,]$day-1

## Fixing the initial watering amount of 0 to be 0.001 or 0.0001
## Not useful right now but maybe later
## water.data[water.data$water.amount < 1 & water.data$dap <= 22, ]$water.amount <- 1

wue.data <- wue.calculate(water.data, traits, genotypes, treatments, "area")
wue.data$group <- paste(wue.data$genotype,'-',wue.data$treatment,sep='')
## wue.data[wue.data$water == 0, ]$water <- wue.data[wue.data$water == 0, ]$water + 1
wue.data$water <- wue.data$water + 1
wue.data$WUE <- wue.data$area / wue.data$water
wue.data$day <- as.integer(wue.data$dap)

wue.data.agg <- aggregate(cbind(area, water, WUE) ~ genotype + treatment + day + group, wue.data, mean)
wue.data.agg$WUE <- wue.data.agg$area / wue.data.agg$water

wue.data <- wue.data[complete.cases(wue.data), ]

wue.delta <- deltaWUE(wue.data, "area")
wue.delta$group <- paste(wue.delta$genotype,'-',wue.delta$treatment,sep='')
wue.delta$delta_water <- wue.delta$delta_water + 1
wue.delta$WUE <- wue.delta$delta_area / wue.delta$delta_water

wue.delta.agg <- aggregate(cbind(delta_area, delta_water, WUE) ~ genotype + treatment + dap + group, wue.delta, mean)

######################################## Relative Growth Rate (RGR) ################################################

numeric.cols <- c(grep("area", names(traits)), grep("height", names(traits)))
non.numeric <- seq_along(traits)[seq_along(traits) %ni% numeric.cols]
names.non.numeric <- names(traits)[non.numeric][names(traits)[non.numeric] %ni% c("dap", "timestamp", "imageday")]
names.numeric.cols <- names(traits)[numeric.cols]
rgr.traits <- ddply(traits, c("day"), function(x) { aggregate(formula(paste0("cbind(",paste(names.numeric.cols, collapse = ", "),") ~ ", paste(names.non.numeric, collapse = " + "))), x, mean) } )
rgr <- ddply(rgr.traits, c("group"), RGR.split)


