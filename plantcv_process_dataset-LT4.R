library(lubridate)

experiments = data.frame(name = c("TM021", "TM022"),
                         planting_date = c("2016-11-17", "2016-11-28"),
                         project_code = c("LT4a", "LT4b"))

for (i in 1:nrow(experiments)) {
  # File names
  snapshot.file <- paste(experiments[i,]$name, "SnapshotInfo.csv", sep = ".")
  barcode.file <- paste(experiments[i,]$name, "barcodes.csv", sep = ".")
  plantcv.file <- paste(experiments[i,]$name, "plantcv.results.csv", sep = ".")
  planting_date = as.POSIXct(experiments[i,]$planting_date)
  project.code <- experiments[i,]$project_code
  
  ####################
  # Build traits table
  ####################
  
  # Read VIS data
  snapshot.data <- read.table(file = snapshot.file, sep = ",", header = TRUE)
  barcode.data <- read.table(file = barcode.file, sep = ",", header = TRUE)
  plantcv.data <- read.table(file = plantcv.file, sep = ",", header = TRUE)
  
  ## Fixing plantcv.data's plantbarcode header
  ## All the datasets require a column named genotype, treatment, and plantbarcode
  #names(plantcv.data)[1] <- "plantbarcode"
  names(barcode.data) <- c("plantbarcode", "genotype", "treatment", "replicate")
  names(snapshot.data)[3] <- "plantbarcode"
  
  ## Cleaning snapshot.data.
  ## Testing massive removal.
  
  snapshot.data <- snapshot.data[, names(snapshot.data) %in% 
                                   c("plantbarcode", "timestamp", "weight.before", 
                                     "weight.after", "water.amount", "measurement.label")]
  
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
  zoom.lm <- lm(zoom.camera ~ zoom, data = data.frame(zoom = c(1, 6000), zoom.camera = c(1, 6)))
  
  # Download data for a reference object imaged at different zoom levels
  if (!file.exists('zoom_calibration_data.txt')) {
    download.file('http://files.figshare.com/2084101/zoom_calibration_data.txt',
                  'zoom_calibration_data.txt')
  }
  
  # Read zoom calibrartion data
  z.data <- read.table(file = "zoom_calibration_data.txt", sep = "\t", header = TRUE)
  
  # Calculate px per cm
  z.data$px_cm <- z.data$length_px / z.data$length_cm
  
  # Calculate area for each row
  z.data$area_cm <- ifelse(z.data$reference == z.data$reference[[1]], (13.2*13.2), (13.2*3.7))
  
  # Calculate px**2 per cm**2
  z.data$px2_cm2 <- z.data$area_px / z.data$area_cm
  
  # Convert LemnaTec zoom units to camera zoom units
  z.data$zoom.camera <- predict(object = zoom.lm, newdata = z.data)
  
  # Zoom correction for area
  area.coef <- coef(nls(log(rel_area) ~ log(a * exp(b * zoom.camera)),
                        z.data, start = c(a = 1, b = 0.01)))
  area.coef <- data.frame(a = area.coef[1], b = area.coef[2])
  area.nls <- nls(rel_area ~ a * exp(b * zoom.camera),
                  data = z.data, start = c(a = area.coef$a, b = area.coef$b))
  
  # Zoom correction for length
  len.poly <- lm(px_cm ~ zoom.camera + I(zoom.camera^2),
                 data = z.data[z.data$camera == 'VIS SV',])
  
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
  
  ############################################
  # Build traits table
  ############################################
  traits <- data.frame(plantbarcode = vis.data$plantbarcode, timestamp = vis.data$timestamp,
                       genotype = vis.data$genotype, treatment = vis.data$treatment, replicate = vis.data$replicate,
                       dap = vis.data$dap, day = vis.data$day, imageday = vis.data$imageday)
  
  # Experimental zoom correct TV and SV area
  # traits$tv_area <- vis.data$tv0_area / vis.data$px2_cm2
  traits$sv_area <- (vis.data$sv0_area / vis.data$px2_cm2) +
    (vis.data$sv90_area / vis.data$px2_cm2)
  # traits$area <- traits$tv_area + traits$sv_area
  traits$area <- traits$sv_area
  
  # Adding in hull area
  # traits$tv_hull.area <- vis.data$tv0_hull.area / vis.data$px2_cm2
  traits$sv_hull.area <- (vis.data$sv0_hull.area / vis.data$px2_cm2) +
    (vis.data$sv90_hull.area / vis.data$px2_cm2)
  # traits$hull.area <- traits$tv_hull.area + traits$sv_hull.area
  traits$hull.area <- traits$sv_hull.area
  
  # Zoom correct height
  traits$height <- ((vis.data$sv0_height_above_bound / vis.data$px_cm) +
                      (vis.data$sv90_height_above_bound / vis.data$px_cm)) / 2
  
  traits$group <- paste(traits$genotype,'-',traits$treatment,sep = '')
  
  barcodes <- unique(sort(traits$plantbarcode))
  
  traits <- traits[with(traits, order(day)), ]
  
  refcols <- c("plantbarcode", "timestamp", "genotype", "treatment", "replicate", "group", "dap", "day", "imageday")
  traits <- traits[, c(refcols, setdiff(names(traits), refcols))]
  
  # Write the traits table out
  write.csv(x = traits, file = paste(experiments[i,]$name, "traits.csv", sep = "."), quote = FALSE,
            row.names = FALSE)
  
  ####################
  # Build water table
  ####################
  water.data <- snapshot.gt.data
  
  ## Add a genotype x treatment group column.
  water.data$group <- paste(water.data$genotype, '-', water.data$treatment, sep = '')
  
  ## Add a column for days after planting since the planting date.
  water.data$date <- ymd_hms(water.data$timestamp)
  water.data$dap = NA
  #water.data$date <- as.POSIXct(water.data$timestamp, origin = "1970-01-01")
  water.data$dap = as.numeric(water.data$date - planting_date)
  
  #water.data$dap = as.numeric(water.data$date - planting_date)
  water.data$day <- as.integer(water.data$dap)
  water.data$imageday <- water.data$day
  write.csv(x = water.data, file = paste(experiments[i,]$name, "water.data.csv", sep = "."),
            quote = FALSE, row.names = FALSE)
}


