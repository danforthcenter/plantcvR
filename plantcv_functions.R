# BEGIN: Function Creation
# #####################################

# BEGIN: Library loading
# #####################################
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
library(lubridate)
# END: Library loading
# #####################################

# Capitalize function
capitalize <- function(x) { paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2))) }

# Negation operator
"%ni%" <- Negate( "%in%" )

# BEGIN: Trait plotting function
# #####################################
plot.trait <- function(dfr, genotype, treatments, x, y, y.units, color = "factor(treatment)",
                       limits = c(),
                       discrete = "Treatment",
                       breaks = c("22.5", "45", "90", "drought", "primer", "recovery"),
                       labels = c("22.5", "45", "Control", "Drought", "Primer", "Recovery"),
                       geno.col = "genotype",
                       treat.col = "treatment",
                       facet = FALSE,
                       facet.col = length(genotypes) / 2,
                       x.text.angle = 0,
                       se = TRUE,
                       y.scale = paste(capitalize(y)),
                       x.scale = paste(capitalize(x)),
                       print.plot = FALSE
) {
  # Plot x vs y
  bplot <- ggplot(dfr[dfr[, geno.col] %in% genotype & (dfr[, treat.col] %in% treatments),], aes_string(x = x, y = y, color = "factor(treatment)"))
  bplot <- bplot + scale_colour_discrete(name = discrete,
                                         breaks = breaks,
                                         labels = labels
  )
  # Plot command
  bplot <- bplot +
    geom_smooth(method = "loess", size = 1, se = se) +
    ggtitle(genotype) +
    scale_x_continuous(name = x.scale) +
    scale_y_continuous(lim = limits, name = paste(y.scale, y.units, sep = " ")) +
    theme_bw() +
    theme(
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    )
  
  if (facet) {
    formula <- as.formula(paste0("~ ", geno.col))
    bplot <- bplot + theme(axis.text.x = element_text(angle = x.text.angle)) + ggtitle("Facet of Genotypes and Treatments")
    bplot <- bplot + facet_wrap(formula, ncol = facet.col)
  }
  
  # Save file
  if ( print.plot ) { print(bplot) }
  return(bplot)
}
# END: Trait plotting function
# #####################################

# BEGIN: Heat Map Function
# #####################################
trait.heatmap <- function(dfr, diff1, diff2, title, measure,
                          type="difference", save = FALSE, save.name = "heatmap",
                          day.col = "imageday",
                          geno.col = "genotype",
                          treat.col = "treatment"
) {
  days = sort(unique(dfr[, day.col]))
  genotypes = levels(factor(dfr[, geno.col]))
  measure.means = c()
  
  for (day in days) {
    for (geno in genotypes) {
      if ( type == "median.diff" | type == "median.ratio") {
        control <- as.numeric(median(dfr[dfr[, geno.col] == geno & dfr[, treat.col] == diff1 &
                                           dfr[, day.col] == day, ][, measure]))
        other <- as.numeric(median(dfr[dfr[, geno.col] == geno & dfr[, treat.col] == diff2 &
                                         dfr[, day.col] == day, ][, measure]))
      } else {
        control = as.numeric(mean(dfr[dfr[, geno.col] == geno & dfr[, treat.col] == diff1 &
                                      dfr[, day.col] == day,][, measure]))
        other = as.numeric(mean(dfr[dfr[, geno.col] == geno & dfr[, treat.col] == diff2 &
                                    dfr[, day.col] == day,][, measure]))
      }
      
      if (is.na(control)) {
        control = 0
      }
      
      if (is.na(other)) {
        other = 0
      }
      
      if (type == "difference") {
        mean_difference <- control - other
        ## mean_difference <- log(abs(mean_difference), 2)
        
        if (is.nan(mean_difference) | is.na(mean_difference) | is.infinite(mean_difference)) {
          mean_difference <-  0
        }
        measure.means <- c(measure.means, mean_difference)
      } else if (type == "ratio") {
        mean_ratio <- as.numeric(control/other)
        
        if (is.infinite(mean_ratio) | is.nan(mean_ratio) | is.na(mean_ratio)) {
          mean_ratio = 0
        }
        
        ## log_ratio=log(abs(mean_ratio),2)
        
        ## if (is.infinite(log_ratio)| is.nan(log_ratio) |is.na(log_ratio)) {
        ##     log_ratio=0
        ## }
        ## measure.means <- c(measure.means, log_ratio)
        measure.means <- c(measure.means, mean_ratio)
        
      } else if (type == "median.diff") {
        median_difference <- control - other
        
        if (is.nan(median_difference) | is.na(median_difference)) {
          median_difference <-  0
        }
        measure.means <- c(measure.means, median_difference)
      } else if (type == "median.ratio") {
        median_ratio <- as.numeric(control/other)
        
        if (is.infinite(median_ratio) | is.nan(median_ratio) | is.na(median_ratio)) {
          median_ratio <-  0
        }
        
        ## log_ratio <- log(abs(median_ratio), 2)
        
        ## if (is.infinite(log_ratio) | is.nan(log_ratio) | is.na(log_ratio)) {
        ##     log_ratio <-  0
        ## }
        ## measure.means <- c(measure.means, log_ratio)
        measure.means <- c(measure.means, median_ratio)
        
      } else if (type == "control") {
        mean.per.day <- as.numeric(control)
        measure.means <- c(measure.means, mean.per.day)
        
      } else if (type == "other") {
        mean.per.day <- as.numeric(other)
        measure.means <- c(measure.means, mean.per.day)
      }
    }
  }
  
  ## measure.means = rescale(measure.means,c(-2,2))
  vis.mat = matrix(measure.means, nrow = length(genotypes), ncol = length(days))
  rownames(vis.mat) = genotypes
  colnames(vis.mat) = days
  vis.mat <- round(vis.mat, 2)
  
  ## print(vis.mat)
  
  ## hc <- hclust(as.dist(1-cor(t(vis.mat))))
  ## plot(hc)
  ## color.palette = colorRampPalette(c("darkorange1","gray","blue", "red", "green", "yellow", "cyan"),space="rgb")
  color.palette <- colorRampPalette(c("darkorange1", "gray", "blue", "green"), space = "rgb")(4 * 100 - 1)
  ## color.palette = colorRampPalette(c("darkorange1","gray","blue"),space="rgb")
  col_breaks <- c(seq(0,3, length = 100), seq(3.1,5, length = 100), seq(5.1, 10, length = 100), seq(10.1,74, length = 100))
  if (save) {
    pdf(file = paste0("pdf_images/",save.name, ".pdf"), width = 8, height = 8, pointsize = 8, useDingbats = FALSE)
  }
  reorg <- heatmap.2(vis.mat,
                     # cellnote = vis.mat,
                     # notecol = "black",
                     main = title,
                     Rowv = TRUE,
                     ## Rowv=as.dendrogram(hc),
                     Colv = FALSE,
                     dendrogram = 'row',
                     scale = 'none',
                     col = color.palette,
                     trace = 'none',
                     symbreaks = FALSE,
                     cex.main = 0.75,
                     ## breaks = col_breaks,
                     keysize = 1.4)
  if (save) {
    dev.off()
    png(file = paste0("pdf_images/aggregate/", save.name, ".png"), width = 1200, height = 900, res = 150, pointsize = 8)
    eval(reorg$call)
    dev.off()
  }
  ## return(list(vis.mat, reorg, hc))
  return(list(vis.mat, reorg))
}
# END: Heat Map Function
# #####################################

# BEGIN: Water Use Efficiency Functions
# #####################################
# Calculate WUE for specified genotypes
wue.calculate <- function(water, vis, genotypes, treatments, trait,
                          geno.col = "genotype",
                          treat.col = "treatment",
                          id.col = "plantbarcode",
                          dap.col = "dap",
                          water.col = "water.amount",
                          group.col = "group",
                          rep.col = "replicate"
) {
  # WUE data vectors
  dap.list <- c()
  plant.list <- c()
  water.list <- c()
  vis.list <- c()
  treatment.list <- c()
  group.list <- c()
  genotype.list <- c()
  replicate.list <- c()
  
  # Get unique barcodes for vis
  barcodes <- unique(vis[(vis[, geno.col] %in% genotypes) & (vis[, treat.col] %in% treatments), id.col])
  for (barcode in barcodes) {
    snapshots <- vis[vis[, id.col] == barcode,]
    snapshots <- snapshots[with(snapshots, order(dap)),]
    for (row in 1:nrow(snapshots)) {
      total.water <- sum(water[water[, dap.col] <= snapshots[row, dap.col] &
                                 water[, id.col] == barcode, water.col])
      dap.list <- c(dap.list, snapshots[row, dap.col])
      plant.list <- c(plant.list, barcode)
      water.list <- c(water.list, total.water)
      vis.list <- c(vis.list, snapshots[row, ][, trait])
      treatment.list <- c(treatment.list, as.character(snapshots[row, treat.col]))
      group.list <- c(group.list,snapshots[row, group.col])
      genotype.list <- c(genotype.list, as.character(snapshots[row, geno.col]))
      replicate.list <- c(replicate.list, snapshots[row, rep.col])
    }
  }
  
  wue.data <- data.frame(plantbarcode = plant.list,
                         dap = dap.list,
                         water = water.list,
                         trait = vis.list,
                         genotype = genotype.list,
                         treatment = treatment.list,
                         replicate = replicate.list
  )
  colnames(wue.data) <- c("plantbarcode", "dap", "water", trait, geno.col, treat.col, rep.col)
  
  return(wue.data)
}

## Change in Water Use Efficiency per day function
deltaWUE <- function(wue.data, trait,
                     id.col = "plantbarcode",
                     geno.col = "genotype",
                     treat.col = "treatment",
                     dap.col = "dap",
                     water.col = "water.amount",
                     group.col = "group",
                     rep.col = "replicate"
) {
  water.list <- c()
  trait.list <- c()
  plant.list <- c()
  dap.list <- c()
  genotype.list <- c()
  treatment.list <- c()
  group.list <- c()
  
  delta_trait <- paste("delta", trait, sep = "_")
  delta_water <- "delta_water"
  barcodes <- unique(sort(wue.data[, id.col]))
  for (barcode in barcodes) {
    gr <- wue.data[wue.data[, id.col] == barcode, ]
    gr[, delta_water] <- c(0, diff(gr$water))
    gr[, delta_trait] <- c(0, diff(gr[, trait]))
    
    dap.list <- c(dap.list, as.numeric(gr$dap))
    plant.list <- c(plant.list, as.character(gr[, id.col]))
    water.list <- c(water.list, gr$delta_water)
    trait.list <- c(trait.list, gr[, delta_trait])
    genotype.list <- c(genotype.list, as.character(gr[, geno.col]))
    treatment.list <- c(treatment.list, as.character(gr[, treat.col]))
  }
  delta.data <- data.frame(plantbarcode = plant.list,
                           dap = dap.list,
                           delta_water = water.list,
                           delta_trait = trait.list,
                           genotype = genotype.list,
                           treatment = treatment.list)
  
  colnames(delta.data) <- c("plantbarcode", "dap", "delta_water", delta_trait, "genotype", "treatment")
  
  delta.data <- delta.data[with(delta.data, order(dap)), ]
  ## delta.data <- delta.data[delta.data$delta_water > 0, ]
  ## delta.data <- delta.data[delta.data[, delta_trait] > 0, ]
  return(delta.data)
}

# WUE plotting function
WUE.plot <- function(dfr, x, y, genotypes, treatments, x.title, y.title,
                     color = "factor(treatment)", limits = c(), color.title = color,
                     main = "All treatments and genotypes",
                     geno.col = "genotype",
                     treat.col = "treatment"
) {
  dfr <- dfr[dfr[, geno.col] %in% genotypes & dfr[, treat.col] %in% treatments, ]
  gplot <- ggplot(dfr, aes_string(x = x, y = y, color = color)) +
    ## geom_point() +
    geom_smooth(method = "loess", se = FALSE) +
    scale_x_continuous(name = x.title) +
    scale_y_continuous(lim = c(limits), name = y.title) +
    theme_bw() +
    theme(
      ## theme(legend.position=c(0.2,0.8),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")) +
    labs(color = color.title)
  gplot <- gplot + ggtitle(main)
  print(gplot)
}

# WUE per genotype function calling the previous one for all treatments
plotWUE <- function(dfr, x, y, genotype, x.title = capitalize(x), y.title = capitalize(y),
                    color = "factor(treatment)", limits = c(), color.title = color) {
  treatments <- c("90", "22.5", "45", "drought", "primer", "recovery")
  bplot <- WUE.plot(dfr, x, y, genotype, treatments, x.title, y.title, color, limits)
  bplot <- bplot + ggtitle(genotype)
  ## # Save file
  ## ## pdf(file=paste(genotype,"_biomass_dap.pdf", sep=''), height=6, width=6, useDingbats=FALSE)
  print(bplot)
  ## dev.off()
}
# END: Water Use Efficiency Functions
# #####################################

# BEGIN: Outlier Testing Function (IAP)
# #####################################
# From IAP
# Function for outlier removal based on Grubbs tests
# fill=TRUE: fill the outlier by mean value
# fill=FALSE: set the outlier by NA
grubbs.test.fill.outliers <- function(x, fill=TRUE) {
  if (is.null(x)) return(x)
  require(outliers)
  if (is.matrix(x)) {
    apply(x, 2, grubbs.test.fill.outliers, fill = fill)
  } else if (is.data.frame(x)) {
    as.data.frame(sapply(x, grubbs.test.fill.outliers, fill = fill))
  } else {
    if (all(is.na(x)) || length(unique(x)) < 2
        || sum(is.na(x)) > length(x)/2
        || sum(!is.na(x)) < 6) return(x)
    while (grubbs.test(x)$p.value < 0.01) {
      if (length(unique(x[!is.na(x)])) < 2) break
      if (fill) { ## filled by mean value
        x[which(x == outlier(x))] <- mean(x[-which(x == outlier(x))], na.rm=TRUE)
      } else {
        ## filled by NA
        x[which(x == outlier(x))] <- NA
      }
    }
    return(x)
  }
}

# check if the column exists
column.check <- function(x, ...) {
  columns <- list(...)
  lapply(columns, function(cn) {
    if (!(cn %in% names(x)))
      stop (paste("Column \"", cn, "\" not found!", sep=""))
  } )
  invisible(TRUE)
}

## Adopt Grubbs?? test to detect outliers based on the assumption of normal distribution
## of phenotypic data points for repeated measures on replicated iplants of a single
## genotype for each trait.
outlier.detection <- function(input, col.day="imageday", col.genotype="genotype",
                              col.treatment="treatment", fill=FALSE) {
  column.check(input, col.day, col.genotype, col.treatment)
  library(tcltk)
  pb <- tkProgressBar("Analysis status", "<Processing>: genotype \"\t\t\" on day \t", 0, 100, 0)
  ind <- 0
  genotypes <- unique(as.character(input[, col.genotype]))
  treatments <- unique(as.character(input[, col.treatment]))
  days <- unique(input[, col.day])
  result <- input
  for (genotype in genotypes) {
    for (day in days) {
      if (length(treatments) > 0) {
        for (treatment in treatments) {
          plant.rows <- which(input[, col.genotype] %in% genotype &
                                input[, col.treatment] %in% treatment &
                                input[, col.day] == day)
          plant.data <- input[plant.rows, trait.vec]
          result[plant.rows, trait.vec] <- grubbs.test.fill.outliers(plant.data, fill=fill)
        }
      } else {
        plant.rows <- which(input[, col.genotype] %in% genotype &
                              input[, col.day] == day)
        plant.data <- input[plant.rows, trait.vec]
        result[plant.rows, trait.vec] <- grubbs.test.fill.outliers(plant.data, fill=fill)
      }
      ## show status in a progress bar
      i <- (ind <- ind + 1)/length(genotypes)/length(days) * 100
      info <- sprintf("%d%% done", round(i))
      msg <- sprintf("<Processing>: genotype \"%s\" on day %s", genotype, day)
      setTkProgressBar(pb, i, info, msg)
    }
  }
  close(pb)
  invisible(result)
}
# END: Outlier Testing Function (IAP)
# #####################################

# BEGIN: Genotype Trait Differences
# #####################################
# For each day, calculate the mean and 95% confidence intervals for each genotype in the control water group
trait.per.day <- function(dfr, trait, treatments = c("90"), day.col = "imageday",
                          geno.col = "genotype", treat.col = "treatment"
) {
  imageday = c()
  genotype = c()
  int.low = c()
  int.up = c()
  est = c()
  
  # Loop through each day and genotype
  for (d in min(as.integer(dfr[, day.col])):max(as.integer(dfr[, day.col]))) {
    if (d %% 2 == 1) {
      for (g in levels(factor(dfr[, geno.col]))) {
        h.data <- dfr[dfr[, geno.col] %in% g & dfr[, treat.col] %in% treatments &
                        as.integer(dfr[, day.col]) == d, trait]
        if (!is.na(sd(h.data))) {
          if (sd(h.data) != 0) {
            imageday <- c(imageday,d)
            genotype <- c(genotype,g)
            test <- t.test(h.data)
            int.low <- c(int.low,test$conf.int[1])
            int.up <- c(int.up,test$conf.int[2])
            est <- c(est,test$estimate)
          }
        }
      }
    }
  }
  
  results <- data.frame(imageday = as.numeric(imageday),
                        genotype = genotype,
                        conf.int.low = int.low,
                        conf.int.up = int.up,
                        mean = est)
  
  return(results)
}
# END: Genotype Trait Differences
# #####################################

# BEGIN: Treatment Trait Differences
# #####################################
# Treatment differences in height
# Statistical analysis for height differences.
analyze.trait <- function(dfr, trait, genotype, treat.1, treat.2, geno.col = "genotype",
                          treat.col = "treatment", day.col = "imageday"
) {
  days = c()
  diff.low = c()
  diff.up = c()
  pvals = c()
  for (day in levels(factor(as.integer(dfr[, day.col])))) {
    day = as.integer(day)
    control = dfr[(as.integer(dfr[, day.col]) == day |
                     as.integer(dfr[, day.col]) == day + 2) &
                    dfr[, geno.col] == genotype &
                    dfr[, treat.col] == treat.1, trait]
    drought = dfr[(as.integer(dfr[, day.col]) == day |
                     as.integer(dfr[, day.col]) == day + 2) &
                    dfr[, geno.col] == genotype &
                    dfr[, treat.col] == treat.2, trait]
    test = t.test(x = control, y = drought)
    days = c(days, day)
    diff.low = c(diff.low, test$conf.int[1])
    diff.up = c(diff.up, test$conf.int[2])
    pvals = c(pvals, test$p.value)
  }
  results = data.frame(imageday = as.numeric(days),
                       conf.int.low = diff.low,
                       conf.int.up = diff.up,
                       pvalue = pvals)
  return(results)
}

# Treatment differences in biomass
# Statistical analysis for trait differences.
drought.response <- function(dfr, genotypes, trait, func,
                             treat1 = "90", treat2 = "22.5", geno.col = "genotype",
                             treat.col = "treatment", day.col = "imageday"
) {
  allowed.funcs <- c("median", "mean")
  if (!func %in% allowed.funcs ) {
    writeLines(paste0("Please choose from: ", paste(allowed.funcs, collapse = ', '), " for the appropriate function.\nThen, rerun the function"))
    return(1)
  }
  # Initialize drought response data frame with days
  drought.resp <- data.frame(day = c(
    seq(min(as.integer(dfr[, day.col])),max(as.integer(dfr[, day.col])),2)))
  
  # Initialize genotypes
  for (g in genotypes) {
    control <- paste(g, 'control', sep = '.')
    drought <- paste(g, 'drought', sep = '.')
    
    # Calculate median trait per treatment per day
    drought.resp[, paste(control)] <- 0
    drought.resp[, paste(drought)] <- 0
    
    ## Estimated Pixel Area Measurements
    for (d in seq(min(as.integer(dfr[, day.col])), max(as.integer(dfr[, day.col])), 2)) {
      ## Using a complete
      control.cmd <- paste0(func, "(dfr[dfr[, geno.col] == g & dfr[, treat.col] == treat1 & as.integer(dfr[, day.col]) == d, trait])")
      drought.cmd <- paste0(func, "(dfr[dfr[, geno.col] == g & dfr[, treat.col] == treat2 & as.integer(dfr[, day.col]) == d, trait])")
      drought.resp[drought.resp$day == d, paste(control)] <- eval(parse(text = control.cmd))
      drought.resp[drought.resp$day == d, paste(drought)] <- eval(parse(text = drought.cmd))
    }
  }
  
  # Calculate trait loss to drought
  days = c()
  response = c()
  genotype = c()
  for (r in 1:nrow(drought.resp)) {
    days = c(days, rep(drought.resp[r, "day"], length(genotypes)))
    for (g in genotypes) {
      control = paste(g, 'control', sep = '.')
      drought = paste(g, 'drought', sep = '.')
      response = c(response, drought.resp[r,paste(drought)] -
                     drought.resp[r,paste(control)])
      genotype = c(genotype, g)
    }
  }
  dr.df = data.frame(day = days, response = response, genotype = as.factor(genotype))
  return(dr.df)
}

# Treatment differences in biomass
treat.diff <- function(dfr, genotypes, trait, func,
                       treat.col = "treatment", day.col = "imageday",
                       treat1 = "90", treat2 = "22.5", geno.col = "genotype"
) {
  allowed.funcs <- c("median", "mean")
  if (!func %in% allowed.funcs ) {
    writeLines(paste0("Please choose from: ", paste(allowed.funcs, collapse = ', '), " for the appropriate function.\nThen, rerun the function"))
    return(1)
  }
  # Initialize drought response data frame with days
  drought.resp <- data.frame(day = c(
    seq(min(as.integer(dfr[, day.col])), max(as.integer(dfr[, day.col])), 2)))
  
  # Initialize genotypes
  for (g in genotypes) {
    control <- paste(g, 'control', sep = '.')
    drought <- paste(g, 'drought', sep = '.')
    
    # Calculate median trait per treatment per day
    drought.resp[, paste(control)] <- 0
    drought.resp[, paste(drought)] <- 0
    
    ## Estimated Pixel Area Measurements
    for (d in seq(min(as.integer(dfr[, day.col])), max(as.integer(dfr[, day.col])), 2)) {
      ## Using a complete
      control.cmd <- paste0(func, "(dfr[dfr[, geno.col] == g & dfr[, treat.col] == treat1 & as.integer(dfr[, day.col]) == d, trait])")
      drought.cmd <- paste0(func, "(dfr[dfr[, geno.col] == g & dfr[, treat.col] == treat2 & as.integer(dfr[, day.col]) == d, trait])")
      drought.resp[drought.resp$day == d, paste(control)] <- eval(parse(text = control.cmd))
      drought.resp[drought.resp$day == d, paste(drought)] <- eval(parse(text = drought.cmd))
    }
  }
  
  # Calculate trait loss to drought
  days = c()
  response = c()
  genotype = c()
  for (r in 1:nrow(drought.resp)) {
    days = c(days, rep(drought.resp[r, "day"], length(genotypes)))
    for (g in genotypes) {
      control = paste(g, 'control', sep = '.')
      drought = paste(g, 'drought', sep = '.')
      response = c(response, drought.resp[r,paste(control)] -
                     drought.resp[r,paste(drought)])
      genotype = c(genotype, g)
    }
  }
  dr.df = data.frame(day = days, response = response, genotype = as.factor(genotype))
  return(dr.df)
}
# END: Treatment Trait Differences
# #####################################

# BEGIN: Other functions
# #####################################
RGR.split <- function(dfr, trait.col = "area", rep.col = "replicate", day.col = "day", geno.col = "genotype", treat.col = "treatment", group.col = "group") {
  dfr.final <- aggregate(as.formula(paste0(trait.col, " ~ ", paste(geno.col, treat.col, day.col, group.col, sep = "+"))), dfr, mean)
  dfr <- dfr[order(dfr[, day.col]), ]
  days <- sort(unique(dfr[, day.col]))
  final <- data.frame(day = sort(unique(dfr[, day.col])))
  final.names <- c("day")
  reps <- sort(unique(dfr[, rep.col]))
  for (rep in reps) {
    df <- dfr[dfr[, rep.col] == rep, ]
    rgr <- log(df[, trait.col])
    if (any(is.infinite(rgr))) {
      next
    }
    df[, "RGR"] <- rgr
    df.final <- merge(df, final, all = TRUE)
    final.names <- c(final.names, paste0("rep.",rep))
    final <- cbind(final, df.final[, "RGR"])
  }
  names(final) <- final.names
  if (ncol(final) == 1) {
    return(NULL)
  } else if ( ncol(final) <= 2 ) {
    final[, "bar.ln.area"] <- final[, -1]
  } else {
    final[, "bar.ln.area"] <- rowMeans(final[, -1], na.rm = TRUE)
  }
  final[, "RGR"] <- c(NA, diff(final[, "bar.ln.area"]))
  dfr.final[, "RGR.area"] <- final[, "RGR"]
  return(dfr.final)
}

tplot.make <- function(plt, filename, ...) {
  
  ## The rest of the passed in arguments go to ggplot
  ## Standard ones are: width = 14, height = 7, units = "in", dpi = 100, pointsize = 10
  ## Making a transparent plot
  tplot <- plt + theme(
    ## panel.border = element_blank(),
    ## legend.key = element_blank(),
    ## axis.ticks = element_blank(),
    ## axis.text.y = element_blank(),
    ## axis.text.x = element_blank(),
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA))
  ggsave(filename, bg = "transparent", ...)
  return(tplot)
}

# Function to create a circle
circle <- function(center = c(0, 0), npoints = 100) {
  r = 1
  tt = seq(0, 2 * pi, length = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

# Function to grab legend from ggplot
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
# END: Other functions
# #####################################

# END: Function Creation
# #####################################