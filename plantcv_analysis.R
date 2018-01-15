#################### PlantCV Output Analysis Script ####################################
## source("plantcv_analysis_helper.R")
snapshot.file <- "SnapshotInfo.csv"
barcode.file <- "barcodes.csv"
plantcv.file <- "plantcv_results.csv"
outlier.removal <- FALSE
filter.zeros <- FALSE
lt1a_planting_date = as.POSIXct("2016-05-05")
lt1b_planting_date = as.POSIXct("2016-05-16")
project.code <- "LT1"
source("~/GitHub/plantcvR/plantcv_data_stub.R")

########################################################################################
# Analyze VIS data
########################################################################################
## The main data frames are: traits (main df), traits.avg (aggregate mean df),
## wue.data (water use efficiency), wue.data.agg (water use efficiency aggregate mean), rgr (Relative Growth Rate)
big <- FALSE
first.day <-traits$day[1]
last.day <- traits$day[length(traits$dap)]
first.num.col <- grep("tv_area", names(traits))
last.num.col <- ncol(traits)
######################################## BEGIN: Trait Plotting Function ###########################################

dir.create("pdf_images")
dir.create("images")
dir.create("pdf_images/aggregate")

## Plot trait for each genotype
trait <- "area"
units <- "cm^2"
units.save <- "cm_2"

## pdf(file=file.path("pdf_images", paste("all", trait, "day.pdf", sep = "_"),height=6,width=6, useDingbats=FALSE)
pdf(file=file.path("pdf_images", paste("all", trait, "day", "independent", "scaling.pdf", sep = "_")),height=7,width=7, useDingbats=FALSE)
for(genotype in genotypes) {
  p <- plot.trait(traits, genotype, treatments, "dap", trait, units, limits = c(), breaks = breaks, labels = labels, x.scale = "Days After Planting", print.plot = TRUE)
  ggsave(paste0("images/",paste(genotype, trait, units.save, sep="_"),".png"),  dpi = 300, width = 7, height = 7, units="in")
}
dev.off()

## Facet plot of all genotype traits.
p <- plot.trait(traits, genotypes, treatments, "dap", trait, paste0("(",units,")"), limits = c(), facet = TRUE, x.text.angle = 0, breaks = breaks, labels = labels, x.scale = "Days After Planting", facet.col = length(genotypes) / 10)
if(big) {
  ggsave(paste0("pdf_images/aggregate/",paste("all",trait, units,"facet",sep="_"),".png"), dpi = 200, width = 7*6, height = 7*3, units = "in", limitsize = FALSE)
} else {
  ggsave(paste0("pdf_images/aggregate/",paste("all",trait, units,"facet",sep="_"),".png"), dpi = 200, width = 7, height = 7*2, units = "in")
}

splits <- split(genotypes, ceiling(seq_along(genotypes)/8))
for( index in seq_along(splits) ) {
  p <- plot.trait(traits, splits[[index]], treatments, "dap", trait, paste0("(",units,")"), limits = c(-5, 270), facet = TRUE, x.text.angle = 90, breaks = breaks, labels = labels, x.scale = "Days After Planting", print.plot = FALSE)
  ggsave(paste0("pdf_images/aggregate/",paste("all", trait, "geno", index, units,"facet",sep="_"),".png"), dpi = 300, width = 14, height = 7, units = "in")
}


######################################## END: Trait Plotting Function #############################################

######################################## BEGIN: Heat Map Plotting #################################################

dir.create("pdf_images")
dir.create("images")

#mean.ratio.area.100to75 <- trait.heatmap(traits.avg, "100", "75", "Ratio of 100 to 75", "area", "ratio",
#  save = TRUE, save.name = "mean.ratio.area.100to75")

#mean.ratio.area.100to50 <- trait.heatmap(traits.avg, "100", "50", "Ratio of 100 to 50", "area", "ratio",
#  save = TRUE, save.name = "mean.ratio.area.100to50")

mean.ratio.area.100to30 <- trait.heatmap(traits.avg, "100", "30", "Ratio of 100 to 30", "area", "ratio",
  save = TRUE, save.name = "mean.ratio.area.100to30")

######################################## END: Heat Map Plotting ###################################################

######################################## BEGIN: Water Use Efficiency Plots ########################################

#################### WUE Plots ####################################

## Plot WUE for each genotype
## dir.create("pdf_images")
pdf(file=file.path("pdf_images","all_wue.pdf"),height=7,width=7, useDingbats=FALSE)
## pdf(file=file.path("pdf_images","all_delta_wue.pdf"),height=6,width=6, useDingbats=FALSE)
treatments <- unique(sort(as.character(traits$treatment)))
genotypes <- unique(sort(as.character(traits$genotype)))
for(genotype in genotypes) {
  p <- plot.trait(wue.data, genotype, treatments, "day", "WUE", "(cm^2 / g)", limits = c(), y.scale = "Water Use Efficiency", x.scale = "Days After Planting")
  p <- p + scale_colour_discrete(name = "Treatment", breaks = breaks, labels = labels, limits = levels(as.factor(traits$treatment)))
  print(p)
  ggsave(paste0("images/",paste(genotype,"wue","cm^2_g",sep="_"),".png"),  dpi = 100, width = 7, height = 7, units="in")
}
dev.off()

## Facet plot of the Water Use Efficiency
trait <- "WUE"
units <- "cm^2/g"
print.units <- gsub("/", "_", units)
p <- plot.trait(wue.data, genotypes, treatments, "dap", trait, paste0("(",units,")"), limits = c(), facet = TRUE, x.text.angle = 0, breaks = breaks, labels = labels, x.scale = "Days After Planting", facet.col = length(genotypes) / 10)
if(big) {
  ggsave(paste0("pdf_images/aggregate/",paste("all",trait, print.units,"facet",sep="_"),".png"), dpi = 200, width = 7*6, height = 7*3, units = "in", limitsize = FALSE)
} else {
  ggsave(paste0("pdf_images/aggregate/",paste("all",trait, print.units,"facet",sep="_"),".png"), dpi = 200, width = 7, height = 7*2, units = "in")
}

## Facet plot of all genotype traits.

splits <- split(genotypes, ceiling(seq_along(genotypes)/8))
for( index in seq_along(splits) ) {
  p <- plot.trait(wue.data, splits[[index]], treatments, "dap", "WUE", "(cm ^2 / g)", limits = c(), facet = TRUE, x.text.angle = 90, breaks = breaks, labels = labels, x.scale = "Days After Planting")
  ggsave(paste0("pdf_images/aggregate/",paste("all","WUE", "geno", index, "cm^2_g","facet",sep="_"),".png"), dpi = 300, width = 14, height = 7, units = "in")
}

######################################## BEGIN: Heat Map Plotting #################################################

dir.create("pdf_images")
dir.create("images")

#mean.ratio.WUE.100to75 <- trait.heatmap(wue.data.agg, "100", "75", "Ratio of 100 to 75", "WUE", "ratio",
#  save = TRUE, save.name = "mean.ratio.WUE.100to75", day.col = "day")

#mean.ratio.WUE.100to50 <- trait.heatmap(wue.data.agg, "100", "50", "Ratio of 100 to 50", "WUE", "ratio",
#  save = TRUE, save.name = "mean.ratio.WUE.100to50", day.col = "day")

mean.ratio.WUE.100to30 <- trait.heatmap(wue.data.agg, "100", "30", "Ratio of 100 to 30", "WUE", "ratio",
  save = TRUE, save.name = "mean.ratio.WUE.100to30", day.col = "day")

######################################## END: Heat Map Plotting ###################################################

######################################## BEGIN: Relative Growth Rate ###############################################

rgr.sub <- rgr[!(is.na(rgr$RGR.area)), ]

pdf(file=file.path("pdf_images","all_rgr.pdf"),height=7,width=7, useDingbats=FALSE)
for(genotype in genotypes) {
  plot.trait(rgr.sub, genotype, treatments, "day", "RGR.area", "(cm^2/day)", limits = c(), breaks = breaks, labels = labels, x.scale = "Days After Planting", print.plot = TRUE)
  ggsave(paste0("images/",paste(genotype,"RGR.area","(cm^2_day)",sep="_"),".png"),  dpi = 100, width = 7, height = 7, units="in" )
}
dev.off()

trait <- "RGR.area"
units <- "(cm^2/day)"
print.units <- "cm^2_day"
## Facet plot of all genotype traits.
p <- plot.trait(rgr.sub, genotypes, treatments, "day", trait, paste0("(",units,")"), limits = c(), facet = TRUE, x.text.angle = 0, breaks = breaks, labels = labels, x.scale = "Days After Planting", facet.col = length(genotypes) / 10)
if(big) {
  ggsave(paste0("pdf_images/aggregate/",paste("all",trait, print.units,"facet",sep="_"),".png"), dpi = 200, width = 7*6, height = 7*3, units = "in", limitsize = FALSE)
} else {
  ggsave(paste0("pdf_images/aggregate/",paste("all",trait, print.units,"facet",sep="_"),".png"), dpi = 200, width = 7, height = 7*2, units = "in")
}

splits <- split(genotypes, ceiling(seq_along(genotypes)/8))
for( index in seq_along(splits) ) {
  p <- plot.trait(rgr.sub, splits[[index]], treatments, "day", "RGR.area", "(cm ^2/day)", limits = c(), facet = TRUE, x.text.angle = 90, breaks = breaks, labels = labels, x.scale = "Days After Planting")
  ggsave(paste0("pdf_images/aggregate/",paste("all","RGR.area", "geno", index, "cm^2_day","facet",sep="_"),".png"), dpi = 300, width = 14, height = 7, units = "in")
}

######################################## BEGIN: Analysis Of Data ##################################################

## Get the mean and confidence interval of height per day
height.per.day.100 <- trait.per.day(traits, "height", treatment = c("100"))
height.per.day.30 <- trait.per.day(traits, "height", treatment = c("30"))

## Get the mean and confidence interval of area per day
area.per.day.100 <- trait.per.day(traits, "area", treatment = c("100"))
area.per.day.30 <- trait.per.day(traits, "height", treatment = c("30"))

## Get the trait differences per two days for each treatment specified
#bd21.height.results <- analyze.trait(traits, "height", "Bd21-0", "90", "22.5")
#bd1.height.results <- analyze.trait(traits, "height", "Bd1-1", "90", "22.5")

## Control for multiple testing by controlling the FDR
#bd21.qvalues.height = p.adjust(bd21.height.results$pvalue, method="fdr")
#bd1.qvalues.height = p.adjust(bd1.height.results$pvalue, method="fdr")

######################################## END: Analysis Of Data #####################################################

dir.create("histogram_gif")
xlim <- c(0, 350)
ylim <- c(0.00, 0.07)
## xlim <- c()
## ylim <- c()
for(day in sort(unique(traits.avg$day))) {
  p <- ggplot(traits.avg[traits.avg$day == day, ], aes(x = area, fill = as.factor(treatment))) + geom_density(alpha=0.3)
  p <- p + ggtitle(paste0("Density Curve of Area By Treatment Per Day (",day,")")) +
    scale_x_continuous( name = "Area", limits = xlim) +
    scale_y_continuous(name = "Density", limits = ylim) +
    theme(axis.title.x=element_text(face="bold"), axis.title.y=element_text(face="bold")) +
    scale_fill_discrete(name = "Treatment", breaks = breaks, labels = labels) +
    theme_minimal()
  ggsave(paste0("histogram_gif/",day,"_histogram_area_v_density",".png"), width = 14, height = 7.88, dpi = 300, units = "in", pointsize = 10)
}

## ImageMagick command from VBox
## convert -background white -alpha remove -layers OptimizePlus -delay 50 -quality 100 *.png all_histogram.gif

## More ridiculous (PCA plots)
dir.create("PCA.images")

## Fancier plots with fviz_pca_biplot
for( day in sort(unique(traits.avg$day)) ) {
  traits.avg.day <- traits.avg[traits.avg$day == day, ]
  rownames(traits.avg.day) <- do.call(paste, c(traits.avg.day[c("genotype", "treatment")], sep = "-"))
  shapes.pca <- PCA(traits.avg.day[, numeric.traits.avg ],graph = F)
  p <- fviz_pca_biplot(shapes.pca, label="all", alpha.var = "contrib",
    habillage=traits.avg.day[, "treatment"],addEllipses = T, invisible = "ind") +
    ## habillage= traits.avg.day[, "treatment"],addEllipses = T) +
    xlim(c(-7,11)) +
    ylim(c(-4,4)) +
    ggtitle(paste0("Biplot of variables and individuals on Day ", day)) +
    theme_minimal()
  ggsave(paste0("PCA.images/", paste("genotypes", "pca", "fviz_pca_biplot_no_ind", paste0("day-", day), sep = "_"), ".png"), width = 14, height = 7.88, units = "in", dpi = 300)
}

#### Jberry PCA Plot
sub <- traits
#*************************************************************************************************
# Making corr plot for shapes data continuous time
#*************************************************************************************************
cor.mat <- round(cor(sub[sub$day==first.day, c(first.num.col:last.num.col)]),3)
corrplot(cor.mat, type="lower", order="original", tl.col="black", tl.srt=45)

#*************************************************************************************************
# Perfoming ANOVA to get variance explained
#*************************************************************************************************
## Forcing the ANOVA to be type 3 because defaults are wtf bbq
options(contrasts = c("contr.sum","contr.poly"))
## mod_lm <- lm(data=sub,area~0+(trait+treatment+treatment:trait):dap+dap)
mod_lm <- lm(data=sub,area~0+(genotype+genotype:treatment):day+day)
## drop1(mod_lm)
mod1 <- Anova(mod_lm,type = 2)
afss <- mod1$`Sum Sq`
mod1$`Pr(>F)` <- round(mod1$`Pr(>F)`,4)
print(cbind(mod1,PctExp=afss/sum(afss)*100))

#*************************************************************************************************
# Getting Slopes
#*************************************************************************************************
zero.day <- 12
slope_mod <- glm(data=sub,height~0+treatment:genotype+(treatment:genotype):(as.numeric(day-zero.day)))
## Getting slopes that interact with dap (days after planting)
interacting.dap <- grep("day", names(slope_mod$coefficient))
num.of.treat <- length(unique(sub$treatment))
slopes.of.interest <- slope_mod$coefficients[interacting.dap]
## slopes <- matrix(slopes.of.interest ,ncol= num.of.treat ,byrow = T, rownames.force = T)
slopes <- as.matrix(data.frame(slopes.of.interest), rownames.force = T)
slopes.spl <- matrix(slopes, ncol = num.of.treat, byrow = T)
rownames(slopes.spl) <- rownames(slopes)[seq(from = 1,to = nrow(slopes), by = num.of.treat)]
slopes <- slopes.spl
## rownames(slopes) <- sapply(sapply(names(slope_mod$coefficients[seq(from=91,to=length(slope_mod$coefficients),by=3)]),function(i) strsplit(i,":")[[1]][2]),function(j) strsplit(j,"name")[[1]][2])
rownames(slopes) <- sapply(sapply(rownames(slopes), function(i) strsplit(i,":")[[1]][2]),function(j) strsplit(j,"type")[[1]][2])
colnames(slopes) <- treatments
slopes

#*************************************************************************************************
# Boxplots of treatments per day
#*************************************************************************************************

for( day in sort(unique(traits$day)) ) {
  p <- ggplot(data=traits[traits$day == day, ],aes(genotype,area)) +
    facet_wrap(~ treatment,scales = "free_x",as.table = F) +
    xlab("Plant Age (days)") +
    ylim(c(-5,475))+
    geom_boxplot(aes(color=genotype))+
    scale_color_discrete(breaks = unique(traits$genotype), labels = unique(traits$genotype)) +
    theme_light()+
    theme(strip.background=element_rect(fill="gray50"),
      strip.text.x=element_text(size=14,color="white"),
      strip.text.y=element_text(size=14,color="white"),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) +
    ggtitle(paste0("Boxplot of genotypes within treatments - Day ", day))
  if(big) {
    p <- p + guides(color=FALSE, fill = FALSE)
  }
  ggsave(paste0("images/", paste("genotypes_within_treatments", paste0("day-", day), sep = "_"), ".png"), width = 14, height = 7.88, units = "in", dpi = 300)
}

#*************************************************************************************************
# Testing water effects for each day and plotting
#*************************************************************************************************
traits.list <- split(traits,traits$genotype)
traits.list <- traits.list[sapply(traits.list, function(x) dim(x)[1] > 0)]
## days <- c(first.day:last.day)
days <- sort(unique(traits$day))
treat.test <- treatments
## days <- days[!days %in% c(17)]

## GTT = genotype treatment time
area_gtt_effects <- data.frame(do.call("rbind",
  lapply(traits.list,function(i){
    sapply(days, function(j){
      sub.gtt <- i[i$day==j & i$treatment %in% treat.test,]
      mod_lm <- lm(data=sub.gtt,area~treatment)
      mod_aov <- anova(mod_lm)
      if(is.na(mod_aov$`Pr(>F)`[[1]])) {
        return(FALSE)
      } else if(mod_aov$`Pr(>F)`[1] < 0.05){
        #TRUE
        mod_coef <- as.numeric(c(mod_lm$coefficients[1],mod_lm$coefficients[1]+mod_lm$coefficients[2]))
        if(all(mod_coef==cummax(mod_coef))|all(mod_coef==cummin(mod_coef))){
          mod_aov$`Pr(>F)`[1]
          #TRUE
        }else{
          FALSE
        }
      }else{
        mod_aov$`Pr(>F)`[1]
        #1
        #FALSE
      }
    })
  })
))

colnames(area_gtt_effects) <- as.character(days)
area_gtt_effects$geno <- rownames(area_gtt_effects)
area_gtt_effects_melt <- melt(area_gtt_effects,id="geno")
area_gtt_effects_melt$adj <- p.adjust(area_gtt_effects_melt$value,method="BH")

p <- ggplot(data=area_gtt_effects_melt,aes(x=variable,y=geno))+
  ggtitle("Area by Watering-Effect Temporal Divergence")+
  ylab("")+
  xlab("Plant Age (Days)")+
  scale_fill_gradientn(colors=rev(colorRampPalette(brewer.pal(9,"Greens"))(10)),"value")+
  geom_tile(data=area_gtt_effects_melt,aes(x=variable,y=geno),fill="white",color="darkgray",size=0.5)+
  geom_tile(data=area_gtt_effects_melt[area_gtt_effects_melt$adj < 0.05,],aes(x=variable,y=geno,fill=rescale(adj,to=c(0,1))),color="darkgray",size=0.5) +
  scale_color_manual(values=c("black","black")) +
  theme_light()+
  theme(strip.background=element_rect(fill="gray50"),
    strip.text=element_text(size=14,color="white"))

if(big) {
  ggsave(paste0("pdf_images/aggregate/",paste("Area_by_Watering-Effect_Temporal_Divergence", sep="_"),".png"), dpi = 200, width = 7*5, height = 7*4, units = "in", limitsize = FALSE)
} else {
  ggsave(paste0("pdf_images/aggregate/",paste("Area_by_Watering-Effect_Temporal_Divergence", sep="_"),".png"), dpi = 200, width = 7*2, height = 7, units = "in")
}

## g <- ggplot_build(p)
## vtiles <- sapply(lapply(g$panel$ranges, "[[", "y.major"), length)
## gt <- ggplot_gtable(g)
## panels <- gt$layout$t[grepl("panel", gt$layout$name)]
## gt$heights[panels] <-lapply(vtiles, unit, "null")
## grid.newpage()
## grid.draw(gt)

## Print out data frame csv for scott, merging all of the data frames.
merge.traits.wue <- merge(traits.avg, wue.data.agg, all = TRUE, by = c("genotype", "treatment", "day", "area"))
#final.merge <- merge(merge.traits.wue, rgr, all = TRUE)
final.merge = merge.traits.wue
final.merge.order <- final.merge[with(final.merge, order(group, day)), ]
final.non.num <- names(final.merge.order)[names(final.merge.order) %ni% grep("area|height|WUE|water", names(final.merge.order), value = T)]
final.merge.order <- final.merge.order[, c(final.non.num, setdiff(names(final.merge.order), final.non.num))]

final.merge.name <- paste0("all_traits_mean_", project.code, ".csv")
write.csv(final.merge.order, final.merge.name, row.names = FALSE)

