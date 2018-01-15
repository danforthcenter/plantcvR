library(data.table)
library(reshape)
library(parallel)
library(RSQLite)
############################## Read data from SQLite database ####################
directory <- getwd()
## Name of the database file, assuming it is in the current working directory
db <- "plantcv.sqlite3"
## Database driver
drv <- dbDriver("SQLite")
# Connect to database
conn <- dbConnect(drv, dbname = db)

## Query all metadata and feature data for VIS images
vis.df <- dbGetQuery(conn = conn, 'SELECT * FROM metadata NATURAL JOIN features WHERE imgtype = "VIS"')
## Remove features with no data
vis.df <- vis.df[,colSums(rbind(apply(vis.df[,seq(1, ncol(vis.df))], 2, function(x) length(unique(x))))) != 1]

## Index the data by timestamp
vis.dt <- data.table(vis.df, key = "timestamp")
## Snapshots are defined as each unique timepoint (containing multiple images)
snapshots <- unique(vis.df[, "timestamp"])

## vis.df/vis.dt are long tables where each row is data for a single image
## For downstream analysis, we need to make a wider, shorter table were each row is a snapshot
options(mc.cores = 4)
plantcv.data <- mclapply(snapshots, function(snap) {
  ## sub <- vis.dt[vis.dt$timestamp == snap, ]
  sub <- vis.dt[J(snap)]
  ## sub[sub$frame == "none", ]$frame <- 0
  sub[frame == "none", frame := 0]
  
  rows <- lapply(1:nrow(sub), function(row.idx) {
    row <- sub[row.idx, ]
    regex <- "^(?!(timestamp|treatment|measurementlabel|plantbarcode))"
    names(row) <- gsub(regex, paste0(tolower(row$camera), row$frame, "_", "\\1"), names(row), perl = T)
    return(row)
  })
  merge_recurse(rows, by = "timestamp")
})
plantcv.data <- do.call(rbind, plantcv.data)

## Write the plantcv data to a file
write.csv(plantcv.data, "plantcv_results.csv", row.names = FALSE, quote = FALSE)

## Close database connection
dbDisconnect(conn = conn)
