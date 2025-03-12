# load R packages
if (!requireNamespace("Mfuzz", quietly = TRUE)) install.packages("Mfuzz")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")

library(Mfuzz)
library(openxlsx)

# read Excel
setwd("/Users/tangyuqiao/Desktop/EAI/downstream/Data Analysis w Tugce/time trend/cluster/")
file_path <- "0yr_dry.xlsx"
data <- read.xlsx(file_path, sheet = 1)

# remove the name column
rownames(data) <- data$Name
data <- data[, -1]

# transfer sheet into ExpressionSet object
data_matrix <- as.matrix(data)
eset <- new("ExpressionSet", exprs = data_matrix)

# filter the data(quality control)
eset <- filter.NA(eset)
eset <- fill.NA(eset, mode = "mean")
eset <- filter.std(eset, min.std = 0.25)

# normalize the data
eset <- standardise(eset)

# select the best number of cluster
set.seed(123)  
m <- mestimate(eset)  
c_range <- 2:10

# calculate the k，and use membership to calculate the density
within_cluster_variation <- sapply(c_range, function(k) {
  cl <- mfuzz(eset, c = k, m = m)
  return(sum(cl$membership^2)) 
})

# select the best k value（max membership²）
best_k <- c_range[which.max(within_cluster_variation)]

#  Mfuzz cluster
cl <- mfuzz(eset, c = best_k, m = m)

# plot
while (dev.cur() > 1) dev.off()


mfuzz.plot(eset, cl = cl, mfrow = c(2, 2), time.labels = colnames(exprs(eset)), min.mem = 0.5)

# output the results in  CSV
output_csv <- file.path(getwd(), "mfuzz_clusters.csv")
cluster_membership <- data.frame(Name = rownames(data_matrix), Cluster = cl$cluster)
write.csv(cluster_membership, output_csv, row.names = FALSE)

print(paste("Mfuzz cluster compelete:", output_csv))
