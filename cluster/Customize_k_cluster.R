# load R packages
if (!requireNamespace("Mfuzz", quietly = TRUE)) install.packages("Mfuzz")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")

library(Mfuzz)
library(openxlsx)

# input data
setwd("/Users/tangyuqiao/Desktop/EAI/downstream/Data Analysis w Tugce/time trend/cluster/")
file_path <- "0yr_dry.xlsx"
data <- read.xlsx(file_path, sheet = 1)

# data clean
rownames(data) <- data$Name
data <- data[, -1]

# transfer ExpressionSet object
data_matrix <- as.matrix(data)
eset <- new("ExpressionSet", exprs = data_matrix)

# data claen
eset <- filter.NA(eset)
eset <- fill.NA(eset, mode = "mean")
eset <- filter.std(eset, min.std = 0.25)

# 数据标准化
eset <- standardise(eset)

# customize the cluster number
user_k <- 4  #number
use_custom_m <- FALSE  
user_m <- 2.0  

# calculate the m if didn't define
set.seed(123)
m <- if (use_custom_m) user_m else mestimate(eset)

#  Mfuzz cluster
cl <- mfuzz(eset, c = user_k, m = m)

# close the plot
while (dev.cur() > 1) dev.off()

# plot
mfuzz.plot(eset, cl = cl, mfrow = c(2, 2), time.labels = colnames(exprs(eset)), min.mem = 0.5)

# output the results to CSV
output_csv <- file.path(getwd(), "mfuzz_clusters_custom.csv")
cluster_membership <- data.frame(Name = rownames(data_matrix), Cluster = cl$cluster)
write.csv(cluster_membership, output_csv, row.names = FALSE)

print(paste("Mfuzz clustering complete:", output_csv))
