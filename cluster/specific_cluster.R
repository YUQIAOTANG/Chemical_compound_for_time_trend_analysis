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

# transfer to the matrix
data_matrix <- as.matrix(data)

# set target_compound
target_compound <- "6PPD-quinone"  

# calculate the person from target compound to others
cor_values <- apply(data_matrix, 1, function(row) cor(row, data_matrix[target_compound, ], use = "complete.obs"))

# sort the person from up to low
cor_values_sorted <- sort(cor_values, decreasing = TRUE)

# selected the top 10
top_n <- 10  
matching_compounds <- names(cor_values_sorted[1:top_n])

# set the threshhold
cor_threshold <- 0.8  
all_matching_compounds <- names(cor_values_sorted[cor_values_sorted >= cor_threshold])

# 
use_all_matching <- FALSE  
selected_compounds <- if (use_all_matching) all_matching_compounds else matching_compounds

# get the extract results
filtered_data_matrix <- data_matrix[selected_compounds, , drop = FALSE]

# transfer intoExpressionSet to cluster
eset <- new("ExpressionSet", exprs = filtered_data_matrix)

# data clean
eset <- filter.NA(eset)
eset <- fill.NA(eset, mode = "mean")
eset <- filter.std(eset, min.std = 0.25)
eset <- standardise(eset)

# calculate the m
set.seed(123)
m <- mestimate(eset)

# set the cluster number(> 1,< the selected number)
user_k <- min(length(selected_compounds), 2)  

# Mfuzz cluster
cl <- mfuzz(eset, c = user_k, m = m)

# make sure only plot once
while (dev.cur() > 1) dev.off()

# plot
mfuzz.plot(eset, cl = cl, mfrow = c(2, 2), time.labels = colnames(exprs(eset)), min.mem = 0.5)

# output results to csv
output_csv <- file.path(getwd(), "filtered_mfuzz_clusters.csv")
cluster_membership <- data.frame(Name = rownames(filtered_data_matrix), Cluster = cl$cluster)
write.csv(cluster_membership, output_csv, row.names = FALSE)

print(paste("cluster finished, saved as", output_csv))
