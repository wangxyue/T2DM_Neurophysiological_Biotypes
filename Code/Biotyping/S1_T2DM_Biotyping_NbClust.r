library(NbClust)
library(factoextra)
library(ggplot2)
library(tidyr)
library(dplyr)


INPUT_DATA_PATH <- "path/to/data/T2DM_deviation_z_matrix.csv"
OUTPUT_DIR <- "path/to/results/Cluster_analysis"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  cat("Created output directory:", OUTPUT_DIR, "\n")
}


set.seed(1234)

data <- read.csv(
  INPUT_DATA_PATH,
  header = FALSE,   
  fileEncoding = "UTF-8-BOM"
)
data_clean <- scale(data) 


# ========== Determine Optimal Number of Clusters using NbClust =============================

nb_results <- NbClust(data_clean, 
                      distance = "euclidean", 
                      min.nc = 2, 
                      max.nc = 10, 
                      method = "kmeans", 
                      index = "all")

p <- fviz_nbclust(nb_results) + 
  theme_minimal() + 
  labs(title = "Optimal Number of Clusters (Winner-take-all across 30 indexes)")
print(p)



# ========== Extract the best k value ==========

best_k <- as.numeric(names(which.max(table(nb_results$Best.nc[1, ]))))
cat("The optimal number of clusters (k) determined by majority vote is: ", best_k, "\n")

# Extract indices supporting the best k value
votes <- nb_results$Best.nc[1, ]
indices_supporting_best_k <- names(votes[votes == best_k])
cat(sprintf("Indices supporting k=%d:\n", best_k))
cat(paste(indices_supporting_best_k, collapse = ", "), "\n")



# ========== Final K-means Clustering ==========

final_km <- kmeans(data_clean, centers = best_k, nstart = 10, iter.max = 1000)
final_data <- cbind(biotype = final_km$cluster, data)
biotype_output_path <- file.path(OUTPUT_DIR, "Biotyping_NbClust_Results.csv")
write.csv(final_data, biotype_output_path, row.names = FALSE)
cat("Clustering complete. Biotyping results saved to: ", biotype_output_path, "\n")



# ========== Export Detailed Voting Report for 30 Indices ==========

index_summary <- as.data.frame(t(nb_results$Best.nc))
index_summary$Index_Name <- rownames(index_summary)
index_summary <- index_summary[, c("Index_Name", "Number_clusters", "Value_Index")]
colnames(index_summary) <- c("Index_Name", "Recommended_K", "Index_Value")

index_summary$Recommended_K <- as.character(index_summary$Recommended_K)
index_summary$Recommended_K[index_summary$Recommended_K %in% c("0", "-Inf")] <- "N/A (Graphical)"

print(index_summary)

report_output_path <- file.path(OUTPUT_DIR, "Detailed_Indices_Report.csv")
write.csv(index_summary, report_output_path, row.names = FALSE)

