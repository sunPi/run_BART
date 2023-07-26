handleRequirements  <- function(pkgs){ # This function checks package requirements and install them if they are missing
  suppressMessages(if (!require("BiocManager", character.only = TRUE)) { # First check via R BioConductior
    install.packages("BiocManager")
    BiocManager::install()
  } else {
    ipkgs <- sapply(pkgs, function(...) require(..., character.only = TRUE))
    if (any(!ipkgs)) {
      BiocManager::install(pkgs[!ipkgs])
      install.packages(pkgs[!ipkgs])
    } else {
      message("\n\nCool! your machine has everything is needed.\n\n")
    }
  })
  
  print("Loading required packages...")
  library(pacman)
  pacman::p_load(pkgs, install = TRUE, character.only = TRUE) # Check via RCran and other repositories
  return(pacman::p_loaded()) # Return loaded packages
}

#---- Package Installation ----
pkgs <- c("metan", "bartMachine", "caret", "here", "readxl", "ggplot2")
handleRequirements(pkgs)

#---- Header ----
"Mesothelioma AI Pipeline - Bayesian Additive Regression Trees

Usage: survival_analysis.R [options]

Options:
  -h --help                    Show this screen.
  -d --data_path               Path to the data frame.
  -o --outfolder               Folder where the results are saved.
  -v --verbose                 If set to 1 prints verbously.
"-> doc

data.path <- "/home/jan1/Documents/Cancer_Studies_PhD/Deciphering_Oncogenes/Workflow/curated_datasets/M53_EMT_DF.xlsx"
#---- Arguments ----
arguments <- docopt(doc, quoted_args = TRUE, help = TRUE)
print(arguments)

pipe.data <- read_excel(arguments$data.path)
pipe.data <- read_excel(data.path)
outfolder <- here(arguments$outfolder, "BART_results")

outfolder <- here("/home/jan1/Documents/Cancer_Studies_PhD/Deciphering_Oncogenes/Workflow/BART", "BART_results")
dir.create(outfolder, recursive = T)

sid   <- pipe.data[[1]]
label <- pipe.data[[2]]
input <- as.data.frame(pipe.data[,-c(1:2)])

# y <- BostonHousing$medv
# df <- within(BostonHousing, rm(medv))

set.seed(42)

ind <- createDataPartition(y = 1:length(label), p = 0.2, list = F)

df_test  <- input[ind, ]
df_train <- input[-ind,]

y_test   <- label[ind]
y_train  <- label[-ind]

paste("Shape of the train data: ")
paste(dim(df_train))
paste("Shape of the test data: ")
paste(dim(df_test))

# Distribution Analysis #
p <- ggplot(pipe.data, aes(x = Response)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Input Features", x = "Response")

pdf(here(outfolder, "histogram.pdf"))
print(p)
dev.off()

# Correlations Heatmap #
input <- apply(input, 2, as.double)
cor.mat <- corr_coef(as.data.frame(input))

pdf(here(outfolder,"heatmap.pdf"))
plot(cor.mat)
dev.off()

# BART Main #
# Set maximum memory available to BART machine
options(java.parameters="-Xmx10000m")

# Main BART Function
bart_machine <- bartMachine(df_train, y_train)
length(y_train)
nrow(df_train)

summary(bart_machine) # Fitting Results

# Check the number of optimum trees (instead of using 50 by default)
rmse_by_num_trees(bart_machine, 
                  tree_list=c(seq(25, 75, by=5)),
                  num_replicates=3)

# Re-train the maxchine with the optimized number of trees
bart_machine <- bartMachine(df_train, y_train, num_trees=70, seed=42)
summary(bart_machine)
plot_convergence_diagnostics(bart_machine)

# Check residual properties
check_bart_error_assumptions(bart_machine)

#Check performance on train vs test data
plot_y_vs_yhat(bart_machine, prediction_intervals = TRUE)
plot_y_vs_yhat(bart_machine, Xtest=df_test, ytest=y_test, prediction_intervals = TRUE)

    # Other Metrics
metrics <- list()

shapiro.test(y_pred)
shapiro.test(y_test)

y_pred <- predict(bart_machine, df_test)
pred.true.cor <- cor.test(y_pred, y_test)
pred.true.cor <- cor.test(y_pred, y_test, method = "spearman")

r2            <- caret::R2(y_pred, y_test)
RMSE          <- caret::RMSE(y_pred, y_test)



library(xgboost)
library(caret)

set.seed(123)  # For reproducibility
# trainIndex <- createDataPartition(data$y, p = 0.7, list = FALSE)
# trainData <- data[trainIndex, ]
# testData <- data[-trainIndex, ]


xgb_model <- xgboost(data = as.matrix(df_train), 
                     label = y_train, 
                     objective = "reg:squarederror",
                     nrounds = 100, 
                     eta = 0.1,
                     max_depth = 3)

predictions <- predict(xgb_model, as.matrix(df_test))

rmse <- sqrt(mean((y_test - predictions)^2))
mae <- mean(abs(y_test$y - predictions))
r_squared <- cor(y_test$y, predictions)^2
shapiro.test(y_test)
shapiro.test(predictions)

cor.test(y_test, predictions)

print(paste("RMSE:", rmse))
print(paste("MAE:", mae))
print(paste("R-squared:", r_squared))

