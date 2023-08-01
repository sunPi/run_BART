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
getFileName <- function(data.path){
  return(tools::file_path_sans_ext(basename(data.path)))
}
checkSanity <- function(x){
  if(shapiro.test(x)$p.value < 0.05){
    sw.res <- "non-normal"
  } else{
    sw.res <- "normal"
  }
  if(ad.test(x)$p.value){
    ad.res <- "non-normal"
  } else{
    ad.res <- "normal"
  } 
  
  nor.test <- data.frame("sw.test"  = rbind("p.value" = signif(shapiro.test(x)$p.value,3),
                                            "normality" = sw.res),
                         "ad.test" = rbind("p.value" = signif(ad.test(x)$p.value,3),
                                           "normality" = ad.res)
  )
  return(nor.test)
}
readPipeData      <- function(data.path) { # This function detects the file format and based on that, it reads in the file
  # Rename this to "readPipeData"  
  # print("-------------------------------")
  # print("Loading the dataset into R...")
  if(!grepl("\\.xlsx$",data.path)){
    pipe.data <- as.data.frame(data.table::fread(data.path))
    attr(pipe.data, "format")      <- ".csv"
    attr(pipe.data, "data.name")   <- basename(data.path)
    return(pipe.data)
  }
  if(grepl("\\.xlsx$", data.path)){
    es <- readxl::excel_sheets(data.path)
    pipe.data <- list()
    for(sheet in es){
      pipe.data <- append(pipe.data, list(as.data.frame(readxl::read_excel(data.path, sheet = sheet))))
    }
    names(pipe.data) <- es
    attr(pipe.data, "format") <- ".xlsx"
    return(pipe.data)
  }
}
buildTrainTestObject <- function(pipe.data, verbose = 1){
  sid   <- pipe.data[[1]]
  label <- pipe.data[[2]]
  input <- as.data.frame(pipe.data[,-c(1:2)])
  
  ind <- createDataPartition(y = 1:length(label), p = 0.2, list = F)
  
  df_test  <- input[ind, ]
  df_train <- input[-ind,]
  
  y_test   <- label[ind]
  y_train  <- label[-ind]
  
  if(verbose){
    paste("Shape of the train data: ")
    paste(dim(df_train))
    paste("Shape of the test data: ")
    paste(dim(df_test))
  }
  
  
  return(list(sid = sid,
              train = list(x = df_train,
                           y = y_train),
              test  = list(x = df_test,
                           y = y_test))
  )
}
writeSummary <- function(model, outfolder = NULL){
  # Define the file path where you want to save the summary
  if(!is.null(outfolder)){
    file_path <- here(outfolder, "bart_summary.txt")
  } else{
    file_path <- here("bart_summary.txt")
  }
  
  # Open a connection to the file
  sink(file_path)
  
  # Print the summary to the file
  summary(model)
  
  # Close the connection to the file
  sink()
}
calculateMetrics <- function(model, test.x, test.y){
  metrics <- list()
  preds <- predict(model, test.x)
  
  if(checkSanity(preds)[[1]][1] < 0.05 || checkSanity(test.y)[[1]][1] < 0.05){
    pred.true.cor <- cor.test(preds, test.y, method = "spearman")
  } else{
    pred.true.cor <- cor.test(preds, test.y)
  }
  metrics$y <- data.frame("predictions"   = preds, 
                          "actual_values" = test.y)
  
  metrics$cor.test <- data.frame("correlation_strength" = pred.true.cor[[4]], 
                                 "p_value"              = signif(pred.true.cor[[3]], 3))
  metrics$r2   <- caret::R2(preds, test.y)
  metrics$RMSE <- caret::RMSE(preds, test.y)
  
  metrics <- lapply(metrics, as.data.frame)
  return(metrics)
}
#---- Package Installation ----
pkgs <- c("metan", "bartMachine", "caret", "here", "readxl", "writexl", "ggplot2", "nortest", "openxlsx", "docopt")
handleRequirements(pkgs)
# s <- colnames(df_train)

# ORIGINAL DATASET with 0.89 correlation
# [1] "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" "YAP1_UP"                                    "YAP1_DN"                                   
# [4] "CORDENONSI_YAP_CONSERVED_SIGNATURE"         "HALLMARK_HYPOXIA"                           "HALLMARK_ANGIOGENESIS"                     
# [7] "HALLMARK_GLYCOLYSIS"                        "HALLMARK_TGF_BETA_SIGNALING"  

#---- Header ----
"Mesothelioma AI Pipeline - Bayesian Additive Regression Trees

Usage: survival_analysis.R [options]

Options:
  -h --help                    Show this screen.
  -d --data_path=<PATH>        Path to the data frame.
  -o --outfolder=<PATH>        Folder where the results are saved.
  -c --do_cv=<BOOLEAN>         If set to T, runs BART with cross-validation on.
  -u --up_seq=<VALUE>          Upper limit for the number of trees to build when, do_cv is 1.
  -n --num_trees=<VALUE>       The number of BART trees to train.
  -v --verbose=<BOOLEAN>       If set to 1 prints verbously.
"-> doc

#---- Arguments ----
arguments <- docopt(doc, quoted_args = TRUE, help = TRUE)
print(arguments)

data.path <- arguments$data_path
outfolder <- arguments$outfolder
do.cv     <- as.logical(arguments$do_cv)
up.seq    <- as.integer(arguments$up_seq)
num.trees <- as.integer(arguments$num_trees)
if(as.integer(arguments$verbose) == 1){
  verbose <- T
} else{
  verbose <- F
}

# do.cv = T
# verbose = T
# up.seq = 75
# num_trees = 50

#---- Read-in the data ----
pipe.data <- readPipeData(data.path)
outfolder <- here(outfolder, "BART_results", getFileName(data.path))
# data.path <- "/home/jan1/Documents/Cancer_Studies_PhD/Deciphering_Oncogenes/Workflow/curated_datasets/medusa150/M150_EMT_CON_DF.CSV"
# pipe.data <- readPipeData(data.path)
# outfolder <- here("/home/jan1/Documents/Cancer_Studies_PhD/Deciphering_Oncogenes/Workflow/BART", "BART_results", getFileName(data.path))
dir.create(outfolder, recursive = T)
# set.seed(42)

#---- Build a BART Training/Testing Object ----
bart.obj <- buildTrainTestObject(pipe.data)

#---- Statistical Exploration ----
# Distribution Analysis #
nor.test <- checkSanity(pipe.data[[2]])

p <- ggplot(pipe.data, aes(x = .data[[colnames(pipe.data)[2]]])) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black") +
  labs(title = "Distribution of the Training Label", x = "Response", y = "Frequency",
       subtitle = paste0("SW test: ", nor.test$sw.test[1], "\n",
                         "AD test: ", nor.test$ad.test[1]))

pdf(here(outfolder, "label_histogram.pdf"))
print(p)
dev.off()

# Correlations Heatmap #
input <- apply(bart.obj$train$x, 2, as.double)
cor.mat <- corr_coef(as.data.frame(input))

pdf(here(outfolder,"heatmap.pdf"))
plot(cor.mat)
dev.off()

#---- BART Main ----
# Set maximum memory available to BART machine
options(java.parameters="-Xmx10000m")

main <- function(outfolder, verbose, do.cv = F, up.seq = 75, num_trees = 50){
  if(do.cv){ # Check the number of optimum trees (instead of using 50 by default)
    bart.model <- bartMachine(bart.obj$train$x, bart.obj$train$y, num_trees = num_trees)
    cv.dir <- here(outfolder, "cv")
    dir.create(cv.dir)
    pdf(here(cv.dir,"rmse_by_num_trees_oos_res.pdf"))
    print(rmse_by_num_trees(bart.model,
                            tree_list=c(seq(25, up.seq, by=5)),
                            num_replicates=5,
                            holdout_pctg = 0.2))
    dev.off()
  } else{
    # run BART main
    bart.model <- bartMachine(bart.obj$train$x, bart.obj$train$y, num_trees = num_trees, serialize = T)
    saveRDS(bart.model, here(outfolder, "bm.RDS"))
    
    if(verbose){
      cat(summary(bart.model)) # Fitting Results
    }
    
    writeSummary(bart.model, outfolder)
    
    # Check for convergence
    pdf(here(outfolder,"convergence.pdf"))
    plot_convergence_diagnostics(bart.model)
    dev.off()
    
    # Check residual properties
    check_bart_error_assumptions(bart.model)
    
    pdf(here(outfolder,"true_vs_predicted.pdf"))
    plot_y_vs_yhat(bart.model, prediction_intervals = TRUE)
    plot_y_vs_yhat(bart.model, Xtest=bart.obj$test$x, ytest=bart.obj$test$y, prediction_intervals = TRUE)
    dev.off()
    
    # Classic Metrics
    metrics <- calculateMetrics(bart.model, bart.obj$test$x, bart.obj$test$y)
    write_xlsx(metrics, here(outfolder, "performance.xlsx"))
    
    if(verbose){
      print(metrics)
    }
    # Create the scatter plot
    sp <- ggplot(metrics$y, aes(x = actual_values, y = predictions)) +
      geom_point(color = "lightcoral", size = 3) +
      geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add the correlation line
      geom_text(aes(x = mean(actual_values), y = mean(predictions), label = paste("Correlation:", round(metrics$cor.test$correlation_strength, 3))),
                color = "black", size = 4, nudge_y = 0.05) +  # Add correlation label
      labs(title = "Predicted vs. Actual Values",
           x = "Actual Values",
           y = "Predicted Values")
    
    pdf(here(outfolder,"true_vs_actual_scatter.pdf"))
    sp
    dev.off()
    
    cat(crayon::green("Done!\n"))
  }
}

# do.cv = T
# verbose = T
# up.seq = 75
# num_trees = 50
main(do.cv = do.cv, outfolder = outfolder, verbose = verbose, up.seq = up.seq)


# library(xgboost)
# library(caret)
# 
# set.seed(123)  # For reproducibility
# # trainIndex <- createDataPartition(data$y, p = 0.7, list = FALSE)
# # trainData <- data[trainIndex, ]
# # testData <- data[-trainIndex, ]
# 
# 
# xgb_model <- xgboost(data = as.matrix(df_train), 
#                      label = y_train, 
#                      objective = "reg:squarederror",
#                      nrounds = 100, 
#                      eta = 0.1,
#                      max_depth = 3)
# 
# predictions <- predict(xgb_model, as.matrix(df_test))
# 
# rmse <- sqrt(mean((y_test - predictions)^2))
# mae <- mean(abs(y_test$y - predictions))
# r_squared <- cor(y_test$y, predictions)^2
# shapiro.test(y_test)
# shapiro.test(predictions)
# 
# cor.test(y_test, predictions)
# 
# print(paste("RMSE:", rmse))
# print(paste("MAE:", mae))
# print(paste("R-squared:", r_squared))

