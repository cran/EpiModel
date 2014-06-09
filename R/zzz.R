.onAttach <- function(...) {
  
  if (!interactive()) return()
  
  vcurr <- utils::packageVersion("EpiModel")
  v <- "1.0"
  
  if (runif(1) > 0.75) {
    msg <- c(
      
      "\n",
      
      paste0("============ Loading EpiModel ", vcurr, " ============"),
      
      "\n"
      
    )
  } else {
    msg <- c(
      
      "\n",
      
      paste0("============ Loading EpiModel ", vcurr, " ============"),
      
      paste0("\nReview recent software updates with:\n", 
             "news(Version == \"", v, "\", package = \"EpiModel\")"),
      
      "\n==============================================",
      
      "\n"
      
    )
  }
  
  packageStartupMessage(msg)
}
