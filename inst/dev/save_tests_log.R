# Load the testthat library
library(testthat)

# Open a file connection to write the logs
log_connection <- file("test_log.txt", open = "wt")

# Redirect both standard output and messages to the file
sink(log_connection)
sink(log_connection, type = "message")

# Run the tests for the specified file.
# Make sure your working directory is the root of your package.
tryCatch({
  test_file("tests/testthat/test-mvnormal-cpp.R")
}, finally = {
  # Stop redirecting output and close the file connection
  # This ensures that the connection is closed even if there's an error
  sink(type = "message")
  sink()
  close(log_connection)
})

# The output is now saved in "test_log.txt" in your working directory.
