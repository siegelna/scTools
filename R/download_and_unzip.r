# Function to download and unzip a remote .gz file
download_and_unzip <- function(download_url, download_dir = "/data/scratch") {
  # Create the target path for the file
  temp_file <- file.path(download_dir, basename(download_url))
  
  # Download the file using wget
  system(paste("wget -O", temp_file, download_url))
  
  # Check if the file was downloaded successfully
  if (!file.exists(temp_file)) {
    stop("Download failed: File does not exist.")
  }
  
  # Unzip the file using gunzip
  system(paste("gunzip", temp_file))
  
  # Return the path to the unzipped file (without the .gz extension)
  unzipped_file <- sub(".gz", "", temp_file)
  return(unzipped_file)
}

## Example usage
# download_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149050/suppl/GSE149050%5FBulk%5FHuman%5FRawCounts.txt.gz"
# unzipped_file <- download_and_unzip(download_url, download_dir = "/data/scratch")

## Print the path to the unzipped file
# print(unzipped_file)
