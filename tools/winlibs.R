# For details see: https://github.com/rwinlib/gdal2
VERSION <- commandArgs(TRUE)
if(!file.exists(sprintf("../windows/gdal2-%s/include/gdal/gdal.h", VERSION))){
  if(getRversion() < "3.3.0") setInternet2()
  download.file(sprintf("https://github.com/rwinlib/gdal2/archive/v%s.zip", VERSION), "lib.zip", quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}
