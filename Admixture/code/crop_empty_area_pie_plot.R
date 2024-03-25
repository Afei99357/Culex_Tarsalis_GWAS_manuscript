library(magick)

rm(list = ls())

# Read the image
img_path <- "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Admixture/admixtrue_pie_map_plot.png"
img <- image_read(img_path)

# Define the cropping region (adjust the values accordingly)
left <- 0
top <- 1800
width <- image_info(img)$width
height <- image_info(img)$height - 3600

# Crop the image
img_cropped <- image_crop(img, geometry = sprintf("%dx%d+%d+%d", width, height, left, top))

img_cropped

# Save the cropped image
output_path <- "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Admixture/admixtrue_pie_map_plot_crop.png"
image_write(img_cropped, path = output_path)
