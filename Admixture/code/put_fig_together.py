### put two png files together by top and bottom, the grids of the two png files should be the same
from PIL import Image
import os
from PIL import Image, ImageDraw, ImageFont

# read the three png files
pie_map_admixture = Image.open(
    "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Admixture/admixtrue_pie_map_plot_crop.png",
    "r",
)

stack_admixture = Image.open(
    "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/"
    "Admixture/stack_bar_plot_admixture_K_4.png",
    "r",
)

PCA = Image.open(
    "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/outlier_analysis/"
    "PCA_RESULTS/pca_region_pc1_pc2_new.png",
    "r",
)

## adding text "A" to the bottom of the pie_map_admixture
# get a font
fnt = ImageFont.truetype("/Library/Fonts/Arial.ttf", 120)
# get a drawing context
d = ImageDraw.Draw(pie_map_admixture)
# draw text, half opacity
d.text((10, 10), "(A)", font=fnt, fill=(0, 0, 0))

## adding text "B" to the bottom of the stack_admixture
# get a drawing context
d = ImageDraw.Draw(stack_admixture)
# draw text, half opacity
d.text((10, 10), "(B)", font=fnt, fill=(0, 0, 0))

## adding text "C" to the bottom of the PCA
# get a drawing context
d = ImageDraw.Draw(PCA)
# draw text, half opacity
d.text((10, 10), "(C)", font=fnt, fill=(0, 0, 0))


# get the size of the pie_map_admixture
width_pie_map_admixture, height_pie_map_admixture = pie_map_admixture.size

# get the size of the stack_admixture
width_stack_admixture, height_stack_admixture = stack_admixture.size

# get the size of the PCA
width_PCA, height_PCA = PCA.size

# create a new image with the size of the pie_map_admixture, stack_admixture and PCA
new_image = Image.new(
    "RGB",
    (width_pie_map_admixture, height_pie_map_admixture + height_stack_admixture + height_PCA),
    (250, 250, 250),
)

# paste the pie_map_admixture on the top of the new image
new_image.paste(pie_map_admixture, (0, 0))

# paste the stack_admixture on the bottom of the new image
new_image.paste(stack_admixture, (0, height_pie_map_admixture))

# paste the PCA on the bottom of the new image
new_image.paste(PCA, (0, height_pie_map_admixture + height_stack_admixture))

# save the new image
new_image.save(
    "/Users/ericliao/Desktop/WNV_project_files/landscape_genetics/Paper_results/Admixture/admixture_pie_map_stack_bar_PCA.png"
)
