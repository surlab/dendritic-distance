data_path = r"/Users/Gregg/Dropbox (MIT)/2021 Gregg rotation/kyle_data"

collect_images_at_path = (
    r"/Users/Gregg/Documents/MIT/SUR_lab/2021 Gregg rotation/annotated_images"
)
# Images and plots will be placed here to easily perform a basic QC
# and ensure that that distance matricies reflect the desired values


precision = "%1.3f"
# formatting string for numpy savetxt to determin how many decimals to include in the output csvs

conversion_factor = 0.09  # pixels/micron for the image

subfolder_name = "dendritic_distance"
# within each session directory, this subfolder will be created and contain
# the detailed output of the sucessful algorithm (distance matricies, and other measurements)
