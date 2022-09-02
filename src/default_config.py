data_path = r"/Users/Gregg/Dropbox (MIT)/2021 Gregg rotation/kyle_data"

collect_summary_at_path = "demo_results"
#Images and plots will be placed here to easily perform a basic QC
#and ensure that that distance matricies reflect the desired values


re_run = False
#This will skip any directories that already contain the dendritic distance subfolder named VVV

subfolder_name = 'demo_dendritic_distance'
#within each session directory, this subfolder will be created and contain
#the detailed output of the sucessful algorithm (distance matricies, and other measurements)


precision = '%1.3f'
#formatting string for numpy savetxt to determin how many decimals to include in the output csvs

conversion_factor = .09 #pixels/micron for the image
