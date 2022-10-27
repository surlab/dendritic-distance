import os
import glob

data_path = "/Users/Gregg/Dropbox (MIT)/2021 Gregg rotation/Katya_data/BM019_longitidunal_good_examples/p28/Cell7/Dend2"#"/Users/Gregg/Dropbox (MIT)/2021 Gregg rotation/kyle_data/ASC27/Cell01/ASC27_cell01_a1-5_tseries-018"#r"/Users/Gregg/Dropbox (MIT)/2021 Gregg rotation/kyle_data"

collect_summary_at_path = r"/Users/Gregg/Documents/MIT/SUR_lab/2021 Gregg rotation/annotated_images"
#Images and plots will be placed here to easily perform a basic QC
#and ensure that that distance matricies reflect the desired values


re_run = True
#This will skip any directories that already contain the dendritic distance subfolder named VVV

walk_dirs = True
#This will walk throught the lower level directories, useful to turn off for debugging

subfolder_name = 'dendritic_distance'
#within each session directory, this subfolder will be created and contain
#the detailed output of the sucessful algorithm (distance matricies, and other measurements)

#Not implemented yet - uses 'center'
#neck_source = 'center' #'center' or 'nearest'
#whether to use the center of a spine ROI of the nearst point on the ROI to the dendrite as the source of the spine neck.


precision = '%1.3f'
#formatting string for numpy savetxt to determin how many decimals to include in the output csvs

conversion_factor = .09 #pixels/micron for the image


#Order of the spine rois
#roi_order = ['spine', 'other', 'dendrite'] #for katya
roi_order = ['spine', 'dendrite', 'other'] #for kyle


#file_identifiers
##for kyle
hide = """
def get_paths_from_data_path(data_path):
    spine_roi_path = None
    shaft_roi_path = None

    spine_roi_paths = []
    rois_zip_path = os.path.join(data_path, "*.*")
    rois_zip_path = glob.glob(rois_zip_path)
    for path in rois_zip_path:
        if "dend" in path and ("tif" not in path):
            if not (shaft_roi_path):
                shaft_roi_path = path
        elif (not (spine_roi_path) and ".zip" in path) and ('old' not in path):
            spine_roi_paths.append(path)

    try:
        try_path = os.path.join(data_path, "*stabilized.tif")
        projection_tif_path = glob.glob(try_path)[0]
    except Exception as E:
        try_path = os.path.join(data_path, "*.tif")
        projection_tif_path = glob.glob(try_path)[0]

    return spine_roi_paths, projection_tif_path, shaft_roi_path
"""

# for katya
def get_paths_from_data_path(data_path):
    spine_roi_path = None
    shaft_roi_path = None

    spine_roi_paths = []
    subdir = os.path.join(data_path, 'ROIs')
    rois_zip_path = os.path.join(subdir, "*.*")
    rois_zip_path = glob.glob(rois_zip_path)
    for path in rois_zip_path:
        if "dend" in path and ("tif" not in path):
            if not (shaft_roi_path):
                shaft_roi_path = path
        elif (not (spine_roi_path) and ".zip" in path) and ('old' not in path):
            spine_roi_paths.append(path)

    subdir = os.path.join(data_path, 'REFS')
    try:
        try_path = os.path.join(subdir, "*stabilized.tif")
        projection_tif_path = glob.glob(try_path)[0]
    except Exception as E:
        try_path = os.path.join(subdir, "*.tif")
        projection_tif_path = glob.glob(try_path)[0]

    return spine_roi_paths, projection_tif_path, shaft_roi_path

