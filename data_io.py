import glob
from config import subfolder_name, collect_images_at_path, precision
from read_roi import read_roi_file, read_roi_zip
import os
import PIL
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd


################
#data in
#################

def get_paths_from_data_path(data_path):

    rois_zip_path = os.path.join(data_path, '*.zip')
    rois_zip_path = glob.glob(rois_zip_path)[0]
    try:
        try_path = os.path.join(data_path, '*stabilized.tif')
        projection_tif_path = glob.glob(try_path)[0]
    except Exception as E:
        try_path = os.path.join(data_path, '*.tif')
        projection_tif_path = glob.glob(try_path)[0]
    shaft_roi_path = None
    try:
        try_path = os.path.join(data_path, '*shaft.roi')
        shaft_roi_path = glob.glob(try_path)[0]
    except Exception as E:
        pass

    return rois_zip_path, projection_tif_path, shaft_roi_path


def load_all(data_path):

    rois_zip_path, projection_tif_path, shaft_roi_path = get_paths_from_data_path(data_path)

    kyle_rois = read_roi_zip(rois_zip_path)

    projection_tif = PIL.Image.open(projection_tif_path)
    shaft_roi = None
    if shaft_roi_path:
        shaft_roi = read_roi_file(shaft_roi_path)

    return kyle_rois, projection_tif, shaft_roi


def seperate_kyle_rois(kyle_rois):
  spine_rois = []
  dend_rois = []
  counter=0
  for name, roi in kyle_rois.items():
    #plt.plot(roi['x'], roi['y'], color=color_list[counter+1])
    if counter==0:
      spine_rois.append(roi)
    if counter ==2:
      counter = 0
      dend_rois.append(roi)
    else:
      counter+=1
  return spine_rois, dend_rois



################
#data out
#################

def convert_pixels_to_um(np_array):
    try:
        um = np_array * .09 #um/pixel
    except Exception as E:
        um = np.array(np_array) * .09 #um/pixel
    return um

def save_distances(spine_dmats, current_data_dir):
#what all do we want to save?
#save the image so we can view the ROIS
#save the 4 distance matricies as CSVs - make sure that they are in microns

  for key, pmat in spine_dmats.items():
    dmat = convert_pixels_to_um(pmat)

    file_name = str(key)+'.csv'
    file_dir = os.path.join(current_data_dir, subfolder_name)
    if not(os.path.isdir(file_dir)):
      os.mkdir(file_dir)
    file_path = os.path.join(current_data_dir, subfolder_name, file_name)
    #print(precision)
    np.savetxt(file_path, dmat, delimiter=",", fmt=precision)



def save_plot(spine_coords, dendrite_roi, current_data_dir):

    plt.figure()
    kyle_rois, projection_tif, shaft_roi = load_all(current_data_dir)
    plt.imshow(projection_tif, cmap='gray')

    color_list = ['b', 'g', 'r', 'c']
    label_list = ['annotated dendrite', 'annotated spine', 'annotated background', 'annotated shaft segment']
    counter = 0
    spine_rois = []
    dend_rois = []
    for name, roi in kyle_rois.items():
        plt.plot(roi['x'], roi['y'], color=color_list[counter+1])
        if counter==0:
            spine_rois.append(roi)
        if counter ==2:
            counter = 0
            dend_rois.append(roi)
        else:
            counter+=1

    for i, _ in enumerate(spine_coords['centers']):
        roi_center = spine_coords['centers'][i]
        spine_stem = spine_coords['spine_stems'][i]
        synapse = spine_coords['synapses'][i]
        #plot the center of the roi as an X
        plt.scatter(roi_center[0], roi_center[1], marker='x', color='y')

        #plot line from center to dendrite
        plt.plot([spine_stem[0], roi_center[0]], [spine_stem[1], roi_center[1]], color='y')
        #plot the location of the spine stem and the synaps
        plt.scatter(spine_stem[0], spine_stem[1], marker='o', color='y')
        plt.scatter(synapse[0], synapse[1], marker='+', color='y')

    den_xs = dendrite_roi['x']
    den_ys = dendrite_roi['y']
    #plot the dendrite
    plt.plot(den_xs, den_ys)

    #save it local with the other data
    file_name = 'annotated_dendrite.png'
    file_dir = os.path.join(current_data_dir, subfolder_name)
    if not(os.path.isdir(file_dir)):
        os.mkdir(file_dir)
    file_path = os.path.join(current_data_dir, subfolder_name, file_name)
    plt.savefig(file_path)

    if collect_images_at_path:
        cell_dir, FOV_name = os.path.split(current_data_dir)
        _, cell_name = os.path.split(cell_dir)
        file_name = cell_name+'_'+FOV_name+'_annotated_dendrite.png'
        if not(os.path.isdir(collect_images_at_path)):
            os.mkdir(collect_images_at_path)
        file_path = os.path.join(collect_images_at_path, file_name)
        plt.savefig(file_path)


def save_den_roi(dendrite_roi, current_data_dir):
    den_xs = dendrite_roi['x']
    den_ys = dendrite_roi['y']
    array = np.array([den_xs, den_ys])

    file_name = 'dendrite_roi.csv'
    file_dir = os.path.join(current_data_dir, subfolder_name)
    if not(os.path.isdir(file_dir)):
        os.mkdir(file_dir)
    file_path = os.path.join(current_data_dir, subfolder_name, file_name)
    np.savetxt(file_path, array, delimiter=",", fmt=precision)



def save_stem_stats(current_data_dir, **kwargs):
    DF = pd.DataFrame(kwargs)
    file_name = 'stem_stats.csv'
    file_dir = os.path.join(current_data_dir, subfolder_name)
    if not(os.path.isdir(file_dir)):
        os.mkdir(file_dir)
    file_path = os.path.join(current_data_dir, subfolder_name, file_name)

    DF.to_csv(file_path, index=False)


def save_summary_plots(**kwargs):
    for name, stat_dict in kwargs.items():
        save_summary_plot(stat_dict, name=name)


def save_summary_plot(stat_dict, name=''):
    plot_dict = {
        'all':[],
        'max':[],
        'median':[]
    }
    all_stat = []
    max_stat = []
    med_stat = []
    for key, stat_list in stat_dict.items():
        plot_dict['all'].extend(stat_list)
        plot_dict['max'].append(max(stat_list))
        plot_dict['median'].append(np.median(stat_list))

    fig, axs = plt.subplots(3, 1, figsize=(15,15))
    for i, (stat_type, stat_list) in enumerate(plot_dict.items()):
        axs[i].hist(stat_list)
        axs[i].set_xlabel('length')
        axs[i].set_ylabel('count')
        axs[i].set_title(name+' '+stat_type)

    file_name = name+'_summary_histograms.png'
    if not(os.path.isdir(collect_images_at_path)):
        os.mkdir(collect_images_at_path)
    file_path = os.path.join(collect_images_at_path, file_name)
    plt.savefig(file_path)



