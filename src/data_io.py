
import logging

try:
    from src import config as cfg
except ImportError as E:
    logging.warning('No custom config found, using default config instead')
    from src import default_config as cfg


# should probably import this whole thing and use the namespace...
from read_roi import read_roi_file, read_roi_zip
import os
import PIL
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


import json
from datetime import datetime as dt

################
# data in
#################

def readfile(path):
    print(f"Reading file at {path}")
    with open(path, "r") as f:
        lines = f.read()
        print(lines)

### for kyle
def get_paths_from_data_path(data_path):
    return cfg.get_paths_from_data_path(data_path)


def load_all(data_path):

    rois_zip_paths, projection_tif_path, shaft_roi_path = get_paths_from_data_path(
        data_path
    )
   # print(rois_zip_path, projection_tif_path, shaft_roi_path)
    kyle_rois = {}
    def load_roi(rois_zip_path):
        roi_dict = _read_roi(rois_zip_path)
        roi_filename = os.path.split(rois_zip_path)[1]
        roi_dict = add_source_file(roi_dict, roi_filename)
        return roi_dict

    def add_source_file(roi_dict, roi_filename):
        new_roi_dict = {}
        for key, value in roi_dict.items():
            new_key = f'{key}_{roi_filename}'
            new_roi_dict[new_key] = value
            new_roi_dict[new_key]['source_file'] = roi_filename

        return new_roi_dict

    for rois_zip_path in rois_zip_paths:
        kyle_rois.update(load_roi(rois_zip_path))

    projection_tif = PIL.Image.open(projection_tif_path)
    shaft_roi = None
    if shaft_roi_path:
        shaft_roi = load_roi(shaft_roi_path)

    return kyle_rois, projection_tif, shaft_roi


def _read_roi(roi_path):
    try:
        roi = read_roi_file(roi_path)
    except Exception as E:
        roi = read_roi_zip(roi_path)
    return roi


def seperate_kyle_rois(kyle_rois):
    spine_rois = []
    dend_rois = []
    other_rois = []
    counter = 0
    for name, roi in kyle_rois.items():
        # plt.plot(roi['x'], roi['y'], color=color_list[counter+1])
        if 'spine' in cfg.roi_order[counter].lower():
            spine_rois.append(roi)
        elif 'dendrite' in cfg.roi_order[counter].lower():
            dend_rois.append(roi)
        else:
            other_rois.append(roi)
        if counter == len(cfg.roi_order)-1:
            counter = 0
        else:
            counter += 1
    return spine_rois, dend_rois, other_rois


################
# data out
#################


def convert_pixels_to_um(np_array):
    try:
        um = np_array * cfg.conversion_factor  # um/pixel
    except Exception as E:
        um = np.array(np_array) * cfg.conversion_factor  # um/pixel
    return um


def convert_area_to_um(np_array):
    return convert_pixels_to_um(convert_pixels_to_um(np_array))


def save_distances(spine_dmats, current_data_dir):
    # what all do we want to save?
    # save the image so we can view the ROIS
    # save the 4 distance matricies as CSVs - make sure that they are in microns

    # don't convert the binary branch point matrix
    for key, pmat in spine_dmats.items():
        if "distance" in key.lower():
            dmat = convert_pixels_to_um(pmat)
        else:
            dmat = pmat

        file_name = str(key) + ".csv"
        file_dir = os.path.join(current_data_dir, cfg.subfolder_name)
        if not (os.path.isdir(file_dir)):
            os.mkdir(file_dir)
        file_path = os.path.join(current_data_dir, cfg.subfolder_name, file_name)
        np.savetxt(file_path, dmat, delimiter=",", fmt=cfg.precision)

def new_ax(ax):
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = None
    return fig, ax


def plot_dendrite(all_segments, ax=None):
    fig, ax = new_ax(ax)
    for segment_key, segment in all_segments.items():
        ax.plot(segment.dend_xs, segment.dend_ys, color="b")
        for endpoint in segment.end_points:
            # mark the endpoints and branch points
            m = "o" if endpoint.branch_point else "x"
            ax.scatter(
                endpoint.end_point_xy[0], endpoint.end_point_xy[1], color="b", marker=m
            )
    return fig, ax

def plot_dendrite_rois(dend_rois, ax=None, legend=False, invert=False):
    fig, ax = new_ax(ax)
    if invert:
        invert = -1
    else:
        invert=1
    color_list = ['b', 'm']
    label_list = ['annotated shaft segment', 'annotated_branch_point']
    counter = 0
    for name, roi in dend_rois.items():
        if "line" in roi["type"]:
            color= color_list[0]
            ax.plot(np.array(roi['x']), invert*np.array(roi['y']), color=color)
        elif "oval" in roi["type"]:
            color= color_list[1]
            #x = roi["left"] + roi["width"] / 2
            #y = roi["top"] + roi["height"] / 2
            roi = convert_roi_to_polygon(roi)
            ax.scatter(roi['center'][0],invert*roi['center'][1], marker="o", color=color) #TODO change the size of this/actually plot the oval
        else:
            roi = convert_roi_to_polygon(roi)
            color = '0'
            ax.plot(roi['center'][0],invert*roi['center'][1], color=color)
    ax.set_aspect("equal", adjustable="box")
    if legend:
        custom_legend(ax, color_list, label_list)
    return fig, ax

def convert_roi_to_polygon(roi):
    if "x" in roi and 'y' in roi:
        roi['center'] = np.array((np.mean(roi["x"]), np.mean(roi["y"])))
        #segment_data["x"], segment_data["y"]
    elif 'oval' in roi['type']:
        #for key, value in roi.items():
        #    print(key, value)
        x = roi["left"] + roi["width"] / 2
        y = roi["top"] + roi["height"] / 2
        roi['center'] = np.array((x, y))
    elif 'freehand' in roi['type']:
        try:
            p1 = np.array([roi['ex1'], roi['ey1']])
            p2 = np.array([roi['ex2'], roi['ey2']])

            center = np.mean([p1,p2], axis=0)
            vect = p2-p1
            distance = np.linalg.norm(vect)
            def get_orthogonal(vect):
                x = np.random.randn(len(vect))  # take a random vector
                x -= x.dot(vect) * vect / np.linalg.norm(vect)**2      # make it orthogonal to k
                x /= np.linalg.norm(x)  # normalize it
                return x
            orth_vect_norm = get_orthogonal(vect)
            orth_vect = orth_vect_norm*(distance*roi['aspect_ratio'])/2
            p3 = center+orth_vect
            p4 = center-orth_vect
            #we need a way to oder these properly, but I'm gonna forge ahead for now

            roi['x'] = [p1[0], p3[0], p2[0], p4[0], p1[0]]
            roi['y'] = [p1[1], p3[1], p2[1], p4[1], p1[1]]
            roi['center'] = np.array(center)
            roi['type'] = 'freehand_oval'
        except KeyError as E:
            for key, value in roi.items():
                print('Here:', key, value)
    elif 'line' in roi['type']:
        #this will only grab pure lines that don't already have x,y coordinates, just x1,y1,x2,y2, nedpoints
        p1 = np.array([roi['x1'], roi['y1']])
        p2 = np.array([roi['x2'], roi['y2']])
        center = np.mean([p1,p2], axis=0)
        roi['x'] = [p1[0], center[0], p2[0]]
        roi['y'] = [p1[1], center[1], p2[1]]
        roi['center'] = np.array(center)
        roi['type'] = 'fiducial'
    else:
        for key, value in roi.items():
            print('Or here: ', key, value)

    return roi

def custom_legend(ax, color_list, label_list):
    legend_elements = [
        Line2D([0], [0], color=color_list[idx], lw=2, label=label_list[idx]) for idx in range(len(label_list))
        ]
    ax.legend(handles=legend_elements)
    return ax

def projection_tif_from_path(filepath, ax=None):
    im = PIL.Image.open(filepath)
    return projection_tif(im, ax)

def projection_tif(im, ax=None):
    fig, ax = new_ax(ax)
    ax.imshow(im, cmap='gray')
    return fig, ax


def plot_kyle_rois(kyle_rois, ax=None, legend=False, invert=False):
    fig, ax = new_ax(ax)
    if invert:
        invert = -1
    else:
        invert=1
    color_list = ['g', 'r', 'c']
    label_list = ['annotated spine', 'annotated background', 'annotated dendrite']
    split_rois = seperate_kyle_rois(kyle_rois)
    counter = 0
    for idx, roi_list in enumerate(split_rois):
        for roi in roi_list:
            roi = convert_roi_to_polygon(roi)
            ax.plot(np.array(roi['x']), invert*np.array(roi['y']), color=color_list[idx])

    ax.set_aspect("equal", adjustable="box")
    if legend:
        custom_legend(ax, color_list, label_list)
    return fig, ax


def plot_spines(spines, ax=None):
    fig, ax = new_ax(ax)
    for spine in spines:
        # plot the ROI in question
        ax.plot(spine.spine_roi[0, :], spine.spine_roi[1, :], color="y")
        # plot the center of the roi as an X
        ax.scatter(
            spine.spine_center_xy[0], spine.spine_center_xy[1], marker="x", color="y"
        )
        # plot line from center to dendrite
        if spine.relative_neck_angle > 0:
            color = 'y'
        else:
            color='chartreuse'
        ax.plot(
            [spine.spine_neck_xy[0], spine.spine_center_xy[0]],
            [spine.spine_neck_xy[1], spine.spine_center_xy[1]],
            color=color,
        )
        # plot the location of the spine stem and the synapse
        ax.scatter(
            spine.spine_neck_xy[0], spine.spine_neck_xy[1], marker="o", color="y"
        )
        ax.scatter(spine.synapse_xy[0], spine.synapse_xy[1], marker="+", color="y")
    return fig, ax


def plot_example_distances(dmats, spines, ax=None):
    fig, ax = new_ax(ax)
    for key, dmat in dmats.items():
        if "center" in key.lower() and "euclidian" in key.lower():
            center_ds = dmat
        if "dendritic" in key:
            dendritic_ds = dmat
    diff = center_ds - dendritic_ds

    def plot_connection(i,j, dendritic_ds, center_ds, ax, c='0'):
        dend_label = "dendritic distance = " + str(
            convert_pixels_to_um(dendritic_ds[i, j])
        )
        eu_label = "euclidian distance = " + str(convert_pixels_to_um(center_ds[i, j]))
        plt.plot(
            [spines[i].spine_center_xy[0], spines[j].spine_center_xy[0]],
            [spines[i].spine_center_xy[1], spines[j].spine_center_xy[1]],
            color="c",
            label=dend_label,
        )
        plt.plot(
            [spines[i].spine_center_xy[0], spines[j].spine_center_xy[0]],
            [spines[i].spine_center_xy[1], spines[j].spine_center_xy[1]],
            color=c,
            label=eu_label,
        )
        return ax

    i, j = np.unravel_index(np.argmax(diff), diff.shape)
    plot_connection(i,j, dendritic_ds, center_ds, ax, c="c")

    i, j = np.unravel_index(np.argmin(diff), diff.shape)
    plot_connection(i,j, dendritic_ds, center_ds, ax, c="m")

    plt.legend()
    return fig, ax

def overlay_plot(tiff_in, all_segments, spines, dmats, ax=None):
    fig, ax = new_ax(ax)

    projection_tif(tiff_in, ax)
    plot_dendrite(all_segments, ax)
    plot_spines(spines, ax)
    plot_example_distances(dmats, spines, ax)
    return fig, ax


def make_and_save_plot(current_data_dir, all_segments, spines, dmats):
    kyle_rois, projection_tif, shaft_roi = load_all(current_data_dir)
    fig, ax=overlay_plot(projection_tif, all_segments, spines, dmats)

    save_plot(fig, current_data_dir)
    return fig, ax

def save_plot(fig, current_data_dir):

    # save it local with the other data
    file_name = "annotated_dendrite.svg"
    file_dir = os.path.join(current_data_dir, cfg.subfolder_name)
    if not (os.path.isdir(file_dir)):
        os.mkdir(file_dir)
    file_path = os.path.join(current_data_dir, cfg.subfolder_name, file_name)
    fig.savefig(file_path)

    if cfg.collect_summary_at_path:
        cell_dir, FOV_name = os.path.split(current_data_dir)
        _, cell_name = os.path.split(cell_dir)
        file_name = cell_name + "_" + FOV_name + "_annotated_dendrite.svg"
        if not (os.path.isdir(cfg.collect_summary_at_path)):
            os.mkdir(cfg.collect_summary_at_path)
        file_path = os.path.join(cfg.collect_summary_at_path, file_name)
        fig.savefig(file_path)


def save_den_roi(dendrite_roi, current_data_dir):
    den_xs = dendrite_roi["x"]
    den_ys = dendrite_roi["y"]
    array = np.array([den_xs, den_ys])

    file_name = "dendrite_roi.csv"
    file_dir = os.path.join(current_data_dir, cfg.subfolder_name)
    if not (os.path.isdir(file_dir)):
        os.mkdir(file_dir)
    file_path = os.path.join(current_data_dir, cfg.subfolder_name, file_name)
    np.savetxt(file_path, array, delimiter=",", fmt=cfg.precision)


def save_stats(current_data_dir, file_name, **kwargs):
    DF = pd.DataFrame(kwargs)
    file_dir = os.path.join(current_data_dir, cfg.subfolder_name)
    if not (os.path.isdir(file_dir)):
        os.mkdir(file_dir)
    file_path = os.path.join(current_data_dir, cfg.subfolder_name, file_name)

    DF.to_csv(file_path, index=False)


def save_summary_plots(**kwargs):
    for name, stat_dict in kwargs.items():
        save_summary_plot(stat_dict, name=name)


def save_summary_plot(stat_dict, name=""):
    plot_dict = {"all": [], "max": [], "median": []}


    max_dict = {key: max(stat_list) for key, stat_list in stat_dict.items()}
    def sort_by_value(pair):
        return pair[1]

    ordered_stat_list = sorted(max_dict.items(), key=sort_by_value, reverse=True)

    num = 5
    max_all = ([None] * num, [-np.inf] * num)
    max_med = ([None] * num, [-np.inf] * num)
    max_max = ([None] * num, [-np.inf] * num)
    min_all = ([None] * num, [np.inf] * num)
    min_med = ([None] * num, [np.inf] * num)
    min_max = ([None] * num, [np.inf] * num)
    for key, stat_list in stat_dict.items():

        # may need to pull out NANs here, but there really shouldn't be any... coming from the angle calculation somehow
        stat_list = np.array(stat_list)
        stat_list = stat_list[~np.isnan(stat_list)]
        plot_dict["all"].extend(stat_list)
        value = min(stat_list)
        if value < max(min_all[1]):
            idx = np.argmax(min_all[1])
            min_all[1][idx] = value
            min_all[0][idx] = key

        value = max(stat_list)
        plot_dict["max"].append(value)
        if value > min(max_all[1]):
            idx = np.argmin(max_all[1])
            max_all[1][idx] = value
            max_all[0][idx] = key
            max_max[1][idx] = value
            max_max[0][idx] = key

        if value < max(min_max[1]):
            idx = np.argmax(min_max[1])
            min_max[1][idx] = value
            min_max[0][idx] = key

        value = np.median(stat_list)
        plot_dict["median"].append(value)
        if value > min(max_med[1]):
            idx = np.argmin(max_med[1])
            max_med[1][idx] = value
            max_med[0][idx] = key

        if value < max(min_med[1]):
            idx = np.argmax(min_med[1])
            min_med[1][idx] = value
            min_med[0][idx] = key

    fig, axs = plt.subplots(3, 1, figsize=(15, 15))
    for i, (stat_type, stat_list) in enumerate(plot_dict.items()):
        axs[i].hist(stat_list)
        axs[i].set_xlabel("length")
        axs[i].set_ylabel("count")
        axs[i].set_title(name + " " + stat_type)

    timestamp = dt.now().strftime("%Y%m%d_%H%M%S")
    file_name = name + "_summary_histograms_" + timestamp + ".png"
    summary_path = os.path.join(cfg.collect_summary_at_path, "summary_plots")
    if not (os.path.isdir(summary_path)):
        os.makedirs(summary_path)
    file_path = os.path.join(summary_path, file_name)
    plt.savefig(file_path)

    check_path_dict = {
        "Ordered highest to lowest": ordered_stat_list,
        "max_all": max_all,
        "max_med": max_med,
        "max_max": max_max,
        "min_all": min_all,
        "min_med": min_med,
        "min_max": min_max,
    }
    file_name = name + "_check_directories_" + timestamp + ".json"
    file_path = os.path.join(summary_path, file_name)
    with open(file_path, "w") as f:
        json.dump(check_path_dict, f, indent=4)


def save_named_iterable_to_json(**kwargs):
    for key, value in kwargs.items():
        timestamp = dt.now().strftime("%Y%m%d_%H%M%S")
        file_name = key +'_' + timestamp + ".json"
        summary_path = os.path.join(cfg.collect_summary_at_path, "summary_plots")
        file_path = os.path.join(summary_path, file_name)
        if not (os.path.isdir(summary_path)):
            os.makedirs(summary_path)
        with open(file_path, "w") as f:
            json.dump(value, f, indent=4)
