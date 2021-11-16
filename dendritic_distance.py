import helper_functions as hf
from config import data_path
import data_io as io
import os



def main():
    for current_data_dir, dirs, files in os.walk(data_path, topdown=False):
        try:
            kyle_rois, projection_tif, shaft_roi = io.load_all(current_data_dir)

            spine_rois, dend_rois = io.seperate_kyle_rois(kyle_rois)

            dendrite_roi = {}
            dendrite_roi['t'], dendrite_roi['x'], dendrite_roi['y'], roi_center_xs, roi_center_ys = hf.infer_dendrite(dend_rois=dend_rois, shaft_roi=shaft_roi)
            spine_coords, spine_stem_t = hf.spine_stem_for_all_rois(spine_rois, dendrite_roi)

            spine_dmats = hf.euclidian_dmats(spine_coords)
            dend_d_mat = hf.dendritic_distance(spine_stem_t, dendrite_roi)

            spine_dmats['dendritic_distance'] = dend_d_mat

            io.save_distances(spine_dmats, current_data_dir)
            io.save_plot(spine_coords, dendrite_roi, current_data_dir)

            io.save_den_roi(dendrite_roi, current_data_dir)
        except Exception as E:
            pass

if __name__=='__main__':
    main()
