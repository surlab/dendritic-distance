import helper_functions as hf
from config import data_path
import data_io as io
import os



def main():
    stem_length_dict = {}
    stem_angle_dict = {}
    dendrite_residuals_dict = {}
    for current_data_dir, dirs, files in os.walk(data_path, topdown=False):
        print('Attempting to find dendritic distance in: '+str(current_data_dir))


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

            stem_lengths = io.convert_pixels_to_um(hf.stem_length(spine_coords))
            stem_angles = hf.stem_angle(spine_coords, spine_stem_t, dendrite_roi)
            dendrite_residuals = io.convert_pixels_to_um(hf.dendrite_residual(dend_rois, dendrite_roi))
            io.save_stem_stats(current_data_dir, stem_lengths=stem_lengths, stem_angles=stem_angles, dendrite_residuals=dendrite_residuals)
            stem_length_dict[current_data_dir] = stem_lengths
            stem_angle_dict[current_data_dir] = stem_angles
            dendrite_residuals_dict[current_data_dir] = dendrite_residuals
        except Exception as E:
            print(E)
            #raise(E)
    io.save_summary_plots(lengths = stem_length_dict, angles = stem_angle_dict, dendrite_residuals=dendrite_residuals_dict)


if __name__=='__main__':
    main()
