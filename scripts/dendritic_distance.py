from src import helper_functions as hf
from src import config as cfg
from src import data_io as io
import os
import traceback
import logging
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)
logger.setLevel(logging.CRITICAL)

def main():
    global_sanity_checks = {}
    asymetric_list = []
    failed_list = []
    missing_annotation = []
    processed = []
    errors = {}
    for current_data_dir, dirs, files in os.walk(cfg.data_path, topdown=False):
        unprocessed = not (cfg.subfolder_name in dirs)
        if not (cfg.subfolder_name in current_data_dir):
            if cfg.re_run or unprocessed:
                print("Attempting to find dendritic distance in: " + str(current_data_dir))
                try:
                    kyle_rois, projection_tif, shaft_roi = io.load_all(current_data_dir)
                    if shaft_roi:
                        spine_rois, dend_rois = io.seperate_kyle_rois(kyle_rois)
                        all_segments, branch_points = hf.initialize_dendrite(shaft_roi)

                        # May want to put something lik this back in if we come up with a good way to find the dendrite
                        # dendrite_roi = {}
                        # dendrite_roi['t'], dendrite_roi['x'], dendrite_roi['y'], roi_center_xs, roi_center_ys = hf.infer_dendrite(dend_rois=dend_rois, shaft_roi=shaft_roi)
                        # spine_coords, spine_stem_t = hf.spine_stem_for_all_rois(spine_rois, dend_rois, dendrite_roi)

                        spines = []
                        sanity_checks = {}
                        for spine_roi, dend_roi in zip(spine_rois, dend_rois):
                            this_spine = hf.Spine(spine_roi, dend_roi, all_segments)
                            spines.append(this_spine)
                            for sanity_check, value in this_spine.sanity_checks.items():
                                if not (sanity_check in sanity_checks):
                                    sanity_checks[sanity_check] = []
                                sanity_checks[sanity_check].append(value)

                        spine_dmats = hf.euclidian_dmats(spines)
                        dend_d_mat, b_mat = hf.dendritic_distance_matrix(spines)

                        asymetry = False
                        for i in range(len(dend_d_mat)):
                            for j in range(len(dend_d_mat)):
                                if not (dend_d_mat[i, j] - dend_d_mat[j, i] < 0.0001):
                                    asymetry = current_data_dir
                        if asymetry:
                            asymetric_list.append(asymetry)

                        spine_dmats["dendritic_distance"] = dend_d_mat
                        spine_dmats["seperated_by_branch_point"] = b_mat.astype(int)

                        io.save_distances(spine_dmats, current_data_dir)
                        io.make_and_save_plot(current_data_dir, all_segments, spines, spine_dmats)

                        # io.save_den_roi(dendrite_roi, current_data_dir) this would need to account for branching
                        for check, value_list in sanity_checks.items():
                            if "length" in check:
                                sanity_checks[check] = io.convert_pixels_to_um(value_list)
                            if "area" in check:
                                sanity_checks[check] = io.convert_area_to_um(value_list)
                            if not (check in global_sanity_checks):
                                global_sanity_checks[check] = {}
                            # putting in the current dict so we can track down errors after the fact
                            global_sanity_checks[check][current_data_dir] = sanity_checks[
                                check
                            ]

                        io.save_stem_stats(current_data_dir, **sanity_checks)
                        processed.append(current_data_dir)
                    else:
                        missing_annotation.append(current_data_dir)
                except Exception as E:
                    failed_list.append(current_data_dir)
                    err_str = f"There was an error processing directory: {current_data_dir}"
                    logger.error(err_str)
                    logger.warning(traceback.format_exc())
                    errors[current_data_dir] = traceback.format_exc()
                    #raise(E)
    io.save_summary_plots(**global_sanity_checks)
    io.save_named_iterable_to_json(asymetric_dmats_dir_list=asymetric_list)
    io.save_named_iterable_to_json(failed_dirs_list=failed_list)
    io.save_named_iterable_to_json(failed_dirs_errors=errors)
    io.save_named_iterable_to_json(dirs_missing_annotation=missing_annotation)
    io.save_named_iterable_to_json(sucessfully_processed_dirs=processed)


if __name__ == "__main__":
    main()
