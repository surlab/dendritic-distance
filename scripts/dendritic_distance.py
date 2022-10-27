from src import helper_functions as hf
from src import config as cfg
from src import data_io as io
import os
import traceback
import logging
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
#logger.setLevel(logging.ERROR)
#logger.setLevel(logging.CRITICAL)

def main():
    globals = {
        'global_sanity_checks': {},
        'asymetric_list': [],
        'failed_list': [],
        'missing_annotation': [],
        'processed': [],
        'errors': {},
    }

    if cfg.walk_dirs:
        for current_data_dir, dirs, files in os.walk(cfg.data_path, topdown=False):
            unprocessed = not (cfg.subfolder_name in dirs)
            if not (cfg.subfolder_name in current_data_dir):
                if cfg.re_run or unprocessed:
                    main_2(current_data_dir, globals)
    else:
        main_2(cfg.data_path, globals)

    io.save_summary_plots(**globals['global_sanity_checks'])
    io.save_named_iterable_to_json(asymetric_dmats_dir_list=globals['asymetric_list'])
    io.save_named_iterable_to_json(failed_dirs_list=globals['failed_list'])
    io.save_named_iterable_to_json(failed_dirs_errors=globals['errors'])
    io.save_named_iterable_to_json(dirs_missing_annotation=globals['missing_annotation'])
    io.save_named_iterable_to_json(sucessfully_processed_dirs=globals['processed'])



def main_2(current_data_dir, globals):
                    print("Attempting to find dendritic distance in: " + str(current_data_dir))
                    try:
                        kyle_rois, projection_tif, shaft_roi = io.load_all(current_data_dir)
                        if shaft_roi:
                            spine_rois, dend_rois, other_rois = io.seperate_kyle_rois(kyle_rois)
                            all_segments, branch_points, fiducial_rois = hf.initialize_dendrite(shaft_roi)

                            # May want to put something lik this back in if we come up with a good way to find the dendrite
                            # dendrite_roi = {}
                            # dendrite_roi['t'], dendrite_roi['x'], dendrite_roi['y'], roi_center_xs, roi_center_ys = hf.infer_dendrite(dend_rois=dend_rois, shaft_roi=shaft_roi)
                            # spine_coords, spine_stem_t = hf.spine_stem_for_all_rois(spine_rois, dend_rois, dendrite_roi)

                            spines = []
                            sanity_checks = {}
                            spine_stats = {}
                            spine_mapping = {}
                            spine_mapping['source_file'] = []
                            spine_mapping['roi_name'] = []
                            for spine_roi, dend_roi in zip(spine_rois, dend_rois):
                                this_spine = hf.Spine(spine_roi, dend_roi, all_segments)
                                spines.append(this_spine)
                                for sanity_check, value in this_spine.sanity_checks.items():
                                    if not (sanity_check in sanity_checks):
                                        sanity_checks[sanity_check] = []
                                    sanity_checks[sanity_check].append(value)
                                for spine_stat, value in this_spine.spine_stats.items():
                                    if not (spine_stat in spine_stats):
                                        spine_stats[spine_stat] = []
                                    spine_stats[spine_stat].append(value)
                                spine_mapping['source_file'].append(this_spine.source_file)
                                spine_mapping['roi_name'].append(this_spine.name)


                            fiducials = []
                            for name, fiducial_roi in fiducial_rois.items():
                                this_fid = hf.Spine(fiducial_roi, fiducial_roi, all_segments)
                                fiducials.append(this_fid)

                            spine_dmats = hf.euclidian_dmats(spines)
                            dend_d_mat, b_mat = hf.dendritic_distance_matrix(spines)
                            spine_dmats["dendritic_distance"] = dend_d_mat
                            spine_dmats["seperated_by_branch_point"] = b_mat.astype(int)

                            fiducial_dmat, b_mat = hf.fiducial_dend_mats(spines, fiducials)
                            spine_dmats["fiducial_distances"] = fiducial_dmat

                            io.save_distances(spine_dmats, current_data_dir)
                            io.make_and_save_plot(current_data_dir, all_segments, spines, spine_dmats)


                            asymetry = False
                            for i in range(len(dend_d_mat)):
                                for j in range(len(dend_d_mat)):
                                    if not (dend_d_mat[i, j] - dend_d_mat[j, i] < 0.0001):
                                        asymetry = current_data_dir
                            if asymetry:
                                globals['asymetric_list'].append(asymetry)

                            # io.save_den_roi(dendrite_roi, current_data_dir) this would need to account for branching
                            def convert_stat_dict_for_output(stat_dict):
                                new_stat_dict = {}
                                for check, value_list in stat_dict.items():
                                    if "length" in check:
                                        new_stat_dict[check] = io.convert_pixels_to_um(value_list)
                                    if "area" in check:
                                        new_stat_dict[check] = io.convert_area_to_um(value_list)
                                return new_stat_dict

                            sanity_checks = convert_stat_dict_for_output(sanity_checks)
                            for check, value_list in sanity_checks.items():
                                if not (check in globals['global_sanity_checks']):
                                    globals['global_sanity_checks'][check] = {}
                                # putting in the current dict so we can track down errors after the fact
                                globals['global_sanity_checks'][check][current_data_dir] = sanity_checks[
                                    check
                                ]

                            stat_dict = {**spine_stats, **sanity_checks, **spine_mapping}

                            io.save_stats(current_data_dir, "spine_stats.csv", **stat_dict)
                            #io.save_stem_stats(current_data_dir, "spine_mapping.csv" ,**spine_mapping)
                            globals['processed'].append(current_data_dir)
                        else:
                            globals['missing_annotation'].append(current_data_dir)
                    except Exception as E:
                        globals['failed_list'].append(current_data_dir)
                        err_str = f"There was an error processing directory: {current_data_dir}"
                        logger.error(err_str)
                        logger.warning(traceback.format_exc())
                        globals['errors'][current_data_dir] = traceback.format_exc()
                        #raise(E)


if __name__ == "__main__":
    main()
