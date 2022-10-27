import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import distance_matrix
import math

from src import data_io as io

def endpoints(xs, ys):
    start_coords = np.array((xs[0], ys[0]))
    end_coords = np.array((xs[-1], ys[-1]))
    return start_coords, end_coords


def roi_to_array(xs, ys):
    array = np.vstack((xs, ys))
    return array


# ^^ this shouldn't really be a function...

# this should actually be in Data_io I think, it isn't actually doing any manipulation
def initialize_dendrite(shaft_rois):

    unassigned_segments = list(shaft_rois.keys())
    all_segments = {}
    unassigned_branch_regions = {}
    fiducials = {}
    for segment_key, segment_data in shaft_rois.items():
        # initialize all the segments here, then we add the branch points later (since we can't possibly know the map beforehand)
        segment_data = io.convert_roi_to_polygon(segment_data)
        if "line" in segment_data["type"]:
            all_segments[segment_key] = Dendrite_segment(segment_key, segment_data)
            start_coords, end_coords = endpoints(segment_data["x"], segment_data["y"])
        elif "oval" in segment_data["type"]:
            unassigned_branch_regions[segment_key] = segment_data
            fiducials[segment_key] = segment_data
        else:
            fiducials[segment_key] = segment_data


    def find_branch_point(unassigned_branch_regions, segment_objects):
        # pick the starting segment
        next_branch_key = next(iter(unassigned_branch_regions))
        next_branch_region = unassigned_branch_regions.pop(next_branch_key)

        #want to make this compatible with rectangles and freehand ovals... should just give them all Xs and ys
        #x = next_branch_region["left"] + next_branch_region["width"] / 2
        #y = next_branch_region["top"] + next_branch_region["height"] / 2
        #branch_center = np.array((x, y))
        marker_center = next_branch_region['center']


        dists = []
        endpoint_mapping = []
        segment_mapping = []
        for segment_key, segment_object in segment_objects.items():
            for coords in segment_object.end_points_coords:
                dists.append(np.linalg.norm(marker_center - coords, axis=0))
                endpoint_mapping.append(coords)
                segment_mapping.append(segment_key)

        connected_enpoints = []
        connected_segments = []
        # look for 3 connected branches - 3 nearest endpoints
        dists = np.array(dists)
        num_segments = min(3, len(all_segments))
        for times in range(num_segments):
            # find nearest
            min_idx = np.argmin(dists)
            # add it to the list
            connected_enpoints.append(endpoint_mapping[min_idx])
            segment = segment_mapping[min_idx]
            assert not (segment in connected_segments)
            connected_segments.append(segment)
            # set the values for the attached segment to a high value so it won't be picked again (same segment should never have both endpounts at a branch - no loops)
            dists[np.array(segment_mapping) == segment] = max(dists)

        branch_point_xy = np.mean(np.array(connected_enpoints), axis=0)
        child_dendrite_list = []
        for segment_key in connected_segments:
            child_dendrite_list.append(segment_objects[segment_key])

        branch_point = End_point(branch_point_xy, child_dendrite_list)

        for segment_object in child_dendrite_list:
            segment_object.assign_endpoint_mapping(branch_point)

        return branch_point

    branch_points = []
    # keep going until all the segments are assigned to at least 1 cluster
    while unassigned_branch_regions:
        branch_points.append(find_branch_point(unassigned_branch_regions, all_segments))

    # assign endpoint objects to the remaining endpoints
    for segment_key, segment in all_segments.items():
        for i, endpoint in enumerate(segment.end_points):
            if endpoint is None:

                end_point = End_point(segment.end_points_coords[i].T, [segment])
                segment.assign_endpoint_mapping(end_point)

    return all_segments, branch_points, fiducials


def _m_distance(point, line):

    num_elements = max(np.shape(line))
    points = np.tile(point, (num_elements, 1))
    # having trouble dealing with arrays going in different directions... so I put in this try
    try:
        diff = line - points  # , (1, num_elements))#.reshape(2,num_elements)
    except ValueError as E:
        diff = line.T - points
    res_sq = np.linalg.norm(diff, axis=1)
    return res_sq


def minimum_distance(point, line):
    res_sq = _m_distance(point, line)
    # find index of minimum distance
    min_dist = np.min(res_sq)
    min_dist_idx = np.argmin(res_sq)
    return min_dist, min_dist_idx


def distance_along_curve(t_1, t_2, xs, ys):
    # should this just take an array is input instead of the roi style dict?
    dist = 0
    start_t = min(t_1, t_2)
    end_t = max(t_1, t_2)
    assert end_t <= len(xs)

    # this is where we could introduce direction
    for i in range(start_t, end_t - 1):
        x_1 = xs[i]
        x_2 = xs[i + 1]
        y_1 = ys[i]
        y_2 = ys[i + 1]
        dist += np.sqrt((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2)
    return dist


class End_point:
    def __init__(self, end_point_xy, child_dendrite_list):
        self.end_point_xy = end_point_xy
        self.child_dendrites = child_dendrite_list

    @property
    def branch_point(self):
        return len(self.child_dendrites) > 2


def dense_interpolate(xs, ys):
    # assuming the coords are in the right order - could be necessary if dendrite loops back near itself
    dense_xs = []
    dense_ys = []
    total_points = (
        1000  # could parameterize this - maybe 2-3x the pixel dimensions of the image?
    )
    # doing linear interpolation since other methods can be problematic in certain cases like sharp angles. if we are manually annotating this should be fine
    # want to make sure that the density is even along the whole dendrite, not dependant on the density of the original manually annotated points
    approx_length = distance_along_curve(0, len(xs), xs, ys)
    points_per_distance = total_points / approx_length
    for i in range(len(xs) - 1):
        seg_xs = xs[i : i + 2]
        seg_ys = ys[i : i + 2]
        segment_length = distance_along_curve(0, len(seg_xs), seg_xs, seg_ys)
        num_points_to_add = round(points_per_distance * segment_length)
        dense_segment_xs = list(np.linspace(seg_xs[0], seg_xs[-1], num_points_to_add))
        dense_xs.extend(dense_segment_xs)
        dense_ys.extend(list(np.linspace(seg_ys[0], seg_ys[-1], num_points_to_add)))
    dense_ts = np.arange(0, len(dense_xs), 1)
    return dense_ts, dense_xs, dense_ys


class Dendrite_segment:
    def __init__(self, segment_key, segment_data):
        self.key = segment_key
        # do the interpolation here to create points all along the annotated line
        self.dend_ts, self.dend_xs, self.dend_ys = dense_interpolate(
            segment_data["x"], segment_data["y"]
        )

        self.coords = roi_to_array(self.dend_xs, self.dend_ys)

        self.end_points_coords = endpoints(segment_data["x"], segment_data["y"])
        self.end_points = [None, None]

        self.segment_length = distance_along_curve(
            self.dend_ts[0], self.dend_ts[-1], self.dend_xs, self.dend_ys
        )
        # ^this could lead to slight errors if the manually annotated segment overshoots the calcullated branch point a bit
        # this will only matter if its the middle segment which means the spines are already a full segment away from eachother so the slight miscalculation shouldn't matter
        # and should be within error from other sources of innacuracy
        # although we could update this every time there is a new enpoint assignment

    def assign_endpoint_mapping(self, segment_endpoint):
        # compute the closest endpoint
        dist1 = np.linalg.norm(
            self.end_points_coords[0] - segment_endpoint.end_point_xy
        )
        dist2 = np.linalg.norm(
            self.end_points_coords[1] - segment_endpoint.end_point_xy
        )

        # make sure its reasonably close
        # assert(min((dist1, dist2)) < self.segment_length/8)
        # ^this doesn't work because sometimes I used very short segment lengths.

        # assign the branch point/endpoint to the appropriate slot
        if dist1 < dist2:
            self.end_points[0] = segment_endpoint
        else:
            self.end_points[1] = segment_endpoint
        # TODO could adjust the segment length measurement here
        # technically we shouldn't just redo the distance calculation - we should rewrite the coords and interpolation so they go to the branch point appropriately...


class Spine:
    def __init__(self, spine_roi, spine_dendrite_roi, dendrite_segment_rois):

        spine_roi = io.convert_roi_to_polygon(spine_roi)
        spine_dendrite_roi = io.convert_roi_to_polygon(spine_dendrite_roi)

        self.parent_dendritic_segment, self.dendrite_residual = find_dendritic_segment(
            spine_dendrite_roi, dendrite_segment_rois
        )

        self.source_file = spine_roi['source_file']
        self.name = spine_roi['name']


        self.spine_roi = roi_to_array(spine_roi["x"], spine_roi["y"])
        self.spine_center_xy = np.mean(self.spine_roi, axis=1)

        #if cfg.neck_source=='center':
        min_dist, self.spine_neck_idx = minimum_distance(
            self.spine_center_xy, self.parent_dendritic_segment.coords
        )
        #elif cfg.neck_source=='nearest':
        # <<<<<<<<<<<<<<<<<<
        # maybe this should be minimizing the pairwise distance on both end?
        # although this might introduce as many errors as it fixes... need to look at examples
        # grab the appropriate location on the dendrite
        self.spine_neck_xy = self.parent_dendritic_segment.coords[
            :, self.spine_neck_idx
        ]

        # find the furthest distance from the neck as the synapse (same steps as above)
        max_dist, max_dist_idx = maximum_distance(self.spine_neck_xy, self.spine_roi)
        self.synapse_xy = self.spine_roi[:, max_dist_idx]

        self.neck_length = np.sqrt(
            (self.spine_neck_xy[0] - self.spine_center_xy[0]) ** 2
            + (self.spine_neck_xy[1] - self.spine_center_xy[1]) ** 2
        )

        # subtract the spine neck coordinates to get vectors from the origin for directions to the other 2
        spine_vect = self.spine_center_xy - self.spine_neck_xy
        try:
            den_vect = (
                self.parent_dendritic_segment.coords[:, self.spine_neck_idx + 2]
                - self.spine_neck_xy
            )
        except Exception as E:
            den_vect = (
                self.parent_dendritic_segment.coords[:, self.spine_neck_idx - 2]
                - self.spine_neck_xy
            )

        # taken from https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
        #will need to figure out how to do the signed version ... need to be taking a consistent direction from the
        #starting fiducial, and also need to be clever how to adapt it to 3d.
        #dot = np.dot(den_vect, spine_vect)
        #den_len = np.linalg.norm(den_vect) + 0.001
        #spine_len = np.linalg.norm(spine_vect) + 0.001
        #angle = dot / den_len / spine_len  # -> cosine of the angle
        #self.relative_neck_angle = np.arccos(np.clip(angle, -1, 1))
        #flat_vect = [[0,0],[1,0]]
        #dot = np.dot(flat_vect, spine_vect)
        #flat_len = np.linalg.norm(flat_vect) + 0.0001
        #angle = dot / flat_len / spine_len  # -> cosine of the angle
        #self.global_neck_angle = np.arccos(np.clip(angle, -1, 1))

        #from https://stackoverflow.com/questions/2150050/finding-signed-angle-between-vectors
        def get_signed_angle_2d(a, b):
            angle = math.atan2( a[0]*b[1] - a[1]*b[0], a[0]*b[0] + a[1]*b[1])
            return angle
        self.relative_neck_angle = get_signed_angle_2d(spine_vect, den_vect)
        flat_vect = np.array([1,0])
        self.global_neck_angle = get_signed_angle_2d(spine_vect, flat_vect)


        self.area_of_roi = 0  # this seems nontrivial... saving for later
        self.estimated_volume = 0

        self.sanity_checks = {
            "neck length": self.neck_length,
            "neck angle error": np.abs(self.relative_neck_angle - np.pi/2),
            "dendrite residual length": self.dendrite_residual,
            "spine circumferance length": self.circumferance,
            "spine area": self.area,
        }

        self.spine_stats = {
            "center_x": self.spine_center_xy[0],
            "center_y": self.spine_center_xy[1],
            "global neck angle": self.global_neck_angle,
            "relative neck angle": self.relative_neck_angle
        }

    @property
    def dist_points_dict(self):
        point_dict = {
            "Spine Center": self.spine_center_xy,
            "Dendrite-Neck Junction": self.spine_neck_xy,
            "Synaptic Cleft": self.synapse_xy,
        }
        return point_dict

    @property
    def circumferance(self):
        curve_dist = distance_along_curve(
            0, len(self.spine_roi), self.spine_roi[0], self.spine_roi[1]
        )
        dist_btw_ends = np.linalg.norm(self.spine_roi[:, 0] - self.spine_roi[:, 1])
        total_circumferance = curve_dist + dist_btw_ends
        return total_circumferance

    @property
    def area(self):
        # NB this will be wrong for any ROIs with acute angles
        start_point = self.spine_roi[:, 0]

        total_area = 0
        for i in range(1, max(self.spine_roi.shape) - 1):
            # compute the area of the triangle formed between the starting point and each edge
            edge_lens = []
            edge_lens.append(np.linalg.norm(start_point - self.spine_roi[:, i]))
            edge_lens.append(np.linalg.norm(start_point - self.spine_roi[:, i + 1]))
            edge_lens.append(
                np.linalg.norm(self.spine_roi[:, i] - self.spine_roi[:, i + 1])
            )
            base = edge_lens.pop(np.argmin(edge_lens))

            height = edge_lens.pop(np.argmin(edge_lens))
            triangle_area = 0.5 * base * height
            total_area += triangle_area
        return total_area


def find_dendritic_segment(spine_dendrite, dendrite_segment_rois):
    min_dist = np.inf
    nearest_segment = None
    spine_dend_center = np.mean(
        roi_to_array(spine_dendrite["x"], spine_dendrite["y"]), axis=1
    )
    for key, segment in dendrite_segment_rois.items():
        dend_dist, min_dist_idx = minimum_distance(spine_dend_center, segment.coords)
        if dend_dist < min_dist:
            min_dist = dend_dist
            nearest_segment = segment
    return nearest_segment, min_dist


def maximum_distance(point, line):
    res_sq = _m_distance(point, line)
    # find index of minimum distance
    max_dist = np.max(res_sq)
    max_dist_idx = np.argmax(res_sq)
    return max_dist, max_dist_idx


def dist_between_2_spines(Spine1, Spine2):
    through_branch = False
    xs = Spine1.parent_dendritic_segment.dend_xs
    ys = Spine1.parent_dendritic_segment.dend_ys
    if Spine1.parent_dendritic_segment == Spine2.parent_dendritic_segment:
        dendritic_distance = distance_along_curve(
            Spine1.spine_neck_idx, Spine2.spine_neck_idx, xs, ys
        )
    else:
        # need to identify path from Spine1 to Spine2...using a for loop and recursion to trace in both directions
        for branch_point in Spine1.parent_dendritic_segment.end_points:
            dendritic_distance = recursive_tree_search(
                Spine1.parent_dendritic_segment, branch_point, Spine2
            )
            if not (dendritic_distance is None):
                through_branch = True
                break
        if dendritic_distance is None:
            return np.nan, through_branch
        end_dist, end_dist_idx = minimum_distance(
            branch_point.end_point_xy, Spine1.parent_dendritic_segment.coords
        )

        spine1_to_first_branch_dist = (
            distance_along_curve(end_dist_idx, Spine1.spine_neck_idx, xs, ys) + end_dist
        )

        dendritic_distance = spine1_to_first_branch_dist + dendritic_distance
    return dendritic_distance, through_branch


def recursive_tree_search(prev_dendrite_segment, prev_branch_point, Spine2):
    # I think that to make this more efficient you would actually want to actually map between all the segments beforehand so you don't have to do it for each spine
    dendritic_distance = None
    # for branch_point in prev_dendrite_segment.end_points:
    #    if not(prev_branch_point == branch_point):
    for child_dendrite in prev_branch_point.child_dendrites:
        if not (child_dendrite == prev_dendrite_segment):
            if child_dendrite == Spine2.parent_dendritic_segment:
                end_dist, end_dist_idx = minimum_distance(
                    prev_branch_point.end_point_xy, child_dendrite.coords
                )
                # ^doing it this way accounts for double counting in case the ends of the manulally annotated segments overlapped a bit, although this is not necessarily handled appropriately in the length of the segment
                dendritic_distance = (
                    distance_along_curve(
                        end_dist_idx,
                        Spine2.spine_neck_idx,
                        child_dendrite.dend_xs,
                        child_dendrite.dend_ys,
                    )
                    + end_dist
                )
                return dendritic_distance
            else:
                # recursion here
                for branch_point in child_dendrite.end_points:
                    if not (prev_branch_point == branch_point):
                        dendritic_distance = recursive_tree_search(
                            child_dendrite, branch_point, Spine2
                        )
                        if not (dendritic_distance is None):
                            dendritic_distance = (
                                child_dendrite.segment_length + dendritic_distance
                            )
                            return dendritic_distance
    return dendritic_distance


def dendritic_distance_matrix(spines):
    return arbitrary_dend_dmat(spines, spines)

def fiducial_dend_mats(spines, fiducials):
    return arbitrary_dend_dmat(spines, fiducials)

def arbitrary_dend_dmat(spines1, spines2):
    my_d_mat = np.zeros((len(spines1), len(spines2)))
    b_mat = np.zeros((len(spines1), len(spines2)))
    for i, Spine1 in enumerate(spines1):
        for j, Spine2 in enumerate(spines2):
            dist_along_dend, through_branch = dist_between_2_spines(Spine1, Spine2)
            my_d_mat[i, j] = dist_along_dend
            b_mat[i, j] = through_branch
    return my_d_mat, b_mat


def euclidian_dmats(spines, visual=False):
    # produce pairwise distance matrix for both euclidean and dendritic distance
    spine_dmats = {}
    plot_types = spines[0].dist_points_dict.keys()
    num_plots = len(plot_types)
    if visual:
        _, plots = plt.subplots(1, num_plots, figsize=(15, 4))
    for i, plot_type in enumerate(plot_types):
        coords = []
        for spine in spines:
            coords.append(spine.dist_points_dict[plot_type])
        try:
            d_mat = distance_matrix(np.array(coords), np.array(coords))
        except Exception as E:
            pass

        new_key = "euclidian distance between " + plot_type + "s"
        spine_dmats[new_key] = d_mat
        if visual:
            plots[i].imshow(d_mat)
            plots[i].title.set_text(new_key)

    return spine_dmats

