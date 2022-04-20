import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import distance_matrix

from config import subfolder_name, collect_images_at_path, precision


#this should actually be in Data_io I think, it isn't actually doing any manipulation
def initialize_dendrite(shaft_rois):

    #first need to figure out where the branch points are

    #dealing with multiple branch points here is tricky
    #also what if there is a 4 way branch instead of 3 way?
    # best way to do it would be to cluster the distances - compute the distance from each segment to each other segment
    #lets assume that each branch point is indicated by 3 segments.
    #actually first just build the code to fail if there are more than 3 - come back to it if you need to
    assert(len(shaft_rois)<4)
    if len(shaft_rois)>1:
        dist_between_end_1 = np.zeroz(len(shaft_rois), len(shaft_rois))
        dist_between_end_2 = np.zeroz(len(shaft_rois), len(shaft_rois))
        #then find the nearest endpoint
        for i, segment_i in enumerate(shaft_rois):
            for j, segment_j in (shaft_rois):
                #have to do this for both endpoints? so for 3 segments there will be 36 distances? that seems so excessive....
                dist_between_end_1[i,j] = np.linalg.norm(diff, axis=0)





        #lastly load them in as segments


class Branch_point():

    def __init__(self):
        branch_point_xy
        child_dendrites = []

class Dendrite_segment():

    def __init__(self, segment_array, branch_point1, branch_point2):
        self.dend_ts
        self.dend_xs
        self.dend_ys
        self.coords =

        #these need to be the branch point objects
        self.end_points = [branch_point1, branch_point2]


        self.segment_length = distance_along_curve(0, len(self.coords), self.coords)

    def other_end_point(self, endpoint_in):
        #wrote this but then included the finctionality in the recursion so I don't think I ever call it
        assert(endpoint_in in self.endpoints)
        for endpoint in self.endpoints:
            if not(endpoint == endpoint_in)
                return(endpoint)
        ValueError('This dendritic seems to have been initialized wrong and'
                   'all endpoints match the one input. '
                   'Circular segments are not allowed.')







def find_dendritic_segment(spine_dendrite, dendrite_segment_rois):
    min_dist = np.inf
    nearest_segment = None
    spine_dend_center = np.mean(roi_to_array(spine_dendrite), axis=1)
    for segment in dendrite_segment_rois:
        dend_dist, min_dist_idx = minimum_distance(spine_dend_center, segment)
        if dend_dist<min_dist:
            min_dist = dend_dist
            nearest_segment = segment
    return neares_segment, min_dist

def roi_to_array(imagej_roi):
    xs = imagej_roi['x']
    ys = imagej_roi['y']
    array = np.vstack((xs, ys))
    return array

class Spine():

    def __init__(self, spine_roi, spine_dendrite_roi, dendrite_segment_rois):

        self.parent_dendritic_segment, self.dendrite_residual = find_dendritic_segment(spine_dendrite_roi, dendrite_segment_rois)

        self.spine_roi = roi_to_array(spine_roi)
        self.spine_center_xy = np.mean(self.spine_roi, axis=1)

        min_dist, self.spine_neck_idx = minimum_distance(self.spine_center_xy, self.parent_dendritic_segment)
    #<<<<<<<<<<<<<<<<<<
        #maybe this should be minimizing the pairwise distance on both end?
        #although this might introduce as many errors as it fixes... need to look at examples
        #grab the appropriate location on the dendrite
        self.spine_neck_xy = self.parent_dendritic_segment[:,self.spine_neck_idx]

        #find the furthest distance from the neck as the synapse (same steps as above)
        max_dist, max_dist_idx = maximum_distance(self.spine_neck_xy, self.spine_roi)
        self.synapse_xy = spine_roi[:,max_dist_idx]

        self.neck_length = np.sqrt((self.spine_neck_xy[0] - self.spine_center_xy[0])**2+(self.spine_neck_xy[1] - self.spine_center_xy[1])**2)

        #subtract the spine neck coordinates to get vectors from the origin for directions to the other 2
        spine_vect = self.spine_center_xy - self.spine_neck_xy
        try:
            den_vect =  self.parent_dendritic_segment[:,self.spine_neck_idx+1] - self.spine_neck_xy
        except Exception as E:
            den_vect =  self.parent_dendritic_segment[:,self.spine_neck_idx-1] - self.spine_neck_xy

        #taken from https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
        c = np.dot(den_vect,spine_vect)/np.linalg.norm(den_vect)/np.linalg.norm(spine_vect) # -> cosine of the angle
        self.neck_angle = np.arccos(np.clip(c, -1, 1))

        #seems like we should start to incorporate the idea that there would be multiple segments here
        #- actually may be over complicating things to include this here....
        #self.endpoints = self.parent_dendritic_segment.endpoints
        #self.dist_to_end_points = []
        #for endpoint in self.endpoints:
        #    self.dist_to_end_points.append()distance_along_curve()

        self.area_of_roi = 0 #this seems nontrivial... saving for later
        self.estimated_volume = 0




    @property
    def dist_points_dict(self):
        point_dict = {
            'Spine Center': self.roi_center,
            'Dendrite-Neck Junction': self.spine_neck,
            'Synaptic Cleft': self.synapse,
        }
        return point_dict


    #def get_dist_to_endpoint(self, endpoint_in):
    #    for i, endpoint in enumerate(self.endpoints):
    #        if endpoint == endpoint_in:
    #            return self.dist_to_end_points[i]
    #    ValueError('The point passed is not one of the endpoints of this segments')


#How do we handle multiple branch points
def dist_between_2_spines(Spine1, Spine2):
    if Spine1.parent_dendritic_segment == Spine2.parent_dendritic_segment:
        dendritic_distance = distance_along_curve(Spine1.spine_neck_t, Spine2.spine_neck_t, Spine1.parent_dendritic_segment.coords)
    else:
        #need to identify path from Spine1 to Spine2...
        #whats the best way to do this? start with the hope that they are on the same?
        #somehow we need to build an understanding of the whole tree structure....
        #just do it recursively with a for loop. trace both directions as far as they go until you get to a match
        for branch_point in Spine1.parent_dendritic_segment.end_points:
            dendritic_distance = recursive_tree_search(Spine1.parent_dendritic_segment, branch_point, Spine2)
            if not(dendritic_distance is None):
                break
        else:
            ValueError('Spine2 was not found on a connected segment')
        end_dist, end_dist_idx = minimum_distance(branch_point.branch_point_xy, Spine1.parent_dendritic_segment.coords)
        spine1_to_first_branch_dist = distance_along_curve(end_dist_idx, Spine1.spine_neck_t, Spine1.parent_dendritic_segment.coords) + end_dist_idx
        dendritic_distance = spine1_to_first_branch_dist + dendritic_distance
    return dendritic_distance

def recursive_tree_search(prev_dendrite_segment, prev_branch_point, Spine2):
    #I think that to make this more efficient you would actually want to actually map between all the segments beforehand so you don't have to do it for each spine
    dendritic_distance = None
    for branch_point in prev_dendrite_segment.end_points:
        if not(prev_branch_point == branch_point):
            for child_dendrite in branch_point.child_dendrites:
                if not(child_dendrite == dendrite_segment):
                    if child_dendrite == Spine2.parent_dendritic_segment:
                        end_dist, end_dist_idx = minimum_distance(branch_point.branch_point_xy, child_dendrite.coords)
                        #^doing it this way accounts for double counting in case the ends of the manulally annotated segments overlapped a bit, although this is not necessarily handled appropriately in the length of the segment
                        dendritic_distance = distance_along_curve(end_dist_idx, Spine2.spine_neck_t, child_dendrite.coords)+end_dist
                        return dendritic_distance
                    else:
                        #recursion here
                        dendritic_distance = recursive_tree_search(child_dendrite, branch_point, Spine2)
                        if not(dendritic_distance is None):
                            dendritic_distance = child_dendrite.segment_length +dendritic_distance
                            return dendritic_distance
    return dendritic_distance


def dendritic_distance_matrix(spines):
    #print(spine_neck_t)

    my_d_mat = np.zeros((len(spines), len(spines)))
    for i, Spine1 in enumerate(spines):
        for j, Spine2 in enumerate(spines):
            dist_along_dend = dist_between_2_spines(Spine1, Spine2)
            my_d_mat[i,j] = dist_along_dend
    #for i, t_1 in enumerate(spine_neck_t):
    #    for j, t_2 in enumerate(spine_neck_t):
    #        dist_along_dend = distance_along_curve(t_1, t_2, dendrite_roi)
    #        #print(dist_along_dend*.09)
    #        my_d_mat[i,j] = dist_along_dend
    return my_d_mat




def distance_along_curve(t_1, t_2, dendrite_roi):
    dist = 0
    start_t = min(t_1,t_2)
    end_t = max(t_1,t_2)
    #print(start_t, end_t)
    #this is where we could introduce direction
    for i in range(start_t,end_t):
        x_1 = dendrite_roi['x'][i]
        x_2 = dendrite_roi['x'][i+1]
        y_1 = dendrite_roi['y'][i]
        y_2 = dendrite_roi['y'][i+1]
        #print(x_1, x_2, y_1, y_2)
        dist += np.sqrt((x_1 - x_2)**2+(y_1-y_2)**2)
    return dist


def _m_distance(point, line):
    #den_coords = np.vstack((np.array(den_xs), np.array(den_ys)))
    num_elements = np.shape(line)[1]
    points = np.tile(point, (num_elements,1)).T
    diff = line - points#, (1, num_elements))#.reshape(2,num_elements)
    res_sq = np.linalg.norm(diff, axis=0)
    return res_sq

def minimum_distance(point, line):
    res_sq = _m_distance(point, line)
    #find index of minimum distance
    min_dist = np.min(res_sq)
    min_dist_idx = np.argmin(res_sq)
    return min_dist, min_dist_idx

def maximum_distance(point, line):
    res_sq = _m_distance(point, line)
    #find index of minimum distance
    max_dist = np.max(res_sq)
    max_dist_idx = np.argmax(res_sq)
    return max_dist, max_dist_idx






def spine_location_on_dendrite(spine_roi, dendrite_roi, plot=False):
    #compute centers of roi
    #should interpolate all these at some point, but for now use centers
    #roi_xs = spine_roi['x']
    #roi_ys = spine_roi['y']
    #roi_center = np.array([np.mean(roi_xs), np.mean(roi_ys)])
    spine_roi = roi_to_array(spine_roi)
    roi_center = np.mean(spine_roi, axis=1)

    #interpolate annotated line
    #dense_den_xs = dendrite_roi['x']
    #dense_den_ys = dendrite_roi['y']
    #dense_den_xs = np.arange(min(den_xs), max(den_xs))
    #dense_den_ys = np.interp(dense_den_xs, den_xs, den_ys)

    #compute distances from all points on annotated line to center
    #den_coords = np.vstack((dense_den_xs, dense_den_ys)) #<- use this one when not in debug mode
    den_coords = roi_to_array(dendrite_roi)

    #den_coords = np.vstack((np.array(den_xs), np.array(den_ys)))
    min_dist, min_dist_idx = minimum_distanct(roi_center, den_coords)
    #num_elements = np.shape(den_coords)[1]
    #roi_centers = np.tile(roi_center, (num_elements,1)).T
    #diff = den_coords - roi_centers#, (1, num_elements))#.reshape(2,num_elements)
    #res_sq = np.linalg.norm(diff, axis=0)

    #find index of minimum distance
    #min_dist = np.argmin(res_sq)#this is the index where the distance was minimized, this also corresponds to the t that we used to parameterize the curve
    spine_neck_t = min_dist_idx

#<<<<<<<<<<<<<<<<<<
    #todo this should really be minimizing the pairwise distance on both ends
    #although this might introduce as many errors as it fixes... need to look at examples

    #grab the appropriate location on the dendrite
    spine_neck = den_coords[:,min_dist]

    #find the furthest distance from the neck as the synapse (same steps as above)
    #roi_coords = np.vstack((roi_xs, roi_ys))

    max_dist, max_dist_idx = maximum_distance(spine_neck, spine_roi)

    #num_elements = np.shape(roi_coords)[1]
    #spine_necks = np.tile(spine_neck, (num_elements,1)).T
    #diff = roi_coords - spine_necks
    #res_sq = np.linalg.norm(diff, axis=0)
    #min_dist = np.argmax(res_sq)
    synapse = spine_roi[:,max_dist_idx]

    #plot this line
    if plot is True:
        #plot the ROI in question
        plt.plot(spine_roi['x'],spine_roi['y'])
        #plot the center of the roi as an X
        plt.scatter(roi_center[0], roi_center[1], marker='x')
        #plot the dendrite
        plt.plot(dense_den_xs, dense_den_ys)
        #plot line from center to dendrite
        plt.plot([spine_neck[0], roi_center[0]], [spine_neck[1], roi_center[1]])
        #plot the location of the spine neck and the synaps
        plt.scatter(spine_neck[0], spine_neck[1], marker='o')
        plt.scatter(synapse[0], synapse[1], marker='+')

    return roi_center, spine_neck, synapse, spine_neck_t


def spine_neck_for_all_rois(spine_rois, dend_rois, dendrite_roi, plot=True):
    spines = []

    #spine_coords = {}
    #spine_coords['centers'] = []
    #spine_coords['spine_necks'] = []
    #spine_neck_t_list = []
    #spine_coords['synapses'] = []

    for spine_roi in spine_rois:
        spines.append()
    #    roi_center, spine_neck, synapse, spine_neck_t = spine_location_on_dendrite(spine_roi, dendrite_roi, plot=plot)
    #    spine_coords['centers'].append(roi_center)
    #    spine_coords['spine_necks'].append(spine_neck)
    #    spine_neck_t_list.append(spine_neck_t)
    #    spine_coords['synapses'].append(synapse)
    return spines#spine_coords, spine_neck_t_list


def euclidian_dmats(spines):
    #produce pairwise distance matrix for both euclidean and dendritic distance
    spine_dmats = {}
    plot_types = spines[0].dist_points_dict.keys()
    num_plots = len(plot_types)
    _, plots = plt.subplots(1, num_plots, figsize=(15, 4))
    for i, plot_type in enumerate(plot_types):
        coords = np.zeros(len(spines)):
        for i, spine in spines:
            coords[i] = spine.dist_points_dict[plot_type]
        try:
            d_mat = distance_matrix(np.array(coords), np.array(coords))
        except Exception as E:
            pass

        new_key = 'euclidian distance between '+plot_type+'s'
        spine_dmats[new_key] = d_mat
        plots[i].imshow(d_mat)
        plots[i].title.set_text(new_key)
    #for i, (key, coords) in enumerate(spine_coords.items()):
    #    try:
    #        d_mat = distance_matrix(np.array(coords), np.array(coords))
    #    except Exception as E:
    #        pass
    ##<<<<<<<<<<<<<<<<<<
    #    new_key = 'euclidian_'+key
    #    spine_dmats[new_key] = d_mat
    #    plots[i].imshow(d_mat)
    #    plots[i].title.set_text(new_key)
    return spine_dmats






def neck_length(spine_coords):
        neck_lengths = []
        for i, _ in enumerate(spine_coords['centers']):
                roi_center = spine_coords['centers'][i]
                spine_neck = spine_coords['spine_necks'][i]

                neck_length = np.sqrt((spine_neck[0] - roi_center[0])**2+(spine_neck[1] - roi_center[1])**2)
                neck_lengths.append(neck_length)
        return neck_lengths

def dendrite_residual(dend_rois, dendrite_roi):
        dendrite_residuals = []
        dense_den_xs = dendrite_roi['x']
        dense_den_ys = dendrite_roi['y']
        den_coords = np.vstack((dense_den_xs, dense_den_ys)) #<- use this one when not in debug mode
        num_elements = np.shape(den_coords)[1]

        roi_center_xs, roi_center_ys = get_dend_path_from_dend_rois(dend_rois)

        for roi_center_x, roi_center_y in zip(roi_center_xs, roi_center_ys):
                roi_center = np.array([np.mean(roi_center_x), np.mean(roi_center_y)])
                #den_coords = np.vstack((np.array(den_xs), np.array(den_ys)))

                roi_centers = np.tile(roi_center, (num_elements,1)).T
                diff = den_coords - roi_centers#, (1, num_elements))#.reshape(2,num_elements)
                res_sq = np.linalg.norm(diff, axis=0)

                #find index of minimum distance
                min_dist = np.argmin(res_sq)#this is the index where the distance was minimized, this also corresponds to the t that we used to parameterize the curve
                spine_neck_t = min_dist

                #<<<<<<<<<<<<<<<<<<
                #grab the appropriate location on the dendrite
                closest_point = den_coords[:,min_dist]

                residual_length = np.sqrt((closest_point[0] - roi_center[0])**2+(closest_point[1] - roi_center[1])**2)
                dendrite_residuals.append(residual_length)
        return dendrite_residuals

def neck_angle(spine_coords, spine_neck_t, dendrite_roi):
        neck_angles = []
        for i, _ in enumerate(spine_coords['centers']):
                roi_center = spine_coords['centers'][i]
                spine_neck = spine_coords['spine_necks'][i]

                #grab the next point alonge the dendrite
                t = spine_neck_t[i] +1
                den_x = dendrite_roi['x'][t]
                den_y = dendrite_roi['y'][t]

                #subtract the spine neck coordinates to get vectors from the origin for directions to the other 2
                spine_vect = roi_center - spine_neck
                den_vect = np.array([den_x,den_y]) - spine_neck

                #taken from https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
                c = np.dot(den_vect,spine_vect)/np.linalg.norm(den_vect)/np.linalg.norm(spine_vect) # -> cosine of the angle
                neck_angle = np.arccos(np.clip(c, -1, 1))

                neck_angles.append(neck_angle)
        return neck_angles







def get_dend_path_from_shaft_roi(shaft_roi):
    roi_center_xs = []
    roi_center_ys = []
    for name, roi in shaft_roi.items():
        for x, y in zip(roi['x'], roi['y']):
            roi_center_xs.append(x)
            roi_center_ys.append(y)
    return roi_center_xs, roi_center_ys

def get_dend_path_from_dend_rois(dend_rois):
        roi_center_xs = []
        roi_center_ys = []
        for dend_roi in dend_rois:
                roi_xs = dend_roi['x']
                roi_ys = dend_roi['y']
                roi_center_xs.append(np.mean(roi_xs))
                roi_center_ys.append(np.mean(roi_ys))
        return roi_center_xs, roi_center_ys

def infer_dendrite(dend_rois=None, shaft_roi=None, ts_per_roi = 100):
        if shaft_roi:
                roi_center_xs, roi_center_ys = get_dend_path_from_shaft_roi(shaft_roi)
        else:
                roi_center_xs, roi_center_ys = get_dend_path_from_dend_rois(dend_rois)


                #num_rois = len(roi_center_ys)
                #dend_ts = np.arange(0,(num_rois+1)*ts_per_roi,1)
                #roi_center_ts = np.arange(ts_per_roi,(num_rois+1)*ts_per_roi,ts_per_roi)
                #dend_ts = np.arange(ts_per_roi,(num_rois-1)*ts_per_roi,1)
                #roi_center_ts = np.arange(0,(num_rois)*ts_per_roi,ts_per_roi)
                #TODOthis shouldn't be an even number of points per ROI, should depend on euclidean distance to next ROI center probably/point on shaft

                #print(roi_center_ts)
                #print(roi_center_xs)
                #print(roi_center_ys)
                #model5 = np.poly1d(np.polyfit(roi_center_xs, roi_center_ys, 5))


                #model5_x = np.poly1d(np.polyfit(roi_center_ts, roi_center_xs, 5))
                #model5_y = np.poly1d(np.polyfit(roi_center_ts, roi_center_ys, 5))

                #It wasn't working to just parameterize by T
                #need to first check which of X or Y has higher varience
                #then use that one as the x for our best fit line
        num_steps = 1000
        if np.var(roi_center_xs) >= np.var(roi_center_ys):
                #start = min(roi_center_xs)
                #end = max(roi_center_xs)
                #extend = int(end-start*.1)
                #start = start-extend
                #end = end+extend
                start = 0
                end = 512
                step = (end - start)/num_steps
                dend_xs = np.arange(start, end, step)
                model5 = np.poly1d(np.polyfit(roi_center_xs, roi_center_ys, 5))
                dend_ys = model5(dend_xs)
                #we want to constrain both the x and y values to the pixels limitations of the image
                mask = np.logical_and(dend_ys>start,dend_ys<end)
                dend_ys = dend_ys[mask]
                dend_xs = dend_xs[mask]
        else:
        #start = min(roi_center_ys)
        #end = max(roi_center_ys)
                start = 0
                end = 512
                step = (end - start)/num_steps
                dend_ys = np.arange(start, end, step)
                model5 = np.poly1d(np.polyfit(roi_center_ys, roi_center_xs, 5))
                dend_xs = model5(dend_ys)    #then paramterize that line by T (over the whole width, min to max)
                #we want to constrain both the x and y values to the pixels limitations of the image
                mask = np.logical_and(dend_xs>start,dend_xs<end)
                dend_ys = dend_ys[mask]
                dend_xs = dend_xs[mask]

        dend_ts = np.arange(0, len(dend_xs), 1)

        #optionally could go back and order the ROI centers, and then re-fit X and Y seperately
        #this would allow us to deal with curves better, but I think picking the right input dimensiont is the better solution
        #since it still won't owrk for branch points.

        #dend_xs = model5_x(dend_ts)
        #dend_ys = model5_y(dend_ts)
        return dend_ts, dend_xs, dend_ys, roi_center_xs, roi_center_ys


