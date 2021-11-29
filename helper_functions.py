import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import distance_matrix


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
        dend_xs = model5(dend_ys)  #then paramterize that line by T (over the whole width, min to max)
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



def spine_location_on_dendrite(spine_roi, dendrite_roi, plot=False):
  #compute centers of roi
  #should interpolate all these at some point, but for now use centers
  roi_xs = spine_roi['x']
  roi_ys = spine_roi['y']
  roi_center = np.array([np.mean(roi_xs), np.mean(roi_ys)])

  #interpolate annotated line
  dense_den_xs = dendrite_roi['x']
  dense_den_ys = dendrite_roi['y']
  #dense_den_xs = np.arange(min(den_xs), max(den_xs))
  #dense_den_ys = np.interp(dense_den_xs, den_xs, den_ys)

  #compute distances from all points on annotated line to center
  den_coords = np.vstack((dense_den_xs, dense_den_ys)) #<- use this one when not in debug mode
  #den_coords = np.vstack((np.array(den_xs), np.array(den_ys)))
  num_elements = np.shape(den_coords)[1]
  roi_centers = np.tile(roi_center, (num_elements,1)).T
  diff = den_coords - roi_centers#, (1, num_elements))#.reshape(2,num_elements)
  res_sq = np.linalg.norm(diff, axis=0)

  #find index of minimum distance
  min_dist = np.argmin(res_sq)#this is the index where the distance was minimized, this also corresponds to the t that we used to parameterize the curve
  spine_stem_t = min_dist

#<<<<<<<<<<<<<<<<<<
  #todo this should really be minimizing the pairwise distance on both ends
  #although this might introduce as many errors as it fixes... need to look at examples

  #grab the appropriate location on the dendrite
  spine_stem = den_coords[:,min_dist]

  #find the furthest distance from the stem as the synapse (same steps as above)
  roi_coords = np.vstack((roi_xs, roi_ys))
  num_elements = np.shape(roi_coords)[1]
  spine_stems = np.tile(spine_stem, (num_elements,1)).T
  diff = roi_coords - spine_stems
  res_sq = np.linalg.norm(diff, axis=0)
  min_dist = np.argmax(res_sq)
  synapse = roi_coords[:,min_dist]

  #plot this line
  if plot is True:
    #plot the ROI in question
    plt.plot(spine_roi['x'],spine_roi['y'])
    #plot the center of the roi as an X
    plt.scatter(roi_center[0], roi_center[1], marker='x')
    #plot the dendrite
    plt.plot(dense_den_xs, dense_den_ys)
    #plot line from center to dendrite
    plt.plot([spine_stem[0], roi_center[0]], [spine_stem[1], roi_center[1]])
    #plot the location of the spine stem and the synaps
    plt.scatter(spine_stem[0], spine_stem[1], marker='o')
    plt.scatter(synapse[0], synapse[1], marker='+')

  return roi_center, spine_stem, synapse, spine_stem_t


def spine_stem_for_all_rois(spine_rois, dendrite_roi, plot=True):
  spine_coords = {}
  spine_coords['centers'] = []
  spine_coords['spine_stems'] = []
  spine_stem_t_list = []
  spine_coords['synapses'] = []

  for spine_roi in spine_rois:
    roi_center, spine_stem, synapse, spine_stem_t = spine_location_on_dendrite(spine_roi, dendrite_roi, plot=plot)
    spine_coords['centers'].append(roi_center)
    spine_coords['spine_stems'].append(spine_stem)
    spine_stem_t_list.append(spine_stem_t)
    spine_coords['synapses'].append(synapse)
  return spine_coords, spine_stem_t_list


def euclidian_dmats(spine_coords):
    #produce pairwise distance matrix for both euclidean and dendritic distance
    spine_dmats = {}
    _, plots = plt.subplots(1, 4, figsize=(15, 4))
    for i, (key, coords) in enumerate(spine_coords.items()):
      try:
        d_mat = distance_matrix(np.array(coords), np.array(coords))
      except Exception as E:
        pass
    #<<<<<<<<<<<<<<<<<<
      #todo distance along dendrite should be calculated differently, to follow the dendrite path
      new_key = 'euclidian_'+key
      spine_dmats[new_key] = d_mat
      plots[i].imshow(d_mat)
      plots[i].title.set_text(new_key)
    return spine_dmats


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

def dendritic_distance(spine_stem_t, dendrite_roi):
    #print(spine_stem_t)
    my_d_mat = np.zeros((len(spine_stem_t), len(spine_stem_t)))
    for i, t_1 in enumerate(spine_stem_t):
      for j, t_2 in enumerate(spine_stem_t):
        dist_along_dend = distance_along_curve(t_1, t_2, dendrite_roi)
        #print(dist_along_dend*.09)
        my_d_mat[i,j] = dist_along_dend
    return my_d_mat

def stem_length(spine_coords):
    stem_lengths = []
    for i, _ in enumerate(spine_coords['centers']):
        roi_center = spine_coords['centers'][i]
        spine_stem = spine_coords['spine_stems'][i]

        stem_length = np.sqrt((spine_stem[0] - roi_center[0])**2+(spine_stem[1] - roi_center[1])**2)
        stem_lengths.append(stem_length)
    return stem_lengths

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
        spine_stem_t = min_dist

        #<<<<<<<<<<<<<<<<<<
        #grab the appropriate location on the dendrite
        closest_point = den_coords[:,min_dist]

        residual_length = np.sqrt((closest_point[0] - roi_center[0])**2+(closest_point[1] - roi_center[1])**2)
        dendrite_residuals.append(residual_length)
    return dendrite_residuals

def stem_angle(spine_coords, spine_stem_t, dendrite_roi):
    stem_angles = []
    for i, _ in enumerate(spine_coords['centers']):
        roi_center = spine_coords['centers'][i]
        spine_stem = spine_coords['spine_stems'][i]

        #grab the next point alonge the dendrite
        t = spine_stem_t[i] +1
        den_x = dendrite_roi['x'][t]
        den_y = dendrite_roi['y'][t]

        #subtract the spine stem coordinates to get vectors from the origin for directions to the other 2
        spine_vect = roi_center - spine_stem
        den_vect = np.array([den_x,den_y]) - spine_stem

        #taken from https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
        c = np.dot(den_vect,spine_vect)/np.linalg.norm(den_vect)/np.linalg.norm(spine_vect) # -> cosine of the angle
        stem_angle = np.arccos(np.clip(c, -1, 1))

        stem_angles.append(stem_angle)
    return stem_angles



