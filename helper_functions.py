import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import distance_matrix




def infer_dendrite(dend_rois=None, shaft_roi=None, ts_per_roi = 100):
  #we parameterized this for T, but this will break if the ROIs are not in the right order...
  #will just have to look
  #need a function to parameterise input rois by T as well
  roi_center_xs = []
  roi_center_ys = []
  if shaft_roi:
    for name, roi in shaft_roi.items():
      for x, y in zip(roi['x'], roi['y']):
        roi_center_xs.append(x)
        roi_center_ys.append(y)
  else:
    for dend_roi in dend_rois:
      roi_xs = dend_roi['x']
      roi_ys = dend_roi['y']
      roi_center_xs.append(np.mean(roi_xs))
      roi_center_ys.append(np.mean(roi_ys))

  num_rois = len(roi_center_ys)
  dend_ts = np.arange(0,(num_rois+1)*ts_per_roi,1)
  roi_center_ts = np.arange(ts_per_roi,(num_rois+1)*ts_per_roi,ts_per_roi)
  #dend_ts = np.arange(ts_per_roi,(num_rois-1)*ts_per_roi,1)
  #roi_center_ts = np.arange(0,(num_rois)*ts_per_roi,ts_per_roi)
  #TODOthis shouldn't be an even number of points per ROI, should depend on euclidean distance to next ROI center probably/point on shaft

  #print(roi_center_ts)
  #print(roi_center_xs)
  #print(roi_center_ys)
  #model5 = np.poly1d(np.polyfit(roi_center_xs, roi_center_ys, 5))


  model5_x = np.poly1d(np.polyfit(roi_center_ts, roi_center_xs, 5))
  model5_y = np.poly1d(np.polyfit(roi_center_ts, roi_center_ys, 5))

  dend_xs = model5_x(dend_ts)
  dend_ys = model5_y(dend_ts)
  return dend_ts, dend_xs, dend_ys, roi_center_xs, roi_center_ys



def spine_location_on_dendrite(spine_roi, dendrite_roi, plot=False):
  #compute centers of roi
  #should interpolate all these at some point, but for now use centers
  roi_xs = spine_roi['x']
  roi_ys = spine_roi['y']
  roi_center = np.array([np.mean(roi_xs), np.mean(roi_ys)])

  #interpolate annotated line
  den_xs = dendrite_roi['x']
  den_ys = dendrite_roi['y']
  dense_den_xs = np.arange(min(den_xs), max(den_xs))
  dense_den_ys = np.interp(dense_den_xs, den_xs, den_ys)

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





