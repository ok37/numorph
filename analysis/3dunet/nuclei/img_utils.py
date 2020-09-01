import numpy as np
import cc3d
import cv2
from scipy.ndimage.measurements import center_of_mass
from scipy import ndimage
# import matplotlib
# matplotlib.use('QT5Agg')
# import matplotlib.pyplot as plt
# from matplotlib.colors import ListedColormap
from skimage.transform import resize
from scipy.spatial import cKDTree


# %% Create weights for patch indexes
def create_3d_weight_patch(chunk_size=[160, 160, 48], overlap=[2, 2, 1]):
    chunk_weight_patch = np.ones(chunk_size)

    # Weights calculated by linear interpolation from border to indicated overlap for each dimension
    # z_top = np.linspace(0, 1, overlap[0] + 2)[1:-1]
    # z_bottom = 1 - z_top
    x_left = np.zeros(overlap[0])
    x_left[int(overlap[0] / 2):] = 1
    x_right = abs(1 - x_left)
    x_string = np.hstack((x_left, np.ones([chunk_size[0] - overlap[0] * 2]), x_right))
    chunk_weight_patch = np.einsum('i,ijk->ijk', x_string, chunk_weight_patch)

    # y_top = np.linspace(0, 1, overlap[1] + 2)[1:-1]
    # y_bottom = 1 - y_top
    y_top = np.zeros(overlap[1])
    y_top[int(overlap[1] / 2):] = 1
    y_bottom = abs(1 - y_top)
    y_string = np.hstack((y_top, np.ones([chunk_size[1] - overlap[1] * 2]), y_bottom))
    chunk_weight_patch = np.einsum('j,ijk->ijk', y_string, chunk_weight_patch)

    # x_left = np.linspace(0, 1, overlap[2] + 2)[1:-1]
    # x_right = 1 - x_left
    z_top = np.zeros(overlap[2])
    z_top[int(overlap[2] / 2):] = 1
    z_bottom = abs(1 - z_top)
    z_string = np.hstack((z_top, np.ones([chunk_size[2] - overlap[2] * 2]), z_bottom))
    chunk_weight_patch = np.einsum('k,ijk->ijk', z_string, chunk_weight_patch)

    return chunk_weight_patch


# %% Pad image chunks so that they match in size
def pad_to_match_sizes(chunk_top, chunk_bottom, chunk_top_positions, chunk_bottom_positions):
    shifts = chunk_top_positions - chunk_bottom_positions

    # Top side
    # Pad bottom chunk
    if shifts[0] < 0:
        chunk_bottom = np.pad(chunk_bottom, ((0, 0), (abs(shifts[0]), 0), (0, 0)), 'constant')
    # Pad top chunk
    elif shifts[0] > 0:
        chunk_top = np.pad(chunk_top, ((0, 0), (abs(shifts[0]), 0), (0, 0)), 'constant')

    # Left side
    # Pad bottom chunk
    if shifts[0] < 0:
        chunk_bottom = np.pad(chunk_bottom, ((0, 0), (0, 0), (abs(shifts[1]), 0)), 'constant')
    # Pad top chunk
    elif shifts[0] > 0:
        chunk_top = np.pad(chunk_top, ((0, 0), (0, 0), (abs(shifts[1]), 0)), 'constant')

    pad_ends = np.array(chunk_top.shape[1:]) - np.array(chunk_bottom.shape[1:])

    # Bottom side
    # Pad bottom chunk
    if pad_ends[0] > 0:
        chunk_bottom = np.pad(chunk_bottom, ((0, 0), (0, abs(pad_ends[0])), (0, 0)), 'constant')
    # Pad top chunk
    elif pad_ends[0] < 0:
        chunk_top = np.pad(chunk_top, ((0, 0), (0, abs(pad_ends[0])), (0, 0)), 'constant')

    # Right side
    # Pad bottom chunk
    if pad_ends[0] > 0:
        chunk_bottom = np.pad(chunk_bottom, ((0, 0), (0, 0), (0, abs(pad_ends[1]))), 'constant')
    # Pad top chunk
    elif pad_ends[0] < 0:
        chunk_top = np.pad(chunk_top, ((0, 0), (0, 0), (0, abs(pad_ends[1]))), 'constant')

    return shifts, chunk_top, chunk_bottom


# %% Display overlayed image
def display_overlayed_images(img1, img2):
    color_full = np.linspace(0, 1, 256)
    color_empty = np.linspace(0, 0, 256)
    new_cmap1 = np.array([color_full, color_empty, color_empty])
    new_cmap2 = np.array([color_empty, color_full, color_empty])

    red = ListedColormap(new_cmap1.T)
    green = ListedColormap(new_cmap2.T)

    height = np.min([img1[0].shape, img2[0].shape])
    width = np.min([img1[1].shape, img2[1].shape])

    plt.imshow(img1[0:height, 0:width], cmap=red, vmin=0, vmax=0.1)
    plt.imshow(img2[0:height, 0:width], cmap=green, vmin=0, vmax=0.1, alpha=0.3)

    return


# %% Display overlayed centroids
def display_overlayed_centroids(img_list, slice, centroids, mask):
    cell_img = cv2.imread(img_list[0][slice], -1)
    centroid_img = np.zeros(cell_img.shape)

    z_range = slice
    centroids = centroids[np.isin(centroids[:, 2], z_range)]

    mask_resized = resize(mask[:, :, slice], cell_img.shape, order=0)
    cell_img[mask_resized < 1] = 0

    print('Displaying', len(centroids), 'detected cells')

    for i in range(len(centroids)):
        [y, x, z] = centroids[i].round().astype(int)
        centroid_img[y, x - 1:x + 2] = 1
        centroid_img[y - 1:y + 2, x] = 1

    color_full = np.linspace(0, 1, 256)
    color_empty = np.linspace(0, 0, 256)
    new_cmap1 = np.array([color_full, color_empty, color_empty])
    new_cmap2 = np.array([color_empty, color_full, color_empty])

    red = ListedColormap(new_cmap1.T)
    green = ListedColormap(new_cmap2.T)

    p1 = plt.imshow(cell_img, cmap='gray', vmin=0, vmax=5000)
    p2 = plt.imshow(centroid_img, cmap=red, vmin=0, vmax=0.1, alpha=0.5)

    return p2


# %% Calculate upper and lower bounds for rescaling intensities
def calculate_rescaling_intensities(img_list, sampling_range=10, low_pct=0, high_pct=1):
    """ Input is image list"""

    start = round(len(img_list) / 2 - sampling_range / 2)
    end = round(len(img_list) / 2 + sampling_range / 2)

    images = [cv2.imread(file, -1) for file in img_list[start:end]]
    images = np.asarray(images, dtype=np.float32)

    n_images, height, width = images.shape

    images_reduced = images[:, round(height * 0.4):round(height * 0.6), round(width * 0.4):round(width * 0.6)]

    low_val = images_reduced.min()
    high_val = images_reduced.max()

    return low_val, high_val


# %% Read measure channel intensities at centroid positions
def measure_colocalization(centroids, img_list):
    images = [cv2.imread(file, -1).astype('float32') for file in img_list]

    kernel = np.ones((3, 3), np.float32) / 9
    images = [cv2.filter2D(img, -1, kernel) for img in images]

    intensities = []
    for i, img in enumerate(images):
        cen_z = centroids[centroids[:, 2] == i, ]
        if len(cen_z) > 0:
            intensities.append([images[i][tuple(cen_z[c, 0:2].astype(int))] for c in range(len(cen_z))])

    intensities = np.concatenate(intensities)

    return intensities


# %% Updated patch based prediction
def patch_wise_prediction2(model, data, overlap=(16, 16, 8), pred_threshold=0.5, int_threshold=None):
    patch_shape = (1,) + model.input_shape[1:]
    data_shape = data.shape

    # Get tiling positions for each axis
    [idx, tiles] = get_axis_locations(data_shape[2:], patch_shape[2:], overlap)

    # Calculate amount of padding and pad
    expanded_patch_size = [idx[0, 1][-1], idx[1, 1][-1], idx[2, 1][-1]]
    padding = expanded_patch_size - np.array(data_shape[2:])
    padding = padding / 2
    padding = tuple(padding.astype(int))
    padded_data = np.pad(data[0, 0], ((padding[0], padding[0]), (padding[1], padding[1]), (padding[2], padding[2])),
                         'median')

    # Predict patches
    prediction = []
    for i in range(tiles[0]):
        for j in range(tiles[1]):
            for k in range(tiles[2]):
                img_chunk = padded_data[idx[0, 0][i]:idx[0, 1][i], idx[1, 0][j]:idx[1, 1][j], idx[2, 0][k]:idx[2, 1][k]]
                img_chunk = np.reshape(img_chunk, (1, 1) + img_chunk.shape)
                prd_chunk = model.predict(img_chunk)
                prediction.append(prediction_to_centroids(prd_chunk, img_chunk, pred_threshold, int_threshold))

    # Trim overlapping regions saving only the 50% closest to the edge
    ov = np.array(tuple(i / 2 for i in overlap), dtype=int)
    idx[0, 0][1:] += ov[0]
    idx[0, 1][:-1] -= ov[0]
    idx[1, 0][1:] += ov[1]
    idx[1, 1][:-1] -= ov[1]
    idx[2, 0][1:] += ov[2]
    idx[2, 1][:-1] -= ov[2]

    # Find corresponding start,end positions in the tiles
    a1 = np.repeat(ov[0], tiles[0]).astype(int)
    a1[0] = 0
    a2 = np.repeat(patch_shape[2] - ov[1], tiles[0]).astype(int)
    a2[-1] = patch_shape[2]
    b1 = np.repeat(ov[1], tiles[1]).astype(int)
    b1[0] = 0
    b2 = np.repeat(patch_shape[3] - ov[1], tiles[1]).astype(int)
    b2[-1] = patch_shape[3]
    c1 = np.repeat(ov[2], tiles[2]).astype(int)
    c1[0] = 0
    c2 = np.repeat(patch_shape[4] - ov[2], tiles[2]).astype(int)
    c2[-1] = patch_shape[4]

    # Save assembled image into an array of the same size as the input data
    prediction_final = np.zeros(padded_data.shape)
    idx2 = 0
    for i in range(tiles[0]):
        for j in range(tiles[1]):
            for k in range(tiles[2]):
                final_chunk = np.squeeze(prediction[idx2])[a1[i]:a2[i], b1[j]:b2[j], c1[k]:c2[k]]
                prediction_final[idx[0, 0][i]:idx[0, 1][i], idx[1, 0][j]:idx[1, 1][j], idx[2, 0][k]:idx[2, 1][k]] = \
                    final_chunk
                idx2 += 1

    prediction_final = prediction_final[padding[0]:-padding[0], padding[1]:-padding[1], padding[2]:-padding[2]]

    prediction_final = remove_touching_centroids(prediction_final)

    return prediction_final


# %% Get locations of tiles based on the data size, patch size, and amount of overlap between tiles
def get_axis_locations(data_shape, patch_shape, overlap):
    idx = []
    tiles = []
    for i in range(3):
        x_tiles = 0
        x = [-overlap[i] / 2]
        x_length = 0
        while x_length < (data_shape[i] + overlap[i] / 2):
            x_length -= overlap[i] / 2
            x_tiles += 1
            x_length += patch_shape[i] - overlap[i] / 2
            x.append(x_length)
        x1 = np.array(x[0:-1], dtype=int)
        x1[0] = 0
        x2 = x1 + patch_shape[i]
        idx.append([x1, x2])
        tiles.append(x_tiles)

    return np.asarray(idx), tiles


# %% Convert chunks to centroids
def prediction_to_centroids(prd_chunk, img_chunk, pred_threshold=0.5, int_threshold=None):
    # Calculate connected components and determine centroid positions
    # Calculate final prediction mask and label connected components
    prd_chunk = np.squeeze(prd_chunk)
    img_chunk = np.squeeze(img_chunk)
    output_thresh = np.where(prd_chunk > pred_threshold, 1, 0)

    labels_out = cc3d.connected_components(output_thresh)


    n_cells = np.max(labels_out)

    # Find centroids
    centroids = center_of_mass(output_thresh, labels_out, index=np.arange(1, n_cells + 1))
    centroids = np.asarray(centroids, dtype=int).round()

    # Remove cells with low intensity
    if int_threshold is not None:
        int_high = [img_chunk[tuple(centroids[c].astype(int))] > int_threshold for c in
                    range(len(centroids))]
        centroids = centroids[int_high]

    centroid_chunk = np.zeros(prd_chunk.shape)
    centroid_chunk[centroids[:, 0], centroids[:, 1], centroids[:, 2]] = 1

    return centroid_chunk


# %% Prune centroids close to each other
def remove_touching_centroids(cen_img, radius=2):
    # Input is a binary centroid image. Get centroid locations by simple thresholding.
    centroids = np.where(cen_img > 0)

    # Create cKDTree object
    tree = cKDTree(np.array(centroids).T)

    # Check for indexes within the indicated range
    near = tree.query_ball_point(np.array(centroids).T, radius)

    # Get which centroids are touching
    jj = []
    for j, n in enumerate(near):
        if len(near[j]) > 1:
            jj.append(n)
    jj = [j[0] for j in jj]
    rm_idx = np.unique(jj)

    # Keep only non-touching centroids
    centroids_new = np.delete(centroids, rm_idx, 1)

    cen_new = np.zeros(cen_img.shape, dtype=cen_img.dtype)
    cen_new[centroids_new[0], centroids_new[1], centroids_new[2]] = 1

    return cen_new


#%% Remove centroids from data frame
def remove_touching_df(cen_df, radius=2):

    # Get centroids positions from data frame
    centroids = np.array(cen_df.iloc[:, :3])

    # Create cKDTree object
    tree = cKDTree(centroids)

    # Check for indexes within the indicated range
    near = tree.query_ball_point(centroids, radius)

    # Get which centroids are touching
    jj = []
    for j, n in enumerate(near):
        if len(near[j]) > 1:
            jj.append(n)
    jj = [j[0] for j in jj]
    rm_idx = np.unique(jj)

    # Keep only non-touching centroids
    centroids = np.delete(centroids, rm_idx, 1)

    return centroids
