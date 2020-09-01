import os
import cv2
import numpy as np
import cc3d
import argparse

from datetime import datetime
from scipy.io import loadmat
from scipy import ndimage
from skimage.transform import resize
from skimage.exposure import rescale_intensity
from scipy.spatial import cKDTree

from unet3d.training import load_old_model
from nuclei.img_utils import create_3d_weight_patch, calculate_rescaling_intensities, measure_colocalization

parser = argparse.ArgumentParser(description='Predict from validation dataset.')
parser.add_argument('--matlab', metavar='matlab', type=str, nargs='+',
                    help='Whether called by matlab function ')
parser.add_argument('--r', metavar='r', type=str, nargs='+',
                    help='Resolution tag')
parser.add_argument('--g', metavar='g', type=str, nargs='+',
                    help='GPU tag')
args = parser.parse_args()











os.environ['CUDA_VISIBLE_DEVICES'] = '0'

input_img_directory = '/media/SteinLab5/WT11L/output/stitched'
model_file = os.path.join(os.environ['PYTHONPATH'], '/nuclei/models/isensee_2017_model2.h5')

mask_file = '/media/SteinLab5/WT11L/output/variables/I_mask.mat'
save_name = 'WT11L_centroids2.csv'
save_name = os.path.join(os.path.split(input_img_directory)[0], save_name)


pred_threshold = 0.5  # Prediction threshold. Default is 0.5
int_threshold = 400  # Minimum intensity of positive cells. Otherwise set to None
apply_opening = False  # Postprocess using morphological opening
measure_coloc = True  # Measure intensity of colocalizaed channels
normalize_intensity = True  # Whether to normalize intensities
n_channels = 3         # Total number of channels

# Image and model configuration parameters
# Dimensions should correspond to [x,y,z]
chunk_size = [112, 112, 32]  # Model chunk size in pixels
overlap = [16, 16, 8]  # Overlap between chunks in pixels

acquired_img_resolution = [1.208, 1.208, 4]  # Resolution of acquired images in um/pixel
trained_img_resolution = [1.208, 1.208, 4]  # Resolution the model was trained in um/pixel
mask_resolution = [10, 10, 10]  # Resolution of mask in um/pixel

#####
# Load model and images
model = load_old_model(model_file)
mask = loadmat(mask_file)
mask = mask['I_mask']

# Determine chunk sizes while considering scaling and chunk overlap
res = [i / j for i, j in zip(acquired_img_resolution, trained_img_resolution)]
mask_res = [i / j for i, j in zip(acquired_img_resolution, mask_resolution)]

load_chunk_size = [i / j for i, j in zip(chunk_size, res)]
load_chunk_size = [round(s) for s in load_chunk_size]

# Taking only channel 1 images (assumed to be cell nuclei)
files = os.listdir(input_img_directory)

# Take img_list for each channel. Nuclei is channel should be first
img_list = []
for i in range(n_channels):
    matching = [s for s in files if "C" + str(i+1) in s]
    matching = sorted(matching)
    # Read only the number of images for each chunk
    img_list.append([input_img_directory + '/' + s for s in matching])

total_slices = len(img_list[0])

# z_pad = np.ceil(total_images / load_chunk_size[2])
# n_chunks = np.ceil((total_images + z_pad) / load_chunk_size[2]).astype(int)

# chunk_start = np.arange(0, n_chunks * (load_chunk_size[2] - overlap[2]), load_chunk_size[2] - overlap[2])
# chunk_end = chunk_start[1:] + overlap[2]
# chunk_end = np.append(chunk_end, total_images - 1)

chunk_start = np.arange(0, total_slices, step=(chunk_size[2] - overlap[2]))
chunk_end = chunk_start + chunk_size[2]
n_chunks = len(chunk_start)

# Resize mask to match the number of slices in input image directory
mask = resize(mask, (mask.shape[0], mask.shape[1], total_slices), order=0)
mask_res[-1] = 1

# Create 3D weight patch
chunk_weight_patch = create_3d_weight_patch(chunk_size, overlap)

# Calculate rescaling intensity values
if normalize_intensity:
    intensity_values = calculate_rescaling_intensities(img_list[0], sampling_range=10)
    if int_threshold is not None:
        int_threshold = rescale_intensity(np.asarray(int_threshold, dtype=np.float32),
                                          in_range=intensity_values)

# Calculate rescaling factor if resolutions are different
if acquired_img_resolution != trained_img_resolution:
    rescale_factor = [i / j for i, j in zip(acquired_img_resolution, trained_img_resolution)]
else:
    rescale_factor = None

# Read first image to get sizes
tempI = cv2.imread(img_list[0][0], -1)
[rows, cols] = tempI.shape

## Take mask index 84 for testing
# mask = mask == 84
prediction_save = np.array([])
x_shift = 0
y_shift = 0
z_shift = 0
total_cells = 0
total_time = datetime.now()

# Begin cell counting
for n in range(n_chunks):
    startTime = datetime.now()
    z_start = chunk_start[n]
    z_end = chunk_end[n]

    # Skip chunk if nothing is present
    if not mask[:, :, z_start:z_end].any():
        continue

    # Take the mask for the respective z positions
    mask_chunk1 = [cv2.resize(mask[:, :, z], (cols, rows), interpolation=0) for z in range(z_start, z_end)]
    mask_chunk1 = np.asarray(mask_chunk1, dtype=bool)
    mask_chunk = np.swapaxes(mask_chunk1, 0, 1).swapaxes(1, 2)

    # Read images
    print('Working on chunk', n + 1, 'out of', n_chunks)
    print('Reading slices', z_start, 'through', z_end)
    images = [cv2.imread(file, -1) for file in img_list[0][z_start:z_end]]
    images = np.asarray(images, dtype=np.float32)
    images = np.swapaxes(images, 0, 1).swapaxes(1, 2)

    # If last chunk, then pad bottom to make it fit into the model
    if z_end > total_slices:
        print('Padding Last Chunk')
        end_pad = chunk_size[2] - images.shape[2]
        images = np.pad(images, ((0, 0), (0, 0), (0, end_pad)), 'mean')
        z_end = total_slices

    # Rescale intensity
    print('Rescaling Intensity...')
    if normalize_intensity:
        images_rescaled = rescale_intensity(images, in_range=intensity_values)
    else:
        images_rescaled = images

    # Rescale image size
    print('Rescaling Size...')
    if rescale_factor is not None:
        images_rescaled = ndimage.zoom(images_rescaled, rescale_factor, order=1)

    [rows, cols, slices] = images_rescaled.shape

    # Pad mask if last slice
    if z_end > total_slices:
        end_pad = chunk_size[2] - images.shape[2]
        mask_chunk = np.pad(mask_chunk, ((0, 0), (0, 0), (0, end_pad)), 'constant')

    # Calculate y and x positions to sample image chunks
    x_positions = np.arange(0, cols, step=(chunk_size[0] - overlap[0]))
    y_positions = np.arange(0, rows, step=(chunk_size[1] - overlap[1]))

    n_chunks_c = len(x_positions)
    n_chunks_r = len(y_positions)

    # Append image and mask chunks to list
    img_chunks = []
    msk_chunks = []
    for i in range(n_chunks_r):
        for j in range(n_chunks_c):
            msk_chunk = mask_chunk[y_positions[i]:y_positions[i] + chunk_size[0],
                        x_positions[j]:x_positions[j] + chunk_size[1], :]
            img_chunk = images_rescaled[y_positions[i]:y_positions[i] + chunk_size[0],
                        x_positions[j]:x_positions[j] + chunk_size[1], :]

            if img_chunk.shape != tuple(load_chunk_size):
                pad_c_right = int(load_chunk_size[1] - img_chunk.shape[1])
                pad_r_bottom = int(load_chunk_size[0] - img_chunk.shape[0])

                msk_chunk = np.pad(msk_chunk, ((0, pad_r_bottom), (0, pad_c_right), (0, 0)), 'constant')
                img_chunk = np.pad(img_chunk, ((0, pad_r_bottom), (0, pad_c_right), (0, 0)), 'mean')

            msk_chunks.append(msk_chunk)
            img_chunks.append(img_chunk)

    print('Images prepared in: ', datetime.now() - startTime)

    # Run prediction
    # The input shape should be 5 dimensions: (m, n, x, y, z)
    # x, y, z represent the image shape, as you would expect. n is the number of
    # channels. In a standard color video image, you would have 3 channels (red,
    # green, blue). In medical imaging these channels can be separate imaging
    # modalities. m is the batch size or number of samples being passed to the
    # model for training.
    startTime = datetime.now()
    output = []
    img_shape = (1, 1) + tuple(chunk_size)
    empty_chunk = np.zeros(img_shape)
    empty_idx = np.zeros(len(img_chunks))
    img_reshaped = []
    msk_reshaped = []
    for idx in range(len(img_chunks)):
        img_reshaped.append(np.reshape(img_chunks[idx], img_shape))
        msk_reshaped.append(np.reshape(msk_chunks[idx], img_shape))
        if msk_reshaped[idx].any():
            output.append(model.predict(img_reshaped[idx]) * msk_reshaped[idx])
        else:
            output.append(empty_chunk)
            empty_idx[idx] = 1


    output = [np.squeeze(chunk) for chunk in output]
    print('Mask prediction time elapsed: ', datetime.now() - startTime)

    # Reassemble prediction image
    startTime = datetime.now()
    a = 0
    cen = []
    for i in range(n_chunks_r):
        for j in range(n_chunks_c):
            if empty_idx[a] != 1:
                # Calculate connected components and determine centroid positions
                # Calculate final prediction mask and label connected components
                output_thresh = np.where(output[a] > pred_threshold, 1, 0)

                labels_out = cc3d.connected_components(output_thresh)
                n_cells = np.max(labels_out)

                # Find centroids
                centroids = ndimage.measurements.center_of_mass(output_thresh, labels_out,
                                                                index=np.arange(1, n_cells + 1))
                centroids = np.asarray(centroids).round()

                # Remove cells with low intensity
                if int_threshold is not None:
                    img_chunk = img_chunks[a]
                    int_high = [img_chunk[tuple(centroids[c].astype(int))] > int_threshold for c in
                                range(len(centroids))]
                    centroids = centroids[int_high]

                if centroids.any():
                    # Remove centroids along borders
                    centroids = centroids[centroids[:, 0] > overlap[0] / 2]
                    centroids = centroids[centroids[:, 0] <= chunk_size[0] - overlap[0] / 2]

                    centroids = centroids[centroids[:, 1] > overlap[1] / 2]
                    centroids = centroids[centroids[:, 1] <= chunk_size[1] - overlap[1] / 2]

                    centroids = centroids[centroids[:, 2] > overlap[2] / 2]
                    centroids = centroids[centroids[:, 2] <= chunk_size[2] - overlap[2] / 2]

                    # Adjust centroid positions
                    centroids[:, 0] += y_positions[i]
                    centroids[:, 1] += x_positions[j]
                    centroids[:, 2] += chunk_start[n]

                    cen.append(centroids)
            a += 1

    if not cen:
        continue

    cent = np.concatenate(cen)
    cent = cent[cent[:, 2].argsort()]

    print('Postprocessing time elapsed: ', datetime.now() - total_time)
    print('Cells counted: ', cent.shape[0])

    total_cells += cent.shape[0]
    print('Total cells counted: ', total_cells)

    print('Measuring colocalization...')

    # Get mask structure id's
    structure_idx = [mask[tuple(np.round(cent[c]*mask_res).astype(int))] for c, cents in enumerate(cent)]
    cent = np.append(cent, np.array(structure_idx)[:, None], axis=1)

    # Remove stray cells
    rm_idx = cent[:, 3] == 0
    cent = cent[~rm_idx, :]
    total_cells += -sum(rm_idx)
    print('Removed ' + str(sum(rm_idx)) + ' empty nuclei')


    if measure_coloc and cent.any():
        cent[:, 2] -= z_start
        for i in range(n_channels):
            intensities = measure_colocalization(cent, img_list[i][z_start:z_end])

            ## Throwing errors
            cent = np.append(cent, intensities[:, None], axis=1)
        cent[:, 2] += z_start

    if n == 0:
        np.savetxt(save_name, cent.round().astype(int), delimiter=",", fmt='%u')
    else:
        with open(save_name, "ab") as f:
            np.savetxt(f, cent.round().astype(int), delimiter=",", fmt='%u')

print('Total cells counted: ', total_cells)
print('Total time elapsed: ', datetime.now() - total_time)
