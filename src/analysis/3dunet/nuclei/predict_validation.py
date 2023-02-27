"""
Run validation cases in ./data/validation

"""
import os
import nibabel as nib
import numpy as np
import argparse
from unet3d.training import load_old_model
from unet3d.prediction import patch_wise_prediction
from nuclei.img_utils import patch_wise_prediction2, prediction_to_centroids

parser = argparse.ArgumentParser(description='Predict from validation dataset.')
parser.add_argument('--r', metavar='m', type=str, nargs='+',
                    help='Prediction resolution tag')
parser.add_argument('--g', metavar='g', type=str, nargs='+',
                    help='GPU tag')
parser.add_argument('--m', metavar='m', type=str, nargs='+',
                    help='Model name')
args = parser.parse_args()

if args.r:
    resolution = args.r[0]
else:
    resolution = '075'

if args.m:
    model_name = args.m[0]
else:
    model_name = '075_121_model.h5'
model_name = model_name + '.h5' if not model_name.endswith('.h5') else model_name

if args.g:
    gpu_idx = args.g[0]
    os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_idx)
else:
    os.environ['CUDA_VISIBLE_DEVICES'] = '0'

model_path = os.path.join(os.getcwd(), 'models', model_name)

# %% Run
def main(resolution, model_path):
    model = load_old_model(model_path)

    # Get images
    print("Loading validation data with fully segmented images")
    images = get_validation_folders(resolution)

    # Run prediction

    print("Running prediction using model " + os.path.basename(model_path))
    predictions = []
    for i in range(len(images)):
        data = images[i].get_fdata()
        data = np.squeeze(data)
        data = np.asarray(data, dtype=np.float64)

        patch_shape = tuple([int(dim) for dim in model.input.shape[-3:]])

        if patch_shape == data.shape[-3:]:
            img_shape = (1, 1) + data.shape
            data_reshaped = np.reshape(data, img_shape)
            prediction = model.predict(data_reshaped)
        else:
            img_shape = (1, 1) + data.shape
            data_reshaped = np.reshape(data, img_shape)
            prediction = patch_wise_prediction2(model=model, data=data_reshaped, overlap=[16, 16, 8],
                                                pred_threshold=0.5, int_threshold=0)

        predictions.append(np.squeeze(prediction))

    # Save results
    save_directory = os.path.join(os.getcwd(), 'predictions')

    if not os.path.isdir(save_directory):
        os.mkdir(save_directory)

    save_prediction_results(images, predictions, save_directory)

    return


# %% Read Images
def get_validation_folders(resolution):
    img_path = os.path.join(os.getcwd(), 'data', 'validation', resolution)
    files = os.listdir(img_path)
    files = [f for f in files if f[0] == "f"]

    images = []
    for f in files:
        images.append(nib.load(os.path.join(img_path, f)))

    return images


# %% Save Results
def save_prediction_results(images, predictions, save_directory):
    for i in range(len(images)):
        save_folder_name = os.path.join(save_directory, 'f' + str(i + 1))

        if not os.path.isdir(save_folder_name):
            os.mkdir(save_folder_name)

        save_name = os.path.join(save_folder_name, 'image.nii')
        nib.save(nib.Nifti1Image(images[i].get_data().astype('uint16'), affine=np.eye(4)), save_name)

        save_name = os.path.join(save_folder_name, 'prediction.nii')
        nib.save(nib.Nifti1Image(predictions[i].astype('uint8'), affine=np.eye(4)), save_name)

    return


# %% Predict individual case
if __name__ == "__main__":
    main(resolution, model_path)
