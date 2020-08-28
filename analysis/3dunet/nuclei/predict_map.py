import os

from train_isensee2017 import config
from unet3d.training import load_old_model
from unet3d.prediction import run_validation_cases

def main():
    prediction_dir = os.path.abspath("prediction")

    #model_file = '/home/ok37/repos/3dunet-centroid/nuclei/isensee_2017_model.h5'
    #model = load_old_model(model_file)

    #img_data = np.asarray(img_chunks, dtype=np.float32)
    #img_shape = [1, 1] + chunk_size

    #startTime = datetime.now()

    #output = []
    #empty = 0nib.Nifti1Image(data, affine)
    #for idx in range(len(img_data)):
     #   img_reshaped = np.reshape(img_data[idx], img_shape)
     #   if msk_chunks[idx].sum() > 0:
     #       output.append(model.predict(img_reshaped))
     #   else:
     #       output.append([])
     #       empty = empty + 1

    run_validation_cases(validation_keys_file=config["validation_file"],
                         model_file=config["model_file"],
                         training_modalities=config["training_modalities"],
                         labels=None,
                         hdf5_file=config["data_file"],
                         output_label_map=True,
                         output_dir=prediction_dir)


if __name__ == "__main__":
    main()
