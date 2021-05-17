function munge_atlas(atlas_file, annotation_file, resolution, orientation, hemisphere, out_resolution)
%--------------------------------------------------------------------------
% Munge custom atlas files and associated annotations and save to /data for
% use. 
%--------------------------------------------------------------------------
% Usage:
% munge_atlas(atlas_file, annotation_file, resolution, orientation,
% hemisphere, out_resolution)
%
%--------------------------------------------------------------------------
% Inputs:
% atlas_file: (string) Path to atlas file.
%
% annotation_file: (string) Path to associated annotations.
%
% resolution: (1x3 numeric) Atlas y,x,z resolution specified as micron per
% voxel.
%
% orientation: (1x3 char) Atlas orientation (a/p,s/i,l/r).
%
% hemisphere: ("left","right","both","none") Which brain hemisphere.
%
% out_resolution: (int) Isotropic resolution of atlas output. (default: 25)
%
%--------------------------------------------------------------------------

if nargin<6
    out_resolution = 25;
end

% Read images
img = read_img(atlas_file);
annotations = read_img(annotation_file);

% Standardize
img = standardize_nii(img, "raw", resolution, orientation, hemisphere,...
    out_resolution, 'ail', hemisphere, 'uint16');

annotations = standardize_nii(annotations, "mask", resolution,...
    orientation, hemisphere, out_resolution, 'ail', hemisphere, 'uint16');

annotationData.annotationVolume = annotations;
annotationData.annotationIndexes = unique(annotationData.annotationVolume);
annotationData.hemisphere = hemisphere;
annotationData.resolution = out_resolution;

% Save files
home_path = fileparts(which('NM_config'));
[~,filename] = fileparts(atlas_file);
filepath = fullfile(home_path,'data','atlas',strcat(filename,'.nii'));
niftiwrite(img,filepath)

filepath = fullfile(home_path,'data','annotation_data',strcat(filename,'.mat'));
save(filepath,'-struct','annotationData')

end