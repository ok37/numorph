function generate_annotations_from_file(config, filepath)
%--------------------------------------------------------------------------
% Generate I_mask.mat from a specified custom annotation file.
%--------------------------------------------------------------------------
% Usage:
% generate_annotations_from_file(filepath)
%
%--------------------------------------------------------------------------
% Inputs:
% config: Configuration structure from analysis stage.
%
% filepath: (string, char) Path to annotation file. (default: leave empty,
% read from config
%
%--------------------------------------------------------------------------

if nargin<2
    filepath = config.annotation_file;
end

% Check if annotations exist
if ~isfile(filepath)
    error("Could not find custom annotation mask %s",filepath)
end

% Read image
try
    I_mask = read_img(filepath);
catch 
    error("Error reading annotation mask file")
end

% Resize to resample resolution 
if ~isempty(config.annotation_resolution) &&...
        config.annotation_resolution ~= config.resample_resolution
    res_adj = config.annotation_resolution/config.resample_resolution;
    I_mask = imresize3(I_mask,res_adj,'Method','nearest');
end

% Permute orientation
%I_mask = permute_orientation(I_mask,config.annotation_orientation,'psl');

% Create .mat file
save_file = fullfile(config.output_directory,'variables',strcat(config.sample_id,'_mask.mat'));
save(save_file,'I_mask','-v7.3')

end