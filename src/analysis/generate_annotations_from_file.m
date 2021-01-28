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
% read from 
%
%--------------------------------------------------------------------------

home_path = fileparts(which('NM_config'));
annot_file = fullfile(home_path,'data','masks',config.use_annotation_mask);

% Check if annotations exist
if ~isfile(annot_file)
    error("Could not find custom annotation mask %s",annot_file)
end

% Read image
try
    I_mask = read_img(annot_file);
catch 
    error("Error reading annotation mask file")
end

% Resize to resample resolution 
if config.annotation_resolution ~= config.resample_resolution
    res_adj = config.annotation_resolution/config.resample_resolution;
    I_mask = imresize3(I_mask,res_adj,'Method','nearest');
end

% Permute orientation
I_mask = permute_orientation(I_mask,config.annotation_orientation,'psl');

% Create .mat file
[a,b] = fileparts(annot_file);
save(fullfile(a,b),'I_mask')

end