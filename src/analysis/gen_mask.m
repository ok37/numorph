function I_mask = gen_mask(structures, hemisphere, orientation, res_adj)
% -------------------------------------------------------------------------
%This function generates a mask from Allen atlas based on CCF ids in a csv
%file
% -------------------------------------------------------------------------

if nargin<4
    res_adj = [];
end

% Load annotation volume and indexes
home_path = fileparts(which('NM_config'));
annotation_path = fullfile(home_path,'data','annotation_data');
load(fullfile(annotation_path,'annotationData.mat'),'annotationVolume')
annotationVolume = reshape_mask(annotationVolume,hemisphere,orientation,res_adj);

% No structures provided
% Return whole brain mask
if nargin < 2 || isequal(structures,"structure_template.csv") ||...
        isequal(structures,"structure_template") ||...
        isempty(structures)
    I_mask = annotationVolume;
    return
end

% Check for default structure
% These have have masks already precomputed
mat_files = dir(fullfile(annotation_path,'*.mat'));
major_structures = arrayfun(@(s) string(s.name(1:end-4)),mat_files);
[~,fname] = fileparts(structures);
idx = ismember(major_structures,fname);
if any(idx)
    if sum(idx) < length(structures)
        error("Structure names not specified correctly")
    end
    
    % Load binary mask
    bw = false(size(annotationVolume));
    s = mat_files(idx);
    for i = 1:length(s)
       load(fullfile(s(i).folder,s(i).name),'bw_mask')
       bw_mask = reshape_mask(bw_mask,hemisphere,orientation,res_adj);
       bw = bw | bw_mask;
    end
    
    % Apply maks and return
    annotationVolume(~bw) = 0;
    I_mask = annotationVolume;
    return
end

% Check for custom structures
% Generate and apply mask, this will take a while
files = dir(fullfile(home_path,'annotations','custom_structures'));
id = cell(1,length(structures));
for i = 1:length(structures)    
    idx = files(arrayfun(@(s) s.name == structures(i),files));
    if length(idx)>1
        error("Multiple annotation files with the same name in the annotations directory")
    elseif isempty(idx)
        error("Could not locate annotation file")
    end
    id{i} = readtable(fullfile(idx.folder,idx.name));
end
id = cat(1,id{:});

% Keep only unique rows
id = unique(id,'rows');

% Read csv file containing region ids
id = id.index;

% Get unique annotations in the volume and their indexes
[C, ~, ic] = unique(annotationVolume(:));

% Set structures not found in the table to 0
C(~ismember(C,id)) = 0;

% Reshape linearized matrix back to its original size
I_mask = reshape(C(ic),size(annotationVolume));

end


function annotationVolume = reshape_mask(annotationVolume, hemisphere, orientation, res_adj)

% Check hemisphere shape, orientation
if isequal(hemisphere,"both")
    img = flip(annotationVolume,2);
    annotationVolume = cat(2,annotationVolume,img);
end
    
% Permute
annotationVolume = permute_orientation(annotationVolume,'sra',char(orientation));

% Adjust resolution if specifed
if ~isempty(res_adj)
    new_size = size(annotationVolume).*res_adj;
    annotationVolume = imresize3(annotationVolume,new_size,'Method','nearest');
end

end
