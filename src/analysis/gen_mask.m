function I_mask = gen_mask(hemisphere, structures, res_adj)
% -------------------------------------------------------------------------
%This function generates a mask from Allen atlas based on CCF ids in a csv
%file
% -------------------------------------------------------------------------

% Load annotation volume and indexes
home_path = fileparts(which('NM_config'));
annotation_path = fullfile(home_path,'data','annotationData.mat');
load(annotation_path,'annotationVolume')

% Default annotations in .mat file are in 'sra' orientation
% Permute to default orientation for registration which is 'air'
%annotationVolume = permute_orientation(annotationVolume, 'sra', 'air');

% Use if no structures provided
if nargin < 2 || isequal(structures,"structure_template.csv") ||...
        isequal(structures,"structure_template") ||...
        isempty(structures)
    I_mask = annotationVolume;
    return
end

% Adjust resolution if specifed
if nargin > 2
    new_size = size(annotationVolume).*res_adj;
    annotationVolume = imresize3(annotationVolume,new_size,'Method','nearest');
end

% Read structure ids
files = dir(fullfile('annotations/**'));

% Get bw structure masks for major structures
mat_files = files(arrayfun(@(s) endsWith(s.name,'.mat'),files));
major_structures = arrayfun(@(s) string(s.name(1:end-4)),mat_files);
idx = ismember(major_structures,structures);
if any(idx)
    if sum(idx) < length(structures)
        error("Structure names not specified correctly")
    end
    
    % Load binary mask
    bw = false(size(annotationVolume));
    s = mat_files(idx);
    for i = 1:length(s)
       load(fullfile(s(i).folder,s(i).name),'bw_mask')
       bw = bw | bw_mask;
    end
    
    % Apply maks and return
    annotationVolume(~bw) = 0;
    I_mask = annotationVolume;
    return
end

% Custom annotations, this may take a while
id = cell(1,length(structures));
for i = 1:length(structures)
    if endsWith(structures,'.csv')
        file = structures(i);
    else
        file = sprintf("%s.csv",structures(i));
    end
    
    idx = files(arrayfun(@(s) s.name == file,files));
    if length(idx)>1
        error("Multiple annotation files with the same name in the annotations directory")
    elseif isempty(idx)
        error("Could not locate annotation file")
    end
    id{i} = readtable(fullfile(idx.folder,idx.name));
end
id = cat(1,id{:});

% Read csv file containing region ids
id = id.index;

% Get unique annotations in the volume and their indexes
[C, ~, ic] = unique(annotationVolume(:));

% Set structures not found in the table to 0
C(~ismember(C,id)) = 0;

% Reshape linearized matrix back to its original size
I_mask = reshape(C(ic),size(annotationVolume));

end