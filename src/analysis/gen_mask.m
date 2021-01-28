function I_mask = gen_mask(hemisphere, structures_of_interest, res_adj)
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
annotationVolume = permute_orientation(annotationVolume, 'sra', 'air');

% Use if no structures provided
if nargin < 2
    I_mask = annotationVolume;
    return
end

% Adjust resolution if specifed
if nargin > 2
    new_size = size(annotationVolume).*res_adj;
    annotationVolume = imresize3(annotationVolume,new_size,'Method','nearest');
end


% Read structure ids
id = cell(1,length(structures_of_interest));
for i = 1:length(structures_of_interest)
    if endsWith(structures_of_interest,'.csv')
        file = structures_of_interest(i);
    else
        file = sprintf("%s.csv",structures_of_interest(i));
    end
    path_to_id = fullfile(home_path,'annotations',file);
    id{i} = readtable(path_to_id);
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