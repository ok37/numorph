function I_mask = gen_mask(hemisphere, structures_of_interest, res_adj)
% -------------------------------------------------------------------------
%This function generates a mask from Allen atlas based on CCF ids in a csv
%file
% -------------------------------------------------------------------------

% Load annotation volume and indexes
annotation_path = fullfile('.','supplementary_data','annotationData.mat');
load(annotation_path,'annotationVolume')
annotationVolume = permute(annotationVolume, [3,1,2]);
if isequal(hemisphere, 'left')
    annotationVolume = flip(annotationVolume,2);
end

% Adjust resolution if requested
if nargin > 2
    new_size = size(annotationVolume).*res_adj;
    annotationVolume = imresize3(annotationVolume,new_size,'Method','nearest');
end

% Use all or some structures
if nargin < 2
    path_to_id = './supplementary_data/structure_template.csv';
    id = readtable(path_to_id);
else
    for i = 1:length(structures_of_interest)
        file = sprintf("%s.csv",structures_of_interest(i));
        path_to_id = fullfile('.','annotations',file);
        if i == 1
            id = readtable(path_to_id);
        else
            id = vertcat(id,readtable(path_to_id));
        end
    end
    
    % Read csv file containing region ids
    id = id.index;

    % Get unique annotations in the volume and their indexes
    [C, ~, ic] = unique(annotationVolume(:));

    % Set structures not found in the table to 0
    C(~ismember(C,id)) = 0;
    
    % Reshape linearized matrix back to its original size
    I_mask = reshape(C(ic),size(annotationVolume));
end

end