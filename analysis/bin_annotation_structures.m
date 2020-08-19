function [output,index] = bin_annotation_structures(input, struct_table, keep_ids)
%-------------------------------------------------------------------------
% This function bins substructure annotations in a mask to a list of
% larger lower level structures defined .csv file structure_path
%-------------------------------------------------------------------------

% Default read harris cortical groupings
if nargin<2 || isequal(struct_table,'cortex')
    struct_table = readtable(fullfile('annotations','harris_cortical_groupings.csv'));
elseif isequal(struct_table,'cortex_large')
    struct_table = readtable(fullfile('annotations','harris_cortical_groupings_large.xls'));
elseif isequal(struct_table,'layers')
    [output,index] = bin_annotation_layers(input);
    return
end

% Whether to keep ids that are not binned
if nargin<3
    keep_ids = 'false';
end

% Get full list of structures
full = readtable(fullfile('supplementary_data','structure_template.csv'));

% Get structure tree paths
paths = cellfun(@(s) strsplit(s,'/'),full.structure_id_path,...
        'UniformOutput',false);

% Get structures of interest
index = struct_table.index;
csv_ids = struct_table.id;

% Get unique annotations in the volume and their indexes
[C, ~, ic] = unique(input(:));

% Get ARA ids for these indexes
ic_ids = full.id(ismember(full.index,C));

% Get new indexes for annotation
C_new = zeros(length(C),1);
for i = 1:length(ic_ids)
    % Find position where id equals ARA id
    idx = full.index(full.id == ic_ids(i))+1;
    % Get structure tree at this position
    p = paths{idx}(cellfun(@(s) ~isempty(s), paths{idx}));
    % Vectorize
    p = double(cellfun(@(s) string(s),p));
    % Get index from structures of interest where id is present
    idx_new = index(ismember(csv_ids,p));
    if ~isempty(idx_new)
        C_new(i) = idx_new;
    elseif isequal(keep_ids,'true')
        C_new(i) = C(i);
    end
end

% Reshape new annotation using new indexes
output = reshape(C_new(ic),size(input));

end