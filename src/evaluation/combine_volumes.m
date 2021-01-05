function df_results = combine_volumes(files, results_path, sum_child_structures)
% Take a registered annotation volume and calculate structure voxel
% volumes. Input is a string array of fullfile locations
resolution = [25,25,25]; % Assumes 25um/voxel resolution
results_path = fullfile(results_path,"TCe_summary_volumes.csv");
new_template_path = './supplementary_data/structure_template.csv';

% Read annotation indexes
df_template = readtable(new_template_path);
annotation_indexes = df_template.index;

% Calculate volume per voxel in mm^3
mm = prod(resolution)/(1E3^3);

% Create matrix for storing volumes
df_volumes = zeros(length(files),length(annotation_indexes));
sample_names = strings(1,length(files));

% Iterate over all f
for i = 1:length(files)
    [~,s] = fileparts(files(i));
    s = strsplit(s,'_');
    sample_names(i) = s(1);
    fprintf(strcat(char(datetime('now')),"\t Loading %s...\n"), sample_names(i));

    load(files(i),'I_mask')
    for j = 1:length(annotation_indexes)
        df_volumes(i,j) = sum(I_mask(:)==annotation_indexes(j))*mm;
    end
end

df_volumes = df_volumes';

% Set top row equal to 0 as this is the background
df_volumes(1,:) = 0;

if isequal(sum_child_structures,'true')
    % Sum according to structure level order, except for background
    ids = df_template.id;
    path = df_template.structure_id_path;
    df_new = zeros(size(df_volumes));
    for i = 2:length(ids)
        idx = cellfun(@(s) contains(s,string("/"+ids(i)+"/")), path);
        df_new(i,:) = sum(df_volumes(idx,:),1);
    end
    df_volumes = df_new;
end

% Create header name
df_header = sample_names + "_" + "Volume";

% Convert to table
volumes_table = array2table(df_volumes,'VariableNames',df_header);

% Concatenate counts to results table
df_results = horzcat(df_template, volumes_table);

% Write new results file
writetable(df_results,results_path)
end