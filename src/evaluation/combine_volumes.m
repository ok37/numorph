function df_results = combine_volumes(sample_files,config)
% Take a registered annotation volume and calculate structure voxel
% volumes. Input is a string array of fullfile locations

vol_results = fullfile(config.results_directory,"NMe_summary_volumes.csv");
temp_file = fullfile(fileparts('NM_config'),'annotations','structure_template.csv');

% Read annotation indexes
df_temp = readtable(temp_file);
annotation_indexes = df_temp.index;

% Create matrix for storing volumes
df_volumes = zeros(length(annotation_indexes),length(sample_files));
for i = 1:length(sample_files)
    tbl = readtable(sample_files(i));
    df_volumes(:,i) = tbl{:,end};
end

% Set top row equal to 0 as this is the background
df_volumes(1,:) = 0;

% Create header name
groups = cat(1,config.groups{:});
df_header = config.samples' + "_" + groups(:,2)'+ "_Volume";

% Convert to table
volumes_table = array2table(df_volumes,'VariableNames',df_header);

% Concatenate counts to results table
df_results = horzcat(df_temp, volumes_table);

% Write new results file
writetable(df_results,vol_results)

end