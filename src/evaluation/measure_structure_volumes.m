function measure_structure_volumes(config)
% Take a registered annotation volume and calculate structure voxel
% volumes. Input is a string array of fullfile locations

% Calculate volume per voxel in mm^3
resolution = repmat(config.resample_resolution,1,3);
mm = prod(resolution)/(1E3^3);

% Read structure template file
df_template = readtable(fullfile(fileparts(which('NM_config')),'annotations','structure_template.csv'));
indexes = df_template.index;

% Measure number of voxels for each structure
load(fullfile(config.output_directory,'variables',strcat(config.sample_id,'_mask.mat')),'I_mask')
total_volume = sum(I_mask(:)>0)*mm;

sums = zeros(1,length(indexes));
counted = histcounts(I_mask(:),'BinMethod','integers');
sums(1:length(counted)) = counted;
df_volumes = sums*mm';

% Set top row equal to 0 as this is the background
df_volumes(1) = 0;

% Sum according to structure level order, except for background
ids = df_template.id;
path = df_template.structure_id_path;
df_new = zeros(size(df_volumes));
for i = 2:length(ids)
    idx = cellfun(@(s) contains(s,string("/"+ids(i)+"/")), path);
    df_new(i) = sum(df_volumes(idx));
end
df_volumes = df_new';

% Create header name
df_header = config.sample_id + "_" + "Volume";

% Convert to table
volumes_table = array2table(df_volumes,'VariableNames',df_header);

% Concatenate counts to results table
df_results = horzcat(df_template(:,[1,end-1,end]), volumes_table);

% Write new results file
save_name = fullfile(config.output_directory,strcat(config.sample_id,'_volumes.csv'));
writetable(df_results,save_name)
fprintf('%s\t Total structure volume: %.2f mm^3\n',datetime('now'),total_volume)

end