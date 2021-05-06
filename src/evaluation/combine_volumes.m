function df_results = combine_volumes(config)
%--------------------------------------------------------------------------
% Measure structure volumes and calculate cell densities if counts are 
% present. Append to full multi-sample dataset for all annotated
% structures.
%--------------------------------------------------------------------------
vol_results = fullfile(config.results_directory,strcat(config.prefix,"_summary_volumes.csv"));
temp_file = fullfile(fileparts('NM_config'),'annotations','structure_template.csv');

% Read annotation indexes
df_temp = readtable(temp_file);
annotation_indexes = df_temp.index;



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





% Create matrix for storing volumes
df_volumes = zeros(length(annotation_indexes),length(sample_files));
for i = 1:length(config.results_path)
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