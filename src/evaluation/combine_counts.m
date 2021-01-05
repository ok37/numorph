function df_results = combine_counts(sample_files,config)
%--------------------------------------------------------------------------
% Measure cell count results and append to full multi-sample dataset for
% all annotated structures
%--------------------------------------------------------------------------
config.results_directory = fullfile(config.results_directory,"TCe_summary_counts.csv");
new_template_path = './supplementary_data/structure_template.csv';
types = config.markers;

% Read template or check for previous results to append to
if ~isequal(config.overwrite,'true')
    % Check for existing results to append to
    if exist(config.results_directory,'file') == 0
       fprintf(strcat(char(datetime('now')),"\t Cannot locate previous "+...
       "results. Creating new cell count summary file\n"))
       df_results = readtable(new_template_path);
    else
       df_results = readtable(config.results_directory);
    end
else
   df_results = readtable(new_template_path);
end

% Annotation indexes
indexes = df_results.index;

for n = 1:length(sample_files)
    % Choose sample
    s_idx = find(arrayfun(@(s) contains(sample_files(n),s),config.sample_names));
    fprintf(strcat(char(datetime('now')),"\t Reading sample %s...\n"),config.sample_names(s_idx));
    assert(length(s_idx) == 1, "Sample name doesn't match just .csv file "+...
        "in sample_directory")

    % Read table
    df_sample = readmatrix(sample_files{n});
    ct = df_sample(:,end);

    % Make sure number of cell-types present equals the number of config.markers
    assert(length(unique(ct)) == length(types),"Number of specified cell types (%d) "+...
        "does not match unique cell type indexes (%d) in sample file.\n",length(unique(ct)),...
        length(types))
    
    % Create count matrix
    counts = zeros(height(df_results), length(types));

    % Apply thresholds and count cells in each structure for each channel
    % First channel cell intensities start at column 5
    for i = 1:length(indexes)
        df_sub = df_sample(df_sample(:,4) == indexes(i),end);
        for j = 1:length(types)
            counts(i,j) = sum(df_sub == j);
        end

    end

    if isequal(config.sum_child_structures,'true')
        % Sum according to structure level order, except for background
        ids = df_results.id;
        path = df_results.structure_id_path;
        for i = 2:length(ids)
            idx = cellfun(@(s) contains(s,string("/"+ids(i)+"/")), path);
            counts(i,:) = sum(counts(idx,:),1);
        end
    end

    % Create header name
    df_header = config.sample_names(s_idx) + "_" + config.markers + "_Counts";

    % Convert to table
    counts_table = array2table(counts,'VariableNames',df_header);

    % Concatenate counts to results table
    df_results = horzcat(df_results, counts_table);
    
    % Print cells counted
    for j = 1:length(types)
        fprintf(strcat(char(datetime('now')),"\t Total %s Count: %d\n"),...
            config.markers(j),counts(2,j));
    end
end

% Write new results file
writetable(df_results,config.results_directory)

end
