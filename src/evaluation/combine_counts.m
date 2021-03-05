function df_results = combine_counts(sample_files,config)
%--------------------------------------------------------------------------
% Measure cell count results and append to full multi-sample dataset for
% all annotated structures
%--------------------------------------------------------------------------
count_results = fullfile(config.results_directory,"NMe_summary_counts.csv");
temp_file = fullfile(fileparts('NM_config'),'annotations','structure_template.csv');

if isequal(config.use_classes,"true")
    classes = config.keep_classes;
    if isempty(classes)
        error("Please specify which classes to count as 'keep_classes'")
    else
        types = config.class_names(classes);
    end
else
    types = "nuclei";
    classes = 1;
end

% Read template or check for previous results to append to
if ~isequal(config.overwrite,'true')
    % Check for existing results to append to
    if exist(count_results,'file') == 0
       fprintf(strcat(char(datetime('now')),"\t Cannot locate previous "+...
       "results. Creating new cell count summary file\n"))
       df_results = readtable(temp_file);
    else
       df_results = readtable(count_results);
    end
else
   df_results = readtable(temp_file);
end

% Annotation indexes
indexes = df_results.index;

% Sum cell counts for each sample
for n = 1:length(sample_files)
    % Choose sample
    s_idx = find(arrayfun(@(s) contains(sample_files(n),s),config.samples));
    fprintf(strcat(char(datetime('now')),"\t Reading sample %s...\n"),config.samples(s_idx));
    assert(length(s_idx) == 1, "Sample name doesn't match just .csv file "+...
        "in sample_directory")

    % Read table
    df_sample = readmatrix(sample_files{n});
    
    % Select classes
    if isequal(config.use_classes,"true")
        ct = df_sample(:,end);
        for i = 1:length(classes)
        
        
        end
        
        % Make sure number of cell-types present equals the number of config.markers
        assert(length(unique(ct)) == length(types),"Number of specified cell types (%d) "+...
            "does not match unique cell type indexes (%d) in sample file.\n",length(unique(ct)),...
            length(types))
    else
        ct = ones(size(df_sample,1),1);
    end
    
    % Create count matrix
    counts = zeros(height(df_results), length(types));

    % Apply thresholds and count cells in each structure for each channel
    % First channel cell intensities start at column 5
    for i = 1:length(indexes)
        df_sub = df_sample(df_sample(:,4) == indexes(i),end);
        if isequal(config.use_classes,"true")
            for j = 1:length(classes)
                counts(i,j) = sum(df_sub == classes(j));
            end
        else
            counts(i) = length(df_sub);
        end
    end

    % Sum according to structure level order, except for background
    ids = df_results.id;
    path = df_results.structure_id_path;
    for i = 2:length(ids)
        idx = cellfun(@(s) contains(s,string("/"+ids(i)+"/")), path);
        counts(i,:) = sum(counts(idx,:),1);
    end

    % Create header name
    df_header = config.samples(s_idx) + "_" + config.groups{s_idx}(2) + "_" + types + "_Counts";

    % Convert to table
    counts_table = array2table(counts,'VariableNames',df_header);

    % Concatenate counts to results table
    df_results = horzcat(df_results, counts_table);
    
    % Print cells counted
    for j = 1:length(types)
        fprintf(strcat(char(datetime('now')),"\t Total %s Count: %d\n"),...
            types(j),counts(2,j));
    end
end

% Write new results file
writetable(df_results,count_results)

end
