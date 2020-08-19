clear
% This script merges analysis results from multiple samples to performs
% pairwise comparison between WT and KO groups. Outputs a merged .csv files
% containing cell counts for each structure. Additional statistical
% comparisons are also performed. 
id = sprintf('TCe_%s.out',char(datetime('now','Format','yyyy-MM-dd-hhmmss')));
config = load('TCe_variables.mat');

%% Count Cell Types
if isequal(config.combine_counts,'true')
    fprintf(strcat(char(datetime('now')),"\t Quantifying cell counts\n"))
    f = dir(config.sample_directory);
    f = f(arrayfun(@(s) contains(s.name, '.csv') && contains(s.name,'centroids'),f));
    sample_files = arrayfun(@(s) string(fullfile(s.folder, s.name)),f);

    % Count cells and merge
    df_counts = combine_counts(sample_files,config);
end

%% Measure Volumes
if isequal(config.combine_volumes,'true')
    fprintf(strcat(char(datetime('now')),"\t Quantifying cell volumes\n"))
    f = dir(config.sample_directory);
    f = f(arrayfun(@(s) contains(s.name, '.csv') && contains(s.name,'volumes'),f));
    if isempty(f)
        fprintf(strcat(char(datetime('now')),"\t No structure volume data found. "+...
            "Searching for mask files to measure structure volumes.\n"))
        f = dir(config.sample_directory);
        f = f(arrayfun(@(s) contains(s.name, '.mat') && contains(s.name,'mask'),f));
        assert(~isempty(f), 'No mask volumes found')
        sample_files = arrayfun(@(s) string(fullfile(s.folder, s.name)),f);
        df_volumes = measure_structure_volumes(sample_files, config.results_directory,...
            config.sum_child_structures);
    else 
        assert(~isempty(f), 'No volume .csv files found')
        sample_files = arrayfun(@(s) string(fullfile(s.folder, s.name)),f);
    
        % Measure volumes and merge
        df_volumes = combine_volumes(sample_files,config.results_directory,...
            config.thresholds,config.sample_names,config.markers,config.sum_child_structures,config.overwrite);
    end
end

%% Run Statistics
if isequal(config.calculate_stats,'true')
    df_results = cell(1,2);
    if isequal(config.use_data,'counts') || isequal(config.use_data,'both')
        if ~exist('df_counts','var')
        % If results table if not defined, attempt to load from results
        % path
        try
            fprintf(strcat(char(datetime('now')),"\t Loading cell counts\n"))
            df_counts = readtable(fullfile(config.results_directory,'TCe_summary_counts.csv'));
        catch ME
            error('Could not locate TCe_summary_counts.csv')
        end
        df_results{1} = df_counts;
        end
    end
    if isequal(config.use_data,'volumes') || isequal(config.use_data,'both')
        if ~exist('df_volumes','var')
        % If results table if not defined, attempt to load from results
        % path
        try
            fprintf(strcat(char(datetime('now')),"\t Loading volume measurements\n"))
            df_volumes = readtable(fullfile(config.results_directory,'TCe_summary_volumes.csv'));
        catch ME
            error('Could not locate TCe_summary_volumes.csv') 
        end
        df_results{2} = df_volumes;
        end
    end

    % Calculate statistics
    df_stats = calc_stats(df_results, config);
end

%% Visualize Results
if isequal(config.visualize_results,'true')
    if ~exist('df_stats','var')
        % If results table if not defined, attempt to load from results
        % path
        try
            df_stats = readtable(fullfile(config.results_directory,'TCe_summary_stats.csv'));
        catch ME
            error('Could not locate TCe_summary_stats.csv') 
        end
    end
    % Create visualization volume
    fprintf(strcat(char(datetime('now')),"\t Loading visualization volume \n"))
    config.compare_structures_by = 'csv';
    vs = create_stats_volume2(df_stats,config.markers,config.structure_depth,...
        config.compare_structures_by, config.orientation);
end


%% Display Slice
% z: z_position
% marker: annotation, ToPro, Ctip2, Cux1
% category: volume, counts, density
% stat: Mean WT, Mean Top1, StDev WT, StDev Top1, Fold Change, p, p.adj,
% sig
key{1} = [30,1,2,5];
key{2} = [30,1,2,7];

%close(gcf)
[p1, p2] = display_slice2(vs,key);

%% Display flat cortex







