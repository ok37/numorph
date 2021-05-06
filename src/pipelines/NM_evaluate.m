function NM_evaluate(config, step, plot_type)
%--------------------------------------------------------------------------
% NuMorph evaluation pipeline to merge analysis results from multiple
% grouped samples and perform pairwise statistical comparisons between WT
% and KO groups.
%--------------------------------------------------------------------------
config.overwrite = "true";
config.temp_file = fullfile(fileparts('NM_config'),'annotations','structure_template.csv');

% Default to run full pipeline
if nargin<2
    step = 'stats';
end

% Make an output directory
if ~isfolder(config.results_directory)
    mkdir(config.results_directory);
end

% Check class names
if isempty(config.class_names) && isequal(config.use_classes,"true")
    config.class_names = "CT" + config.keep_classes;
end

%% Run single step and return if specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(step,'stats') 
    perform_stats(config);
    return
elseif isequal(step,'visualize')
    perform_stats(config);
    return
end

%% Count Cell Types
count_results = fullfile(config.results_directory,"NMe_summary_counts.csv");
if isequal(config.combine_counts,"update") ||...
        (isequal(config.combine_counts,"true") && ~isfile(count_results))
    
    n_samples = length(config.samples);
    fprintf('%s\t Quantifying cell counts for %d samples\n',datetime('now'),n_samples)
    
    % Get sample files
    sample_files = repmat("",1,n_samples);
    for i = 1:n_samples
        f = dir(config.centroids_directory(i));
        if isequal(config.use_classes,"true")
            f = f(arrayfun(@(s) endsWith(s.name, '.csv') && contains(s.name,strcat(config.samples(i),'_classes')),f));
        else
            f = f(arrayfun(@(s) endsWith(s.name, '.csv') && contains(s.name,strcat(config.samples(i),'_centroids')),f));
        end
        if length(f)>1
            f = f(end);
        end
        sample_files(i) = string(fullfile(config.centroids_directory(i),f.name));
    end

    % Count cells and merge
    combine_counts(sample_files,config);
elseif isequal(config.combine_counts,"true") && isfile(count_results)
        fprintf('%s\t Cell counts already quantified \n',datetime('now'))
end

%% Measure Volumes
vol_results = fullfile(config.results_directory,"NMe_summary_volumes.csv");
if isequal(config.combine_volumes,"update") ||...
        (isequal(config.combine_volumes,"true") && ~isfile(vol_results))
    
    n_samples = length(config.samples);
    fprintf('%s\t Quantifying structure volumes for %d samples \n',datetime('now'),n_samples)
    
    % Get sample files
    sample_files = repmat("",1,n_samples);
    run = true;
    for i = 1:n_samples
        sample_files(i) = string(fullfile(config.centroids_directory(i),...
            strcat(config.samples(i),'_volumes.csv')));
        if ~isfile(sample_files(i))
            fprintf("Could not locate volume data for sample %s. Skipping "+...
                "volume calculations \n",string(config.samples{i}));
            run = false;
        end
    end
    
    % Measure volumes and merge
    if run
        combine_volumes(sample_files,config);
    end
elseif isequal(config.combine_counts,"true") && isfile(count_results)
    fprintf('%s\t Volumes already quantified \n',datetime('now'))
end

%% Run Statistics
% Check for groups
if isempty(config.groups)
    fprintf('%s\t Only 1 sample so skipping statistics...\n',datetime('now'))
    fprintf('%s\t Evaluation completed! \n',datetime('now'))
    return
end

if isequal(config.compare_structures_by,'csv')
    [~,a] = fileparts(config.structure_csv);
    stats_results = fullfile(config.results_directory,strcat('NMe_',a,'_stats.csv'));
else
    stats_results = fullfile(config.results_directory,strcat('NMe_summary_stats.csv'));
end

if isequal(config.calculate_stats,'update') ||...
        (isequal(config.calculate_stats,'true') && ~isfile(stats_results))
    df_results = cell(1,2);
    if isequal(config.use_data,'counts') || isequal(config.use_data,'both')
        if ~exist('df_counts','var')
        % If results table if not defined, attempt to load from results
        % path
        try
            fprintf(strcat(char(datetime('now')),"\t Loading cell counts\n"))
            df_counts = readtable(fullfile(config.results_directory,'NMe_summary_counts.csv'));
        catch ME
            error('Could not locate NMe_summary_counts.csv')
        end
        df_results{1} = df_counts;
        end
        
        % Get sample info from 
        c =  cellfun(@(s) strsplit(s,'_'),df_results{1}.Properties.VariableNames(10:end),'UniformOutput',false);        
        c = cat(1,c{:});
        
        % Edit this, samples and groups should be read directly from column
        % names
        idx = ismember(config.samples,string(unique(c(:,1))))';
        config.samples = config.samples(idx);
        config.groups = config.groups(idx);
        config.markers = string(unique(c(:,2),'stable'));
    end
    if isequal(config.use_data,'volumes') || isequal(config.use_data,'both')
        if ~exist('df_volumes','var')
        % If results table if not defined, attempt to load from results
        % path
        try
            fprintf(strcat(char(datetime('now')),"\t Loading volume measurements\n"))
            df_volumes = readtable(fullfile(config.results_directory,'NMe_summary_volumes.csv'));
        catch ME
            error('Could not locate NMe_summary_volumes.csv') 
        end
        df_results{2} = df_volumes;
        end
    end
    
     % Check that number of each group is the same for each sample
    try
        groups = cat(1,config.groups{:});
    catch
       error("Each sample should have the same number of group assigments") 
    end
    
    % Get sub-group
    n_groups = size(groups,2);
    if n_groups>1
        sub_groups = unique(groups(:,2));
        for i = 1:length(sub_groups)
            fprintf('%s\t Found %d samples in sub-group %s\n',datetime('now'),...
                sum(groups(:,2) == sub_groups(i)),sub_groups(i))
        end
    end

    % Get additional categorical covariates
    if n_groups>2
        for i = 3:n_groups
            fprintf('%s\t Found %d categorical group with %d unique values \n',datetime('now'),...
                i-2, length(unique(groups(:,i))))
        end
    end
    % Calculate statistics
    df_stats = calc_stats(df_results, config);
else
    fprintf('%s\t Loading calculated statistics \n',datetime('now'))
end

fprintf('%s\t Evaluation completed! \n',datetime('now'))
return

%% Visualize Results
if isequal(config.orientation,"coronal")
    vol_path = fullfile(fileparts(which('NM_config')),'data','coronalSliceVol.mat');
elseif isequal(config.orientation,"sagittal")
    vol_path = fullfile(fileparts(which('NM_config')),'data','sagittalSliceVol.mat');
elseif isequal(config.orientation,"axial")
    vol_path = fullfile(fileparts(which('NM_config')),'data','axialSliceVol.mat');
end

if isequal(config.visualize_results,'update') ||...
        (isequal(config.visualize_results,'true') && ~isfile(vol_path))
    
    % Read stats table
    if isfile(stats_results)
        df_stats = readtable(stats_results,'PreserveVariableNames',true);
    else
        error("Could not locate stats table %s",stats_results)
    end

    % Create visualization volume
    fprintf('%s\t Generating visualization volume\n',datetime('now'))
    markers = config.class_names(config.keep_classes);
    vs = create_stats_volume2(df_stats,markers,config.orientation,config.custom_comp);
else
    fprintf('%s\t Loading visualization volume\n',datetime('now'))
    load(vol_path,'vs')
end

%% Display Slice
% z: z_position
% marker: annotation, c1,c2,c3...
% category: volume, counts, density, custom
% stat: Mean WT, Mean Top1, StDev WT, StDev Top1, Fold Change, p, p.adj,
% sig
key{1} = [56,4,4,5];
key{2} = [56,4,4,7];

close(gcf)
[p1, p2] = display_slice2(vs,key);
end


function perform_stats(config)

% Check if 1 or more samples
samples = config.samples;
n_samples = length(config.samples);
if n_samples > 1
    multi = true;
    config.prefix = config.groups{1}(1);
else
    multi = false;
    config.prefix = config.samples(1);
end

% Measure cell-types
fprintf('%s\t Quantifying cell counts for %d samples\n',datetime('now'),n_samples)    
%combine_counts(config);

% Measure volumes and densities
fprintf('%s\t Quantifying structure volumes for %d samples\n',datetime('now'),n_samples)    
combine_volumes(config);


% Run statistics if more than 1 sample
if multi
    % Check for groups
    if isempty(config.groups)
        fprintf('%s\t Only 1 sample so skipping statistics...\n',datetime('now'))
        fprintf('%s\t Evaluation completed! \n',datetime('now'))
        
        
        %%% Print tables
        return
    end

    if isequal(config.compare_structures_by,'csv')
        [~,a] = fileparts(config.structure_csv);
        stats_results = fullfile(config.results_directory,strcat('NMe_',a,'_stats.csv'));
    else
        stats_results = fullfile(config.results_directory,strcat('NMe_summary_stats.csv'));
    end

    if isequal(config.calculate_stats,'update') ||...
            (isequal(config.calculate_stats,'true') && ~isfile(stats_results))
        df_results = cell(1,2);
        if isequal(config.use_data,'counts') || isequal(config.use_data,'both')
            if ~exist('df_counts','var')
            % If results table if not defined, attempt to load from results
            % path
            try
                fprintf(strcat(char(datetime('now')),"\t Loading cell counts\n"))
                df_counts = readtable(fullfile(config.results_directory,'NMe_summary_counts.csv'));
            catch ME
                error('Could not locate NMe_summary_counts.csv')
            end
            df_results{1} = df_counts;
            end

            % Get sample info from 
            c =  cellfun(@(s) strsplit(s,'_'),df_results{1}.Properties.VariableNames(10:end),'UniformOutput',false);        
            c = cat(1,c{:});

            % Edit this, samples and groups should be read directly from column
            % names
            idx = ismember(config.samples,string(unique(c(:,1))))';
            config.samples = config.samples(idx);
            config.groups = config.groups(idx);
            config.markers = string(unique(c(:,2),'stable'));
        end
        if isequal(config.use_data,'volumes') || isequal(config.use_data,'both')
            if ~exist('df_volumes','var')
            % If results table if not defined, attempt to load from results
            % path
            try
                fprintf(strcat(char(datetime('now')),"\t Loading volume measurements\n"))
                df_volumes = readtable(fullfile(config.results_directory,'NMe_summary_volumes.csv'));
            catch ME
                error('Could not locate NMe_summary_volumes.csv') 
            end
            df_results{2} = df_volumes;
            end
        end

         % Check that number of each group is the same for each sample
        try
            groups = cat(1,config.groups{:});
        catch
           error("Each sample should have the same number of group assigments") 
        end

        % Get sub-group
        n_groups = size(groups,2);
        if n_groups>1
            sub_groups = unique(groups(:,2));
            for i = 1:length(sub_groups)
                fprintf('%s\t Found %d samples in sub-group %s\n',datetime('now'),...
                    sum(groups(:,2) == sub_groups(i)),sub_groups(i))
            end
        end

        % Get additional categorical covariates
        if n_groups>2
            for i = 3:n_groups
                fprintf('%s\t Found %d categorical group with %d unique values \n',datetime('now'),...
                    i-2, length(unique(groups(:,i))))
            end
        end
        % Calculate statistics
        df_stats = calc_stats(df_results, config);
    else
        fprintf('%s\t Loading calculated statistics \n',datetime('now'))
    end
end

fprintf('%s\t Evaluation completed! \n',datetime('now'))


end

