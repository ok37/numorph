function NM_evaluate(config, step, varargin)
%--------------------------------------------------------------------------
% NuMorph evaluation pipeline to merge analysis results from multiple
% grouped samples and perform pairwise statistical comparisons between WT
% and KO groups.
%--------------------------------------------------------------------------
config.overwrite = "true";
config.home_path = fileparts(which('NM_config'));
config.temp_file = fullfile(fileparts(which('NM_config')),'annotations','structure_template.csv');

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
    if ~isempty(config.custom_class)
        config.class_names = [config.class_names,config.custom_class];
    end
else
    assert(length(config.class_names) == length(config.keep_classes)+length(config.custom_class),...
        "Number of class names does not match number of classes selected")
end

% Check if 1 or more samples
n_samples = length(config.samples);
if n_samples > 1
    multi = true;
    config.prefix = config.groups{1}(1);
else
    multi = false;
    config.prefix = config.samples(1);
end

% Name summary file
if isequal(config.compare_structures_by,'table')
    [~,a] = fileparts(config.structure_table);
    config.stats_results = fullfile(config.results_directory,config.prefix,'tables',sprintf('%s_%s_stats.xls',config.prefix,a));
else
    config.stats_results = fullfile(config.results_directory,config.prefix,'tables',strcat(config.prefix,'_summary_stats.xls'));
end

%% Run single step and return if specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(step,'stats') 
    perform_stats(config,multi);
    fprintf('%s\t All statistics generated! \n',datetime('now'))
    return
elseif isequal(step,'plot')
    % Using apps as main method for plotting
    return
    %visualize_results(config,varargin);
end

end


function perform_stats(config,multi)

% Start combining regions statistics
config.tab_directory = fullfile(config.results_directory,config.prefix,'tables');
if ~isfolder(config.tab_directory)
    mkdir(config.tab_directory)
end

config.vox_directory = fullfile(config.results_directory,config.prefix,'voxels');
if ~isfolder(config.vox_directory)
    mkdir(config.vox_directory)
end

config.flat_directory = fullfile(config.results_directory,config.prefix,'flatmaps');
if ~isfolder(config.flat_directory)
    mkdir(config.flat_directory)
end

n_samples = length(config.samples);

% Measure cell-types
count_results = fullfile(config.tab_directory,strcat(config.prefix,"_summary_counts.csv"));
if ~isfile(count_results)
    fprintf('%s\t Quantifying cell counts for %d samples\n',datetime('now'),n_samples)    
    combine_counts(config,count_results);
else
    fprintf('%s\t Cell count results already quantified \n',datetime('now'))
end

% Measure volumes
vol_results = fullfile(config.tab_directory,strcat(config.prefix,"_summary_volumes.csv"));
if ~isfile(vol_results)
    fprintf('%s\t Quantifying structure volumes for %d samples\n',datetime('now'),n_samples)
    combine_volumes(config,vol_results);
else
    fprintf('%s\t Structure volumes already quantified \n',datetime('now'))
end

% Measure densities
dense_results = fullfile(config.tab_directory,strcat(config.prefix,"_summary_densities.csv"));
if ~isfile(dense_results)
    fprintf('%s\t Quantifying cell densities for %d samples\n',datetime('now'),n_samples)
    combine_densities(config,count_results,vol_results,dense_results);
else
    fprintf('%s\t Cell densities already quantified \n',datetime('now'))        
end

% Measure cortex 
cortex_results = fullfile(config.tab_directory,strcat(config.prefix,"_summary_ctxVolSATH.csv"));
if isequal(config.measure_cortex,"true")
    if ~isfile(cortex_results)
        fprintf('%s\t Measuring cortical structure for %d samples\n',datetime('now'),n_samples)
        ctx_idx = arrayfun(@(s) {who('-file',s)},config.results_path);
        ctx_idx = cellfun(@(s) any(ismember(s,{'cortex'})),ctx_idx);
        if ~all(ctx_idx)
            measure_cortical_sa_th(config.results_path(~ctx_idx));
        end
        combine_cortexes(config,cortex_results);
        
    else
        fprintf('%s\t Cortical structures already quantified. Loading measurments  \n',datetime('now'))        
        
    end
end

% Statistics run on mutliple samples
if ~multi || isempty(config.groups)
    fprintf('%s\t Only 1 sample so skipping statistics...\n',datetime('now'))
    return
end

% Read counts and volume data directly from csv file
df_results = cell(1,3);
if isfile(count_results)
    fprintf('%s\t Loading cell counts \n',datetime('now'))    
    df_results{1} = readtable(count_results,'PreserveVariableNames',true);
else
    error('Could not locate %s in results directory',count_results)
end
if isfile(vol_results)
    fprintf('%s\t Loading structure volumes \n',datetime('now'))    
    df_results{2} = readtable(vol_results);
else
    error('Could not locate %s in results directory',vol_results)
end
if isequal(config.measure_cortex,"true")
    if isfile(cortex_results)
        fprintf('%s\t Loading cortex measurements \n',datetime('now'))    
        df_results{3} = readtable(cortex_results);
    else
        error('Could not locate %s in results directory',cortex_results)
    end
end

% Get sample info from 
c =  cellfun(@(s) strsplit(s,'_'),df_results{1}.Properties.VariableNames(10:end),'UniformOutput',false);        
c = cat(1,c{:});

% Edit this, samples and groups should be read directly from column
% names
idx = ismember(config.samples,string(unique(c(:,1))))';
config.samples = config.samples(idx);
config.groups = config.groups(idx);
config.markers = string(unique(c(:,3),'stable'));

% Get sub-group
groups = cat(1,config.groups{:});
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

% Calculate region statistics
config.groups = config.groups';
config.samples = config.samples';
calc_table_stats(df_results, config);

% Calculate voxel image
%calc_voxel_stats(config)

% Calculate flatmap stats
calc_flatmap_stats(config)

end


function visualize_results(config,plot_params)
% Visualize Results

plot_type = plot_params{1};
if ismember(plot_type,{'coronal','saggital','axial','voxel'})
    main_plot = 'voxel';
    if isequal(string(plot_type),"voxel")
        plot_type = 'coronal';
    end
elseif ismember(plot_type,{'thickness','flatmap'})
    main_plot = 'flatmap';
else
    main_plot = plot_type;
end

% Check for previous structure
plot_path = fullfile(config.results_directory,config.prefix,sprintf("%s_plotData.mat",config.prefix));
if ~isfile(plot_path)
    plotData = struct;
    save(plot_path,'-struct','plotData')
    var_names = {};
else
    var_names = who('-file',plot_path);
end

% Check if pre-loading existing munged data
new_data = false;

% Pick what to plot
switch main_plot
    case 'voxel'
        % Create plot data structure and save for later
        if ~new_data && ismember(plot_type,var_names)
            fprintf('%s\t Loading visualization volume\n',datetime('now'))
            vox = load(plot_path,plot_type);
            vox = vox.(plot_type);
            new_data = false;
        else
            new_data = true;
        end

        if new_data
            fprintf('%s\t Generating visualization volume\n',datetime('now'))
            % Read stats table
            if isfile(config.stats_results)
                df_stats = munge_stats(config.stats_results);
            else
                error("Could not locate stats table %s",config.stats_results)
            end
            
            % Set markers
            markers = config.class_names;
            if isequal(config.sum_all_classes,"true")
                markers = [markers,"all"];
            end
            
            % Read voxel images
            markers2 = arrayfun(@(s) strrep(s,'./',''),markers);
            for i = 1:length(markers)
                voxel_imgs(i).marker = markers(i);
                voxel_imgs(i).data =  niftiread(fullfile(config.results_directory,...
                    config.prefix,'voxels',...
                    config.prefix + "_" + markers2(i) + "_voxels.nii"));
            end
            
            % Create visualization volume
            vox_results.(plot_type) = create_stats_volume(df_stats,voxel_imgs,...
                markers,string(plot_type));
            save(plot_path,'-append','-struct','vox_results',plot_type)
            vox = vox_results.(plot_type);
        end

        % Display Slice
        % z: z_position
        % marker: annotation, c1,c2,c3...
        % category: volume, counts, density, voxel
        % stat: Mean WT, Mean KO, StDev WT, StDev KO, Fold Change, p, p.adj, sig
        % [30,50,70,90] 
        
        
        
        
    case 'cloud'
        
    case 'flatmap'
        % Create plot data structure and save for later
        %if ~new_data && ismember(plot_type,var_names)
        %    fprintf('%s\t Loading visualization volume\n',datetime('now'))
        %    load(plot_path,'flatmap')
        %else        
        %    % Read saved flatmap images
        %    class_names2 = arrayfun(@(s) strrep(s,'./',''),config.class_names);
        %    for i = 1:length(config.class_names)
        %        flatmap(i).marker = config.class_names(i);
        %        flatmap(i).data =  niftiread(fullfile(config.results_directory,...
        %            config.prefix + "_" + class_names2(i) + "_flatmap.nii"));
        %    end
        %    save(plot_path,'-append','flatmap')
        %end

        % Visualize flat cortex
        key = [1,1];
        plot_flat_cortex(flatmap, key)

    case 'volume'
        
    otherwise
        error("Unrecognized plot type")        
end


end

