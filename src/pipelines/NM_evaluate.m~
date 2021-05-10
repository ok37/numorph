function NM_evaluate(config, step, plot_type)
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

if nargin<3
    plot_type = 'coronal';
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

% Check for nuclear channel
%if isequal(config.contains_nuclear_channel,"true")
%    config.class_names = ["nuclei",config.class_names];
%end

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
    config.stats_results = fullfile(config.results_directory,sprintf('%s_%s_stats.csv',config.prefix,a));
else
    config.stats_results = fullfile(config.results_directory,strcat(config.prefix,'_summary_stats.csv'));
end

%% Run single step and return if specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(step,'stats') 
    perform_stats(config,multi);
    fprintf('%s\t All statistics generated! \n',datetime('now'))
    return
elseif isequal(step,'plot')
    visualize_results(config,plot_type);
    return
end

end


function perform_stats(config,multi)

n_samples = length(config.samples);
% Measure cell-types
count_results = fullfile(config.results_directory,strcat(config.prefix,"_summary_counts.csv"));
if ~isfile(count_results)
    fprintf('%s\t Quantifying cell counts for %d samples\n',datetime('now'),n_samples)    
    combine_counts(config,count_results);
else
    fprintf('%s\t Cell count results already quantified \n',datetime('now'))        
end

% Measure volumes
vol_results = fullfile(config.results_directory,strcat(config.prefix,"_summary_volumes.csv"));
if ~isfile(vol_results)
    fprintf('%s\t Quantifying structure volumes for %d samples\n',datetime('now'),n_samples)
    combine_volumes(config,vol_results);
else
    fprintf('%s\t Structure volumes already quantified \n',datetime('now'))        
end

% Measure densities
dense_results = fullfile(config.results_directory,strcat(config.prefix,"_summary_densities.csv"));
if ~isfile(dense_results)
    fprintf('%s\t Quantifying cell densities for %d samples\n',datetime('now'),n_samples)
    combine_densities(config,count_results,vol_results,dense_results);
else
    fprintf('%s\t Cell densities already quantified \n',datetime('now'))        
end

% Statistics run on mutliple samples
if ~multi || isempty(config.groups)
    fprintf('%s\t Only 1 sample so skipping statistics...\n',datetime('now'))
    return
end

% Read counts and volume data directly from csv file
df_results = cell(1,2);
fname = fullfile(config.results_directory,strcat(config.prefix,"_summary_counts.csv"));
if isfile(fname)
    fprintf('%s\t Loading cell counts \n',datetime('now'))    
    df_results{1} = readtable(fname,'PreserveVariableNames',true);
else
    error('Could not locate %s in results directory',fname)
end
fname = fullfile(config.results_directory,strcat(config.prefix,"_summary_volumes.csv"));
if isfile(fname)
    fprintf('%s\t Loading structure volumes \n',datetime('now'))    
    df_results{2} = readtable(fname);
else
    error('Could not locate %s in results directory',fname)
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

% Calculate statistics
config.groups = config.groups';
config.samples = config.samples';
%calc_stats(df_results, config);

% Calculate voxel image
calc_voxel_stats(config)

end


function visualize_results(config,plot_type)
% Visualize Results

if ismember(plot_type,{'coronal','saggital','axial','voxel'})
    main_plot = 'voxel';
    if isequal(string(plot_type),"voxel")
        plot_type = 'coronal';
    end
else
    main_plot = plot_type;
end

% Check for previous structure
plot_path = fullfile(config.results_directory,sprintf("%s_plotData.mat",config.prefix));
if ~isfile(plot_path)
    plotData = struct;
    save(plot_path,'-struct','plotData')
    var_names = {};
else
    var_names = who('-file',plot_path);
end

% Check if pre-loading existing munged data
new_data = false;
if isequal(config.update_plotdata,"true")
    new_data = true;
end

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
                df_stats = readtable(config.stats_results,'PreserveVariableNames',true);
            else
                error("Could not locate stats table %s",stats_results)
            end
            
            % Read voxel images
            for i = 1:length(config.class_names)
                voxel_imgs(i).marker = config.class_names(i);
                voxel_imgs(i).data =  niftiread(fullfile(config.results_directory,...
                    config.prefix + "_" + config.class_names(i) + "_voxels.nii"));
            end

            % Create visualization volume
            vox_results.(plot_type) = create_stats_volume2(df_stats,voxel_imgs,...
                config.class_names,string(plot_type));
            save(plot_path,'-append','-struct','vox_results',plot_type)
            vox = vox_results.(plot_type);
        end

        % Display Slice
        % z: z_position
        % marker: annotation, c1,c2,c3...
        % category: volume, counts, density, voxel
        % stat: Mean WT, Mean KO, StDev WT, StDev KO, Fold Change, p, p.adj, sig
        % [30,50,70,90]
        key{1} = [90,4,4,5];
        key{2} = [90,4,4,7];

        close(gcf)
        [p1, p2] = display_slice2(vox,key);
    case 'cloud'
        
    case 'flatmap'
        % Create plot data structure and save for later
        if ~new_data && ismember(plot_type,var_names)
            fprintf('%s\t Loading visualization volume\n',datetime('now'))
            load(plot_path,'flatmap')
        else        
            % Read saved flatmap images
            for i = 1:length(config.class_names)
                flatmap(i).marker = config.class_names(i);
                flatmap(i).data =  niftiread(fullfile(config.results_directory,...
                    config.prefix + "_" + config.class_names(i) + "_flatmap.nii"));
            end
            save(plot_path,'-append','flatmap')
        end
        
        % Visualize flat cortex
        key = [1,1];
        visualize_flat_cortex(flatmap, key)

    case 'volume'
        
    otherwise
        error("Unrecognized plot type")        
end


end

