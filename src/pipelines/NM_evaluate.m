function NM_evaluate(config, step, varargin)
%--------------------------------------------------------------------------
% NuMorph evaluation pipeline to merge analysis results from multiple
% grouped samples and perform pairwise statistical comparisons between WT
% and KO groups.
%--------------------------------------------------------------------------
config.home_path = fileparts(which('NM_config'));

% Generate evaluate config if provided analysis config for one sample
if isequal(config.main_stage, 'analyze')
    config = NM_config('evaluate', config.sample_id);
end

% Default to run full pipeline
if nargin<2
    step = 'stats';
end

% Make an output directory
if ~isfolder(config.results_directory)
    mkdir(config.results_directory);
end

% Check for annotation table template
if isequal(config.compare_structures_by, "table")


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
