function NM_evaluate(config, step, varargin)
%--------------------------------------------------------------------------
% NuMorph evaluation pipeline to merge analysis results from multiple
% grouped samples and perform pairwise statistical comparisons between WT
% and KO groups.
%--------------------------------------------------------------------------

% Generate evaluate config if provided analysis config for one sample
if isequal(config.main_stage, 'analyze')
    config2 = NM_config('evaluate', config.sample_id);
    
    % Specify classes to keep and default class names
    if isequal(config2.use_classes,"true")
        has_classes = load_results(config,'classes',[],true);
        if has_classes
            config2.keep_classes = config.keep_classes;
            classes = unique(load_results(config,'classes'));
            config2.class_names = ["CT" + num2str(classes)]';
            if ~isempty(config2.keep_classes)
                config2.class_names = config2.class_names(config2.keep_classes);
            else
                config2.keep_classes = classes;
            end
        else
            error('Results structure does not contain cell classes')
        end
    end
    config = config2;
end
config.home_path = fileparts(which('NM_config'));

%%%%% Some defaults here
config.measure_cortex = "false"; % Set this to false for the time being
config.sum_all_classes = "true"; % We should always want to count total nuclei
config.paired = "false";         % This worked previously but buggy now. Easier to calculate state outside of matlab

% Default to run full pipeline
if nargin<2
    step = 'stats';
end

% Make an output directory
if ~isfolder(config.results_directory)
    mkdir(config.results_directory);
end

% Check for annotation table template
if isequal(config.compare_structures_by, "table") && ... 
    ~isfile(fullfile(config.home_path, 'annotations', config.template_file))
    fprintf('%s\t Structure template %s was not found in /annotations. Comparing by index. \n',datetime('now'))
    config.compare_structures_by = "index";
end

% Check class names
if isequal(config.use_classes,"true") && isempty(config.class_names)
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
    [~,b] = fileparts(config.structure_table);
else
    [~,b] = fileparts(config.template_file);
end
b = erase(b,"_structure_template");
config.stats_results = fullfile(config.results_directory,config.prefix,'tables',sprintf('%s_%s_summary_stats.xls',config.prefix,b));

%% Run single step and return if specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(step,'stats') 
    perform_stats(config,multi);
    fprintf('%s\t All statistics generated! \n',datetime('now'))
    return
elseif isequal(step,'plot')
    % Using apps as main method for plotting
    visualize_results(config,varargin);
    return
end

end
