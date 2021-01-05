function [config, path_table] = NM_config(stage, sample, run)
%--------------------------------------------------------------------------
% NM_config Save configured parameters for running a pipeline.
%
% Syntax:  config = NM_config(stage, sample, run)
% 
% Inputs:
%   stage - 'process', 'analyze', 'evaluate'
%   sample - (string) sample or group identifier(s) in NM_samples
%   run - (optional, logical) whether to run respective pipeline
%
% Output:
%   config - parameter configuration structure
%--------------------------------------------------------------------------

if nargin<3
    run = false;
end
    
% Elastix paths
% Uncomment and specify your path to elastix bin and library
 elastix_path_bin = '/Users/Oleh/Programs/elastix-5.0.0-mac/bin';
 elastix_path_lib = '/Users/Oleh/Programs/elastix-5.0.0-mac/lib'; 

% Make sure path is set
addpath(genpath('.'))

% Load variables
switch stage
    case {'process','stitch','align','intensity'}
        NMp_template
        main_stage = 'process';
    case 'analyze'
        NMa_template
        main_stage = 'analyze';
    case 'evaluate'
        NMe_template
        main_stage = 'evaluate';
    otherwise
        error("Invalid stage input")
end

% Save config structure
home_path = fileparts(which('NM_config.m'));
cd(home_path)
save(fullfile('templates', 'NM_variables.mat'),'-mat')

% Load and append sample info
if nargin > 1 && ~isequal(main_stage,'evaluate')
    [img_directory, output_directory] = NM_samples(sample, true);
elseif nargin > 1 && isequal(main_stage,'evaluate')
    fid = fopen('./templates/NM_samples.m');
    c = textscan(fid,'%s');
    sample_idx = c{:}(find(cellfun(@(s) isequal(s,'case'),c{:}))+1);
    for i = 1:length(sample_idx)
        [~, centroids_directory(i), group(i)] = NM_samples(sample_idx{i}(2:end-1),false);
    end
    fclose(fid);
    output_directory = results_directory;
    use_processed_images = "false";
    save('./templates/NM_variables.mat','centroids_directory','group','output_directory','-mat','-append')
else
    error("Sample information is unspecified. Set 'sample' variable.")
end

% Update image directory if using processed or analyzed images
if ~isequal(use_processed_images,"false")
    img_directory = fullfile(output_directory,use_processed_images);
    if ~exist(img_directory,'dir')
        error("Could not locate processed image directory %s\n",img_directory)
    else
        save(fullfile('templates','NM_variables.mat'),'img_directory','-mat','-append')
    end
end

% Check variable lengths for some variables
check_variable_lengths(main_stage)

% Add elastix paths if present
if exist('elastix_path_bin','var')
    path1 = getenv('PATH');
    if ~contains(path1,'elastix')
        path1 = [path1, ':', elastix_path_bin];
        setenv('PATH',path1)
    end
end
if exist('elastix_path_lib','var')
    if ismac
        path2 = getenv('DYLD_LIBRARY_PATH');
        if ~contains(path2,'elastix')
            path2 = [path2, ':', elastix_path_lib];
            setenv('DYLD_LIBRARY_PATH',path2)
        end
    else
        path2 = getenv('LD_LIBRARY_PATH');
        if ~contains(path2,'elastix')
            path2 = [path2, ':', elastix_path_lib];
            setenv('LD_LIBRARY_PATH',path2)
        end
    end
end

% Check if MATLAB toolboxes exist
try
    gcp('nocreate');
catch
    warning("Could not load Parallel Computing Toolbox. It's recommended "+...
        "that this toolbox be installed to speed up analysis.")
end

% Try setting up vl_feat if not permanently installed
if ~(exist('vl_version','file') == 3)
    try 
        vl_setup
    catch
    end
end

% Reload config
if nargout > 0
    config = load(fullfile('templates','NM_variables.mat'));
    config = orderfields(config);
end

% Run
if run
    % Make an output directory
    if exist(output_directory,'dir') ~= 7
        mkdir(output_directory);
    end

    % Make a variables directory
    var_directory = fullfile(output_directory,'variables');
    if exist(var_directory,'dir') ~= 7
        mkdir(var_directory);
    end
    
    % Copy variables to output destination
    copyfile(fullfile('templates', 'NM_variables.mat'),var_directory)
    
    % Run pipeline
    switch stage
        case 'process'
            NM_process(var_directory)
        case 'analyze'
            NM_analyze(var_directory)
        case 'evaluate'
            NM_evaluate(var_directory)
        case 'stitch'
            NM_process(var_directory,'stitch',true)
        case 'align'
           NM_process(var_directory,'align',true)
        case 'intensity'
           NM_process(var_directory,'intensity',true)
    end
elseif nargout >= 1
    % Update img_directory if using processed images
    if ~isequal(config.use_processed_images,"false")
        config.img_directory = fullfile(config.output_directory,config.use_processed_images);
    end
end

if nargout == 2
    path_table = path_to_table(config);
end

end


function check_variable_lengths(stage)
% Check user input variable lengths to make sure they're the correct length
% and match number of markers in most cases

switch stage
    case 'process'
        % Variables to check
        variable_names = {'markers','single_sheet','ls_width','laser_y_displacement','blending_method',...
            'elastix_params','rescale_intensities','subtract_background','Gamma','smooth_img','smooth_sigma',...
            'DoG_img','DoG_minmax','DoG_factor','darkfield_intensity', 'resolution','z_initial',...
            'adjust_tile_shading','adjust_tile_position','lowerThresh','upperThresh','signalThresh'};
        load(fullfile('templates','NM_variables.mat'),variable_names{:});
        
        if exist('markers') ~= 1  || isempty(markers)
            error("Must provide marker names for channel");end
        if exist('single_sheet','var') == 1 && length(single_sheet) == 1
            single_sheet = repmat(single_sheet,1,length(markers));end
        if exist('ls_width','var') == 1 && length(ls_width) == 1
            ls_width = repmat(ls_width,1,length(markers));end
        if exist('laser_y_displacement','var') == 1 && length(laser_y_displacement) == 1
            laser_y_displacement = repmat(laser_y_displacement,1,length(markers));end
        if exist('blending_method','var') == 1 && length(blending_method) == 1
            blending_method = repmat(blending_method,1,length(markers));end
        if exist('elastix_params','var') == 1 && length(elastix_params) == 1
            elastix_params = repmat(elastix_params,1,length(markers)-1);end
        if exist('rescale_intensities','var') == 1 && length(rescale_intensities) == 1
            rescale_intensities = repmat(rescale_intensities,1,length(markers));end
        if exist('subtract_background','var') == 1 && length(subtract_background) == 1
            subtract_background = repmat(subtract_background,1,length(markers));end
        if exist('Gamma','var') == 1 && length(Gamma) == 1 || isempty(Gamma)
            Gamma = ones(1,length(markers));end
        if exist('smooth_img','var') == 1 && length(smooth_img) == 1
            smooth_img = repmat(smooth_img,1,length(markers));end
        if exist('smooth_sigma','var') == 1 && isempty(smooth_sigma)
            smooth_sigma = ones(1,length(markers));end
        if exist('smooth_sigma','var') == 1 && length(smooth_sigma) == 1
            smooth_sigma = repmat(smooth_sigma,1,length(markers));end
        if exist('DoG_img','var') == 1 && length(DoG_img) == 1
            DoG_img = repmat(DoG_img,1,length(markers));end
        if exist('DoG_minmax','var') == 1 && isempty(DoG_minmax) == 1
            DoG_minmax = [0.8,1.5];end
        if exist('DoG_factor','var') == 1 && length(DoG_factor) == 1
            DoG_factor = repmat(DoG_factor,1,length(markers));end
        if exist('darkfield_intensity','var') == 1 && length(darkfield_intensity) == 1
            darkfield_intensity = repmat(darkfield_intensity,1,length(markers));end
        if exist('resolution','var') == 1
            if ~iscell(resolution)
                resolution = {resolution};
            end
            resolution = repmat(resolution,1,length(markers));
        end
        if exist('z_initial','var') == 1 
            if isempty(z_initial) 
                z_intial = zeros(1,length(markers));
            elseif length(z_initial) == 1
                z_initial = [0,repmat(z_initial,1,length(markers)-1)];
            end
        end
        if exist('adjust_tile_shading','var') == 1 && length(adjust_tile_shading) == 1
            adjust_tile_shading = repmat(adjust_tile_shading,1,length(markers));end
        if exist('adjust_tile_position','var') == 1 && length(adjust_tile_position) == 1
            adjust_tile_position = repmat(adjust_tile_position,1,length(markers));end
        if exist('lowerThresh','var') == 1 && ~isempty(lowerThresh)
            assert(length(lowerThresh) == length(markers),"Lower threshold values need to be "+...
                "specified for all markers or left empty")
            if isempty(lowerThresh) && ~isempty(upperThresh)
                upperThresh = zeros(0,1,length(markers));
            end
        end
        if exist('upperThresh','var') == 1
            if ~isempty(upperThresh)
                assert(length(upperThresh) == length(markers),"Upper threshold values need to be "+...
                    "specified for all markers or left empty")
            end
            if isempty(upperThresh) && ~isempty(lowerThresh)
                upperThresh = repmat(65535,1,length(markers));
            end
        end
        if exist('signalThresh','var') == 1
            if ~isempty(signalThresh)
                assert(length(signalThresh) == length(markers),"Signal threshold values need to be "+...
                    "specified for all markers or left empty")
            elseif ~isempty(lowerThresh) || ~isempty(upperThresh)
                error("signalThresh must specified if lowerThresh and/or upperThresh "+...
                    "is also specified")
            end
        end
    case 'analyze'
        % Variables to check
        variable_names = {};
        load(fullfile('templates','NM_variables.mat'),variable_names{:});
end

% Resave these variables
save(fullfile('templates', 'NM_variables.mat'),variable_names{:},'-mat','-append')

end