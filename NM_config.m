function config = NM_config(stage, sample, run)
%--------------------------------------------------------------------------
% NM_config Save configured parameters for running a pipeline.
%
% Syntax:  config = NM_config(stage, sample, run)
% 
% Inputs:
%   stage - 'process', 'analyze', 'evaluate'
%   sample - (string) sample or group identifier(s) in NMsamples
%   run - (optional, logical) whether to run respective pipeline
%
% Output:
%   config - parameter configuration structure
%--------------------------------------------------------------------------

% Elastix paths
% Uncomment and specify your path to elastix bin and library
 elastix_path_bin = '/Users/Oleh/Programs/elastix-5.0.0-mac/bin';
 elastix_path_lib = '/Users/Oleh/Programs/elastix-5.0.0-mac/lib';

% Make sure path is set
addpath(genpath('.'))

% Load variables
switch stage
    case 'process'
        NMp_template
    case 'analyze'
        NMa_template
    case 'evaluate'
        NMe_template
    otherwise
        error("Invalid stage input")
end

% Save config structure
home_path = fileparts(which('NM_config.m'));
cd(home_path)
save(fullfile('templates', 'NM_variables.mat'),'-mat')

% Load and append sample info
if nargin > 1 && ~isequal(stage,'evaluate')
    [img_directory, output_directory] = NMsamples(sample, true);
elseif nargin > 1 && isequal(stage,'evaluate')
    fid = fopen('./templates/NMsamples.m');
    c = textscan(fid,'%s');
    sample_idx = c{:}(find(cellfun(@(s) isequal(s,'case'),c{:}))+1);
    for i = 1:length(sample_idx)
        [~, centroids_directory(i), group(i)] = NMsamples(sample_idx{i}(2:end-1),false);
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
check_variable_lengths(stage)

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

% Reload config
if nargout == 1
    config = load(fullfile('templates','NM_variables.mat'));
end

% Run
if nargin>2 && run
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
    end
end
end


function check_variable_lengths(stage)
% Check user input variable lengths to make sure they're the correct length
% and match number of markers in most cases

switch stage
    case 'process'
        % Variables to check
        variable_names = {'markers','single_sheet','ls_width','laser_y_displacement','blending_method',...
            'param_folder','rescale_intensities','subtract_background','gamma','smooth_img','smooth_sigma',...
            'DoG_img','DoG_minmax','DoG_factor','darkfield_intensity', 'update_intensity_channels'};
        load(fullfile('templates','NM_variables.mat'),variable_names{:});

        for i = 1:length(variable_names)
            if exist('single_sheet','var') == 1 && length(single_sheet) == 1 && i == 1
                single_sheet = repmat(single_sheet,1,length(markers));
            elseif exist('ls_width','var') == 1 && length(ls_width) == 1 && i == 2
                ls_width = repmat(ls_width,1,length(markers)); idx(3) = false;
            elseif exist('laser_y_displacement','var') == 1 && length(laser_y_displacement) == 1 && i == 3
                laser_y_displacement = repmat(laser_y_displacement,1,length(markers));
            elseif exist('blending_method','var') == 1 && length(blending_method) == 1 && i == 4
                blending_method = repmat(blending_method,1,length(markers));
            elseif exist('param_folder','var') == 1 && length(param_folder) == 1 && i == 5
                param_folder = repmat(param_folder,1,length(markers)-1);
            elseif exist('rescale_intensities','var') == 1 && length(rescale_intensities) == 1 && i == 6
                rescale_intensities = repmat(rescale_intensities,1,length(markers));
            elseif exist('subtract_background','var') == 1 && length(subtract_background) == 1 && i == 7
                subtract_background = repmat(subtract_background,1,length(markers));
            elseif exist('gamma','var') == 1 && length(gamma) == 1 || isempty(gamma) && i == 8
                gamma = ones(1,length(markers));
            elseif exist('smooth_img','var') == 1 && length(smooth_img) == 1 && i == 9
                smooth_img = repmat(smooth_img,1,length(markers));
            elseif exist('smooth_sigma','var') == 1 && isempty(smooth_sigma) && i == 10
                smooth_sigma = ones(1,length(markers));
            elseif exist('smooth_sigma','var') == 1 && length(smooth_sigma) == 1 && i == 10
                smooth_sigma = repmat(smooth_sigma,1,length(markers));
            elseif exist('DoG_img','var') == 1 && length(DoG_img) == 1 && i == 11
                DoG_img = repmat(DoG_img,1,length(markers));
            elseif exist('DoG_minmax','var') == 1 && isempty(DoG_minmax) == 1 && i == 12
                DoG_minmax = [0.8,1.5];
            elseif exist('DoG_factor','var') == 1 && length(DoG_factor) == 1 && i == 13
                DoG_factor = repmat(DoG_factor,1,length(markers));
            elseif exist('darkfield_intensity','var') == 1 && length(darkfield_intensity) == 1 && i == 14
                darkfield_intensity = repmat(darkfield_intensity,1,length(markers));
            elseif exist('update_intensity_channels','var') == 1 && isempty(update_intensity_channels) && i == 15
                update_intensity_channels = 1:length(markers);
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