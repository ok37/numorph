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
check_variable_lengths

% Make an output directory
if exist(output_directory,'dir') ~= 7
    mkdir(output_directory);
end

% Make a variables directory
if exist(fullfile(output_directory,'variables'),'dir') ~= 7
    mkdir(fullfile(output_directory,'variables'))
end

% Reload config
if nargout == 1
    config = load(fullfile('templates','NM_variables.mat'));
end

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
        path2 = [path2, ':', elastix_path_lib];
        if ~contains(path2,'elastix')
            setenv('DYLD_LIBRARY_PATH',path2)
        end
    else
        path2 = getenv('LD_LIBRARY_PATH');
        path2 = [path2, ':', elastix_path_lib];
        if ~contains(path2,'elastix')
            setenv('LD_LIBRARY_PATH',path2)
        end
    end
end

% Run
if nargin>2 && run
    switch stage
        case 'process'
            NM_process
        case 'analyze'
            NM_analyze
        case 'evaluate'
            NM_evaluate
    end
end
end


function check_variable_lengths
% Check user input variable lengths to make sure they're the correct length
% and match number of markers in most cases

% Variables to check
variable_names = {'markers','single_sheet','ls_width','laser_y_displacement','blending_method'};
load(fullfile('templates','NM_variables.mat'),variable_names{:});

for i = 1:length(variable_names)
    if exist('single_sheet','var') == 1 && length(single_sheet) == 1
        single_sheet = repmat(single_sheet,1,length(markers));
    elseif exist('ls_width','var') == 1 && length(ls_width) == 1
        ls_width = repmat(ls_width,1,length(markers));
    elseif exist('laser_y_displacement','var') == 1 && length(laser_y_displacement) == 1
        laser_y_displacement = repmat(laser_y_displacement,1,length(markers));
    elseif exist('blending_method','var') == 1 && length(blending_method) == 1
        blending_method = repmat(blending_method,1,length(markers));
    end
end

% Resave these variables
save(fullfile('templates', 'NM_variables.mat'),variable_names{:},'-mat','-append')

end