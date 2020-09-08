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

% Make sure path is set
addpath(genpath('.'))

% Load variables
switch stage
    case 'process'
        NMp_template
        temp_path = fileparts(which('NMp_template'));
        addpath(genpath(fullfile(temp_path,'..')))
        home_path = fileparts(which('NM_process.m'));
    case 'analyze'
        NMa_template
        temp_path = fileparts(which('NMa_template'));
        addpath(genpath(fullfile(temp_path,'..')))
        home_path = fileparts(which('NM_analyze.m'));
    case 'evaluate'
        NMe_template
        temp_path = fileparts(which('NMe_template'));
        addpath(genpath(fullfile(temp_path,'..')))
        home_path = fileparts(which('NM_evaluate.m'));
    otherwise
        error("Invalid input")
end

% Save config structure
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
        [~, centroids_directory(i), group(i)] = NMsamples(sample_idx{i}(2:end-1));
    end
    fclose(fid);
    output_directory = results_directory;
    save(fullfile('templates','NM_variables.mat'),{'centroids_directory','group','output_directory'},'-mat','-append')
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