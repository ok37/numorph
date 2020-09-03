function predict_centroids_3dunet(config)
%--------------------------------------------------------------------------
% Wrapper function to run nuclei centroid predictions via 3D-Unet. Requires
% installation of the necessary python packages located in environment.yml
% via conda. 
%--------------------------------------------------------------------------

% Save matlab's config structure in the directory
save_path = fullfile(config.home_path,'analysis','3dunet','config.mat');
config = convertContainedStringsToChars(config);
save(save_path,'config')

% Set pythonpath
pythonpath = fullfile(config.home_path, 'analysis','3dunet');
setenv('PYTHONPATH', pythonpath)

% Run prediction
command = sprintf("source activate 3dunet-centroid; "+...
    "export PYTHONPATH=%s; "+...
    "python %s --mat %s",...
    pythonpath,fullfile(pythonpath,'nuclei','generate_chunks.py'),save_path);

%command = sprintf("source activate 3dunet-centroid; echo $PATH");

[status] = system(command,'-echo');

end
