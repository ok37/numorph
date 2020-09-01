function predict_centroids_3dunet(config)
%--------------------------------------------------------------------------
% Wrapper function to run nuclei centroid predictions via 3D-Unet. Requires
% installation of the necessary python packages located in environment.yml
% via conda. 
%--------------------------------------------------------------------------

% Set pythonpath
pythonpath = fullfile(config.home_path, 'analysis','3dunet');
setenv('PYTHONPATH',pythonpath)

% Run prediction
command = sprintf("source activate 3dunet-centroid; python %s -m",...
    fullfile(pythonpath,'nuclei','predict_validation.py'));
[status,cmdout] = system(command,'-echo');
if status ~= 0 
    error(cmdout)
end


end