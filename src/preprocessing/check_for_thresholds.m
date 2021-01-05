function config = check_for_thresholds(config,path_table)
%--------------------------------------------------------------------------
% Check for intensity thresholds (i.e. thresholds.mat) structure and attach
% to config.
%--------------------------------------------------------------------------
% Usage: 
% config = check_for_thresholds(config,path_table)
% 
% Inputs:
% config - config structure from NM_process.
%
% path_table - table of file paths.
%
% Outputs:
% config = config structure with thresholds attached.
%--------------------------------------------------------------------------

var_file = fullfile(config.output_directory,'variables','thresholds.mat');

% First check if adj_params exists in variables folder
if exist(var_file,'file') == 2
    load(var_file,'thresholds')
    if all(thresholds.markers == config.markers) && all(thresholds.img_directory == config.img_directory)
        % Save into config
        config.lowerThresh = thresholds.lowerThresh;
        config.signalThresh = thresholds.signalThresh;
        config.upperThresh = thresholds.upperThresh;
        return
    end
end

fprintf("%s\t Lower and upper intensity thresholds are unspecified but are required "+...
    "for processing. Measuring these now... \n",datetime('now'));
for i = 1:length(config.markers)
    [thresholds.lowerThresh(i), thresholds.upperThresh(i), thresholds.signalThresh(i)] = measure_images(config,path_table,i,true);
end
thresholds.img_directory = config.img_directory;
thresholds.markers = config.markers;
save(var_file,'thresholds')

% Save into config
config.lowerThresh = thresholds.lowerThresh;
config.signalThresh = thresholds.signalThresh;
config.upperThresh = thresholds.upperThresh;

end