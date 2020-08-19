function [adj_params, config, lowerThresh, upperThresh] = check_adj_parameters(adj_params, config, lowerThresh_measured, upperThresh_measured)
% This function checks for user defined adjustment parameters and will
% override any measured. Some of these variables may be required for
% subsequent steps such as image stitching

if nargin <3
    lowerThresh_measured = adj_params.lowerThresh;
    upperThresh_measured = adj_params.upperThresh;
    gamma = [];
end
        
% Check for user-defined intensity threshold values
lowerThresh = config.lowerThresh;
upperThresh = config.upperThresh;

% If empty, set to predicted threshold values
if isempty(lowerThresh)
    fprintf('%s\t Using measured lowerThresh of %s \n',datetime('now'),num2str(lowerThresh_measured));
    adj_params.lowerThresh = lowerThresh_measured/65535;
else
    adj_params.lowerThresh = lowerThresh/65535;
end

if isempty(upperThresh)
    fprintf('%s\t Using measured upperThresh of %s \n',datetime('now'),num2str(upperThresh_measured));    
    adj_params.upperThresh = upperThresh_measured/65535;
else
    adj_params.upperThresh = upperThresh/65535;
end

% Set gamma equal to 1 if undefined
if isempty(config.gamma)
    adj_params.gamma = ones(1,length(config.markers));
else
    adj_params.gamma = config.gamma; 
end

% Set default nucleus diameter
if isempty(config.nuc_radius) && isequal(config.subtract_background,'true')
    fprintf('%s\t No nucleus radius indicated, using default 10um \n',datetime('now'));
    config.nuc_radius = round(10/config.resolution(1));
end

end