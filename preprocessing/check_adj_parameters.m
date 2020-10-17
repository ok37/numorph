function [adj_params, config, lowerThresh, upperThresh] = check_adj_parameters(adj_params, config, lowerThresh_measured, upperThresh_measured)
% This function checks for user defined adjustment parameters and will
% override any measured. Some of these variables may be required for
% subsequent steps such as image stitching

if nargin <3
    lowerThresh_measured = adj_params.lowerThresh;
    upperThresh_measured = adj_params.upperThresh;
end
        
% Check for user-defined intensity threshold values
lowerThresh = config.lowerThresh;
upperThresh = config.upperThresh;

% If given is 16-bit integer, rescale to unit interval
if lowerThresh>1
    lowerThresh = lowerThresh/65535;
end
if upperThresh>1
    upperThresh = upperThresh/65535;
end
if lowerThresh_measured>1
    lowerThresh_measured = lowerThresh_measured/65535;
end
if upperThresh_measured>1
    upperThresh_measured = upperThresh_measured/65535;
end

% If empty, set to predicted threshold values
if isempty(lowerThresh)
    fprintf('%s\t Using measured lowerThresh of %s \n',datetime('now'),num2str(round(lowerThresh_measured*65535)));
    adj_params.lowerThresh = lowerThresh_measured;
    config.lowerThresh = lowerThresh_measured;
else
    adj_params.lowerThresh = lowerThresh;
end

if isempty(upperThresh)
    fprintf('%s\t Using measured upperThresh of %s \n',datetime('now'),num2str(round(upperThresh_measured*65535)));    
    adj_params.upperThresh = upperThresh_measured;
    config.upperThresh = upperThresh_measured;
else
    adj_params.upperThresh = upperThresh;
end

% Set gamma equal to 1 if undefined
if isempty(config.gamma)
    adj_params.gamma = ones(1,length(config.markers));
else
    adj_params.gamma = config.gamma; 
end

% Set default nucleus diameter
if isempty(config.nuc_radius) && isequal(config.subtract_background,'true')
    fprintf('%s\t No nucleus radius indicated, using default 15um \n',datetime('now'));
    config.nuc_radius = round(15/config.resolution(1));
end

end