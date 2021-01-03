function [adj_params, config, lowerThresh, upperThresh, signalThresh] = check_adj_parameters(adj_params, config, lowerThresh_measured, upperThresh_measured, signalThresh_measured)
%--------------------------------------------------------------------------
% Check for user defined adjustment parameters that will override measured
% values.
%--------------------------------------------------------------------------
% 
% Inputs:
% adj_params - structure containing adjustment parameters.
%
% config - config structure from NM_process.
%
% lowerThresh_measured - 
%
% upperThresh_measured - 
%
% signalThresh_measured - 
%
% Ouputs:
%--------------------------------------------------------------------------


if nargin <3
    lowerThresh_measured = adj_params.lowerThresh;
    upperThresh_measured = adj_params.upperThresh;
    signalThresh_measured = adj_params.signalThresh;
end
        
% Check for user-defined intensity threshold values
lowerThresh = config.lowerThresh;
signalThresh = config.signalThresh;
upperThresh = config.upperThresh;

% If given is 16-bit integer, rescale to unit interval
if lowerThresh>1
    lowerThresh = lowerThresh/65535;
end
if upperThresh>1
    upperThresh = upperThresh/65535;
end
if signalThresh>1
    signalThresh = signalThresh/65535;
end
if lowerThresh_measured>1
    lowerThresh_measured = lowerThresh_measured/65535;
end
if upperThresh_measured>1
    upperThresh_measured = upperThresh_measured/65535;
end
if signalThresh_measured>1
    signalThresh_measured = signalThresh_measured/65535;
end

% If empty, set to predicted threshold values
if isempty(lowerThresh)
    fprintf('%s\t Using measured lowerThresh of %s \n',datetime('now'),num2str(round(lowerThresh_measured*65535)));
    adj_params.lowerThresh = lowerThresh_measured;
    config.lowerThresh = lowerThresh_measured;
else
    adj_params.lowerThresh = lowerThresh;
    config.lowerThresh = lowerThresh;
end

if isempty(upperThresh)
    fprintf('%s\t Using measured upperThresh of %s \n',datetime('now'),num2str(round(upperThresh_measured*65535)));    
    adj_params.upperThresh = upperThresh_measured;
    config.upperThresh = upperThresh_measured;
else
    adj_params.upperThresh = upperThresh;
    config.upperThresh = upperThresh;
end

if isempty(signalThresh)
    fprintf('%s\t Using measured signalThresh of %s \n',datetime('now'),num2str(round(signalThresh_measured*65535)));    
    adj_params.signalThresh = signalThresh_measured;
    config.signalThresh = signalThresh_measured;
else
    adj_params.signalThresh = signalThresh;
    config.signalThresh = signalThresh;
end

% Set Gamma equal to 1 if undefined
if isempty(config.Gamma)
    adj_params.Gamma = ones(1,length(config.markers));
else
    adj_params.Gamma = config.Gamma; 
end

% Set default nucleus diameter (15um)
if isempty(config.nuc_radius)
    config.nuc_radius = round(15/config.resolution(1));
end

% Update which adjustment to apply based on current configs
adj_params.adjust_tile_shading = config.adjust_tile_shading;
adj_params.adjust_tile_position = config.adjust_tile_position;

end