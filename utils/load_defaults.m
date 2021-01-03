function load_defaults(stage)
%--------------------------------------------------------------------------
% Overwrite template files with default configuration.
%--------------------------------------------------------------------------
% Inputs:
%
% stage - NuMorph stage to load (i.e.'process', 'analyze', 'evaluate').
%--------------------------------------------------------------------------

% Prompt user what stage
if nargin<1
    stage = input("Choose which stage to load (i.e. 'process',"+...
        " 'analyze', 'evaluate'): \n");    
end

switch stage
    case 'process'
        
        
    case 'analyze'
        
        
    case 'evaluate'
        
        
    otherwise
        error("Unrecognized stage selected")
end

usr_confirmation = input(sprintf("Defaults for stage %s have been loaded. Type 'yes' "+...
    "to continue and overwrite or press enter to quit.\n",stage));  

    
end