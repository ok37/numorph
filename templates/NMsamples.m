function [img_directory, output_directory, group] = NMsamples(sample, save_flag)
%--------------------------------------------------------------------------
% Record sample specific information here. Additionally, any default
% processing and/or analysis parameters can be overwritten as this function
% is called after defaults are loaded.
%--------------------------------------------------------------------------

switch sample
    case 'TEST1'
        %img_directory(1) = "/users/Oleh/test_images/output/stitched";  % Input image directory (required)
        %img_directory = "/Users/Oleh/Documents/MATLAB/numorph/test1.csv";  % Input image directory (required)
        %img_directory(2) = "/users/Oleh/test_images/Ctip2-downsampled";
        img_directory = "/users/Oleh/test_images/Ctip2-ToPro";

        output_directory = "/users/Oleh/test_images/output";    % Ouput directory (required)
        sample_name = "TEST1";                                  % Sample name/id (required)
        group = "TEST";                                         % Group name/id (optional: used for evaluation steps)
        position_exp = ["[\d*", "\d*]","Z\d*"];                 % 3 element string containing regular expression to find positions. Row(y), column(x), slice(z) (required)
        %channel_num = ["C1","C2"];                            % Channel id. Should correspond to imaged marker (optional: only if multiple channels are mixed in the same directory)
        channel_num = ["C01","C00"];
        
        markers = ["topro", "ctip2"];                           % Name of markers present (required). If directories contain multiple channels and channel ids are not unique, markers must be present in filename
        resolution = [1.21, 1.21, 4];                           % Image reolution in um/voxel (required)
        ls_width = [50 50 50];                                  % Light sheet width (optional: only required from measured intensity adjustment)
        overlap = 0.15;                                         % Overlap between tiles (optional: only required for stitching)
    otherwise
        error("Sample %s does not exist in NMsamples.",sample)
end

%--------------------------------------------------------------------------
% Do not edit
% Append sample info to variable structure
if save_flag
    save(fullfile('templates','NM_variables.mat'),'-mat','-append')
end
end
