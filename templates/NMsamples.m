function [img_directory, output_directory, group] = NMsamples(sample, save_flag)
%--------------------------------------------------------------------------
% Record sample specific information here. Additionally, any default
% processing and/or analysis parameters can be overwritten as this function
% is called after defaults are loaded.
%--------------------------------------------------------------------------

switch sample
    case 'TEST1'
        img_directory = "/ssd2/userdata/ok37/test/Ctip2-ToPro";
        output_directory = "/ssd2/userdata/ok37/test/output";
        sample_name = "TEST1";
        group = "TEST1";
        markers = ["ToPro", "Ctip2"];
        channel_num = ["C01","C00"];
        resolution = [1.21, 1.21, 4];
        ls_width = [50 50 40];
        overlap = 0.15;
    case 'WT1L'
        img_directory(1) = "/media/SteinLab3/WT1L/High/Ctip2-ToPro";
        img_directory(2) = "/media/SteinLab3/WT1L/High/Cux1";
        output_directory = "/media/SteinLab3/WT1L/output";
        sample_name = "WT1L";
        group = 'WT';
        markers = ["ToPro","Ctip2","Cux1"];
        channel_num = ["C01", "C00", "C00"];
        resolution = [1.21, 1.21, 4]; 
        ls_width = [50 50 40];             
        laser_y_displacement = [0 0 -0.025];
        overlap = 0.15;
        points_file = 'native_landmarks1.csv';
    case 'WT11L'
        img_directory(1) = "/media/SteinLab3/WT11L/Ctip2-ToPro";
        img_directory(2) = "/media/SteinLab3/WT11L/Cux1";
        output_directory = "/media/SteinLab3/WT11L/output";
        sample_name = "WT11L";
        group = 'WT';
        markers = ["ToPro","Ctip2","Cux1"];
        channel_num = ["C01", "C00", "C00"];
        resolution = [1.21, 1.21, 4]; 
        ls_width = [50 50 40];             
        laser_y_displacement = [0 0 0];
        overlap = 0.15;
        points_file = 'native_landmarks11.csv';
    case 'WT8R'
        img_directory(1) = "/media/SteinLab3/WT8R/High/Ctip2-ToPro";
        img_directory(2) = "/media/SteinLab3/WT8R/High/Cux1";
        output_directory = "/media/SteinLab3/WT8R/output";
        sample_name = "WT8R";
        group = 'WT';
        markers = ["ToPro","Ctip2","Cux1"];
        channel_num = ["C01", "C00", "C00"];
        resolution = [1.21, 1.21, 4]; 
        ls_width = [50 50 30];             
        laser_y_displacement = [0 0 -0.025];
        overlap = 0.15;
        points_file = 'native_landmarks8.csv';
    case 'WT7R'
        img_directory(1) = "/media/SteinLab4/WT7R/High/Ctip2-ToPro";
        img_directory(2) = "/media/SteinLab4/WT7R/High/Cux1";
        output_directory = "/media/SteinLab4/WT7R/output";
        sample_name = "WT7R";                   
        group = 'WT';
        markers = ["topro","ctip2","cux1"];   
        channel_num = ["C01", "C00", "C00"];
        resolution = [1.21, 1.21, 4]; 
        ls_width = [50 50 30];            
        laser_y_displacement = [0 0 0];
        overlap = 0.15;
        points_file = 'native_landmarks7.csv';
    case 'TOP110R'
        img_directory(1) = "/media/SteinLab4/TOP110R/High/Ctip2-ToPro";
        img_directory(2) = "/media/SteinLab4/TOP110R/High/Cux1";
        output_directory = "/media/SteinLab4/TOP110R/output";
        sample_name = "TOP110R";
        group = 'TOP1';
        markers = ["ToPro","Ctip2","Cux1"];
        channel_num = ["C01", "C00", "C00"];
        resolution = [1.21, 1.21, 4]; 
        ls_width = [50 50 30];            
        laser_y_displacement = [0 0 0];
        overlap = 0.15;
        points_file = 'native_landmarks110.csv';
    case 'TOP11L'
        img_directory(1) = "/media/SteinLab4/TOP11L/High/Ctip2-ToPro";
        img_directory(2) = "/media/SteinLab4/TOP11L/High/Cux1";
        output_directory = "/media/SteinLab4/TOP11L/output";
        sample_name = "TOP11L";
        group = 'TOP1';
        markers = ["ToPro","Ctip2","Cux1"];
        channel_num = ["C01", "C00", "C00"];
        resolution = [1.21, 1.21, 4]; 
        ls_width = [50 50 40];
        laser_y_displacement = [0 0 0];
        overlap = 0.15;
        points_file = 'native_landmarks11.csv';
    case 'TOP14R'
        img_directory(1) = "/media/SteinLab4/TOP14R/High/Ctip2-ToPro";
        img_directory(2) = "/media/SteinLab4/TOP14R/High/Cux1";
        output_directory = "/media/SteinLab4/TOP14R/output";
        sample_name = "TOP14R";
        group = 'TOP1';
        markers = ["ToPro","Ctip2","Cux1"];
        channel_num = ["C01", "C00", "C00"];
        resolution = [1.21, 1.21, 4]; 
        ls_width = [50 50 40];
        laser_y_displacement = [0 0 0];
        overlap = 0.15;
        points_file = 'native_landmarks14.csv';
    case 'TOP16R'
        img_directory(1) = "/media/SteinLab4/TOP16R/High/ToPro";
        img_directory(2) = "/media/SteinLab4/TOP16R/High/Ctip2";
        img_directory(3) = "/media/SteinLab4/TOP16R/High/Cux1";
        output_directory = "/media/SteinLab4/TOP16R/output";
        sample_name = "TOP16R";
        group = 'TOP1';
        markers = ["topro","ctip2","cux1"];
        ls_width = [50 50 40];
        laser_y_displacement = [0 0 0];
        resolution = [1.21, 1.21, 4]; 
        channel_num = ["C00", "C00", "C00"];
        overlap = 0.15;
        points_file = 'native_landmarks16.csv';
    otherwise
        error("Sample %s does not exist in NMsamples.",sample)
end

% Append sample info to variable structure
if save_flag
    save(fullfile('templates','NM_variables.mat'),'-mat','-append')
end
end
