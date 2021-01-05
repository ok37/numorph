% For each Top1 cKO, measure DICE and CoV for all cortical structures
cd('/home/ok37/tc_pipeline')
addpath(genpath(pwd))
load 'TCa_variables.mat'
cd(home_path)
config = load(fullfile(home_path,'TCa_variables.mat'));

s = 0.4;
dices = zeros(1,8);
use_inverse = 'false';
config.registration_method = 'p';
num_points = 200;

I_mask = gen_mask(config.hemisphere, config.structures_of_interest);
I_mask2 = imresize3(I_mask,0.4/s,'Method','nearest');

%save('I_mask.mat','I_mask')

vals = unique(I_mask);

for n = 1:4
    switch n 
        case 1
            output_directory = '/media/SteinLab4/TOP110R/output';
            config.points_file = 'native_landmarks110.csv';
        case 2
            output_directory = '/media/SteinLab4/TOP11L/output';
            config.points_file = 'native_landmarks11.csv';
        case 3
            output_directory = '/media/SteinLab4/TOP14R/output';
            config.points_file = 'native_landmarks14.csv';
        case 4
            output_directory = '/media/SteinLab4/TOP16R/output';
            config.points_file = 'native_landmarks16.csv';
        case 5
            output_directory = '/media/SteinLab4/WT7R/output';
                config.points_file = 'native_landmarks7.csv';
            %n = n -4;
        case 6
            output_directory = '/media/SteinLab5/WT11L/output';
                config.points_file = 'native_landmarks11.csv';
                %n = n -4;
        case 7
            output_directory = '/media/SteinLab5/WT1L/output';
                config.points_file = 'native_landmarks1.csv';
            %n = n -4;
        case 8
            output_directory = '/media/SteinLab5/WT8R/output';
                config.points_file = 'native_landmarks8.csv';
            %n = n -4;
    end

    if any(n == 1:4)
        use_inverse = 'true';
    else
        use_inverse = 'false';
    end
        disp(use_inverse)

    config.output_directory = output_directory;
    if isequal(use_inverse,'true')
        config.direction = 'img_to_atlas';
    else
        config.direction = 'atlas_to_img';
    end    

    %% Read Filename Information
    path_cell = {};
    % Load images from resampled directory
    fprintf('%s\t Reading image filename information from resampled directory \n',datetime('now'))
    if exist(fullfile(output_directory,'resampled'),'dir') == 7
        path_cell{1} = dir(fullfile(output_directory,'resampled'));
        location = 'resampled';
        path_table_resampled = path_to_table(path_cell,markers,channel_num,location,sample_name);
    else
        msg = sprintf('%s\t Cannot locate resampled directory in the specified image directory',datetime('now'));
        error(msg);
    end

    mov_img_path = path_table_resampled.file{1};
    config.save_registered_image = 'false';

    %% Evaluate Registration Accuracy
    % This registers image to atlas and would apply to registration parameters
    % to the true positive, manually traced mask. Not totally representitive of
    % true DICE but faster. To get true overlap of mask on native sampple, need
    % to calculate inverse transformation parameters and apply these to the
    % annotation image. Then measure DICE on manually traced mask without
    % transforming it. For points registration, specify how many points to test
    % in num_points.
    for nn = 1:length(num_points)
        % Run registration
        fprintf('\t Working on image: %s \n',mov_img_path) 

        % Find true positive, manually traced mask. This is at 10um resolution
        mask_location = fullfile(output_directory,'resampled');
        files = dir(mask_location);
        file_idx = arrayfun(@(s) contains(s.name,'cortex_mask.tif'), files);

        % Load manually traced mask and resize
        I_true = loadtiff(fullfile(files(file_idx).folder,files(file_idx).name));
        I_true = imresize3(I_true,0.4/s,'Method','nearest');

        % Perform registration
        reg_params = register_to_atlas(mov_img_path, config, num_points);
        
        % Adjust sizes and spacing
        size1 = reg_params.native_img_size;
        for j = 1:length(reg_params.atlas_to_img.TransformParameters)
            reg_params.atlas_to_img.TransformParameters{j}.FinalBSplineInterpolationOrder = 0;
            reg_params.atlas_to_img.TransformParameters{j}.Size = size1;
            reg_params.atlas_to_img.TransformParameters{j}.Spacing = [s, s, s];
        end

        I_mask = transformix(I_mask2,reg_params.atlas_to_img,[s, s, s], []);

        % Calculate DICE
        dices(n) = dice(logical(I_true),logical(I_mask));
        disp(dices)
	save(fullfile(output_directory,'variables','I_mask.mat'),'I_mask')
	%I_mask = bin_annotation_structures(I_mask);
        ssizes(:,n)  = histcounts(I_mask, max(vals)+1);
    end
    
    %cov = ssizes(all(ssizes>0,2),:);
    
    %cov2 = std(cov,[],2)./mean(cov,2);
    disp(mean(dices))
    disp(std(dices))
    %disp(mean(cov2))
    %disp(std(cov2))    
end
writematrix(ssizes, 'CoV_AB2.csv')
