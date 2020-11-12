function [z_displacement,ave_score] = z_align_channel(config,path_mov,path_ref,channel_idx)
%--------------------------------------------------------------------------
% Align .tif series from 2 channels along z dimension.
%--------------------------------------------------------------------------

% Set number of peaks and subixel value
peaks = 3;
usfac = 1;
max_shift = 0.1;
verbose = false;

% Unpack variables
z_positions = config.z_positions;
z_window = config.z_window;
lowerThresh = config.lowerThresh(channel_idx);
z_initial = config.z_initial(channel_idx);

% Adjust for resolution
res_adj = config.resolution{1}./config.resolution{channel_idx};
res_equal = all(res_adj == 1);

tempI = imread(path_ref.file{1});
shift_threshold = max(size(tempI))*max_shift;

if ~res_equal
    res_adj = size(tempI).*res_adj(1:2);
    shift_threshold = max(res_adj)*max_shift;
end

% Check lowerThresh
if lowerThresh<1
    lowerThresh = lowerThresh*65535;
end

% Read number of reference/moving images
nb_imgs_mov = height(path_mov);    
nb_imgs_ref = height(path_ref);

% Check file lengths
if nb_imgs_mov ~= nb_imgs_ref
%    error("Number of images do not match up")
end

% Pick reference images in range
z = linspace(0,nb_imgs_ref,z_positions+2);
z = round(z(2:end-1));
path_ref = path_ref(z,:);

% Initial matrix of cross correlation values and perform registration using
% cross correlation
cc = zeros(z_positions,z_window*2+1);
max_signal = zeros(1,length(z));
a=1;
for i = 1:length(z)
    % Read reference image and specify z range to look based on user-defined
    % window
    z_range = z(i)-z_window+z_initial:1:z(i)+z_window+z_initial;
    if any(z_range<1) || any(z_range>nb_imgs_mov)
        continue
    end
    
    ref_img = imread(path_ref.file{i});

    % Resample the reference image if the resolution is different
    % Note: we're resizing reference image instead of moving image for
    % speed and because moving image is usually at lower resolution
    if ~res_equal
        ref_img = imresize(ref_img, res_adj);
    end

    signal = zeros(1,length(z_range));
    b = 1;
    for j = z_range
        % Read moving image
        mov_img = imread(path_mov.file{j});

        % Calculate number of bright pixels
        signal(b)= numel(mov_img(mov_img>lowerThresh))/numel(mov_img);
        
        % If no bright pixels, continue to next image
        if signal(b)<0.01
            cc(a,b) = 0;
            continue
        end
        
        % Do simple crop in case image sizes don't match up
        if any(size(ref_img) ~= size(mov_img))
            mov_img = crop_to_ref(ref_img,mov_img);
        end    
        
        %Calculate phase correlation
        [~,~,~,type,cc(a,b)] = calculate_phase_correlation(mov_img,ref_img,peaks,usfac,shift_threshold);        

        %If pc didn't work, set cc to 0
        if type == 0
            cc(a,b) = 0;
        end

        %Check results
        %[pc_img,ref_img2,~,type,cc(a,b)] = calculate_phase_correlation(mov_img,ref_img,peaks,usfac,shift_threshold);
        %disp(cc)
        %if i == 1
            %figure
            %imshowpair(ref_img2*30,uint16(pc_img)*25)
        %end
        b = b+1;
    end
    max_signal(i) = max(signal);
    a = a+1;
end

% Display cross-correlation matrix
if verbose
    disp(cc)
end

% Generate score by multiplying intensity value by cross-correlation. This
% will lower effect caused by noise in images with no features
[val, pos] = max(cc'.*max_signal);
pos_idx = unique(pos);

% Sum scores for all unqiue positions
score = zeros(1,length(pos_idx));
for i = 1:length(pos_idx)
    score(i) = sum(val(pos==pos_idx(i)));    
end

% Max score determines final z displacement
[~,k] = max(score);

% Adjust final z displacement based on initial z position
z_displacement = pos_idx(k)-(z_window+1)+z_initial;
ave_score = mean(val);

end
