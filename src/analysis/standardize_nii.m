function img = standardize_nii(img, img_type, resolution, orientation, hemisphere, out_type)

% Defaults:
register_resolution = 25;
register_orientation = 'ail';

if nargin<6 
    out_type = class(img);
end

res_adj = resolution/register_resolution;
% Image is flipped and rotated if double or single compared to uint16
% Standardized image is double
if ~isequal(img_type,"mask")
    if isequal(class(img),'double') || isequal(class(img),'single')
        img = imrotate(img,90);
        img = flip(img,1);
    end
    img = double(img);
end
    
% Transform atlas_img to match sample orientation
if isequal(img_type,"atlas") || isequal(img_type,"mask")
    if isequal(hemisphere, "right")
        img = flip(img,1);
    elseif isequal(hemisphere,"both")
        atlas_img2 = flip(img,3);
        img = cat(3,img,atlas_img2);
    elseif isequal(hemisphere,"left") 
    elseif ~isequal(hemisphere,"none")
        error("Unrecognized sample hemisphere value selected")
    end
end
    
% Permute and resize
img = permute_orientation(img,char(orientation),register_orientation);

if ~isequal(img_type,"mask")
    img = imresize3(img,round(res_adj.*size(img)));
    
    % Rescale intensities
    min_int = min(img,[],'all');
    max_int = max(img,[],'all');
    img = 65535*(img - min_int)/(max_int-min_int);
    
    % Return output type
    if isequal(out_type,'double')
    elseif isequal(out_type,'uint16')
        img = uint16(img);
    elseif isequal(out_type,'single')
        img = single(img);
    end
else
    img = imresize3(img,round(res_adj.*size(img)),'Method','nearest');
end

end