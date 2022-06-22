function img = standardize_nii(img, in_res, in_or, in_hem, is_mask,...
    out_res, out_or, out_hem, out_type)
%--------------------------------------------------------------------------
% Standardize 3D image to match specified orientation. Additionally,
% convert between left/right/both hemispheres.
%--------------------------------------------------------------------------
% Usage:
% img = standardize_nii(img, in_res, in_or, in_hem, is_mask, out_res,...
%                       out_or, out_hem, out_type)
%
%--------------------------------------------------------------------------
% Inputs:
% img: 3D image matrix.
%
% in_res: Input voxel resolution.
%
% in_or: 3 character vector specifying input orientation. 
%
% in_hem: Input hemisphere. ('left', 'right', 'both')
%
% is_mask: Logical specifying whether input is a mask.
%
% out_res: Output voxel resolution. (default: 25)
%
% out_or: Output orientation. (default: 'ail')
%
% out_hem: Ouput hemisphere. (default: "left")
%
% out_type: Output data type. (default: same input data)
%
%--------------------------------------------------------------------------
% Outputs:
% img: Output 3D image.
%
%--------------------------------------------------------------------------

in_hem = string(in_hem); out_hem = string(out_hem);

% Defaults:
if nargin<6 || isempty(out_res) 
    out_res = 25;
end
if nargin<7 || isempty(out_or)
    if isequal(in_hem,'left')
        out_or = 'ail';
    elseif isequal(in_hem, 'right')
        out_or = 'pir';
    elseif isequal(in_hem, 'both')
        out_or = 'pls';
    else
        out_or = in_or;
    end
end
if nargin<8 || isempty(out_hem)
    out_hem = in_hem;
end
if nargin<9 || isempty(out_type)
    out_type = class(img);
end

res_adj = in_res/out_res;

% Image is flipped and rotated if double or single compared to uint16
% Standardized image is double
%if is_mask
    %if isequal(class(img),'double') || isequal(class(img),'single')
    %    img = imrotate(img,90);
    %    img = flip(img,1);
    %end
%    img = double(img);
%else
%    img = double(img);
%end
img = double(img);

% Transform atlas_img to match sample orientation
if ~isequal(in_hem, out_hem)
    img = adjust_hemisphere(img, in_hem, out_hem, in_or);
end

% Permute and resize
img = permute_orientation(img,char(in_or),out_or);

if ~is_mask
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


function img = adjust_hemisphere(img, in_hem, out_hem, orientation)
% Transform 3D hemisphere to match desired orientation

orientation = char(orientation);
mid_axis = find(orientation == 'l' | orientation == 'r');
if isequal(in_hem,"left")
    if isequal(out_hem,"right")
        img = flip(img,mid_axis);
    elseif isequal(out_hem,"both")
        img2 = flip(img,mid_axis);
        img = cat(3,img,img2);
    end
end

if isequal(in_hem,"right")
    if isequal(out_hem,"left")
        img = flip(img,mid_axis);
    elseif isequal(out_hem,"both")
        img2 = flip(img,mid_axis);
        img = cat(3,img,img2);
    end
end

if isequal(in_hem,"both")
    dims = size(img);
    mid = round(dims(mid_axis)/2);
    S.subs = repmat({':'},1,3);
    S.type = '()';
    if isequal(out_hem,"left")
        S.subs{mid_axis} = 1:mid;
        img = subsasgn(img,S,[]);
    elseif isequal(out_hem,"right")
        S.subs{mid_axis} = mid+1:dims(mid_axis);
        img = subsasgn(img,S,[]);
    end
end

end