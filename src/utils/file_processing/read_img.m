function [img, metadata] = read_img(filepath,index,tiff_flag)
%--------------------------------------------------------------------------
% Wrapper function to read images of various types.
%--------------------------------------------------------------------------
% Usage: 
% img = read_img(filepath,index,tiff_flag);
%
%--------------------------------------------------------------------------
% Inputs:
% filepath: (string or cell) Path to file location.
%
% index: (1x4 int) Channel, z, y, x position in multi-image table.
% (default if table: [1,1,1,1])
%
% tiff_flag: (logical) Use loadtiff instead of MATLAB's native imread.
% loadtiff may work faster for larger file sizes. (default: false)
%
%--------------------------------------------------------------------------
% Outputs:
% img: (numeric) Target image.
%
% metadata: (structure) Image metadata availabel for certain file types. 
%
%--------------------------------------------------------------------------

if istable(filepath) && nargin<2
    index = ones(1,4);
elseif istable(filepath) && length(index)<4
    index = [index, ones(1,4-length(index))];
else
    index = [];
end

if nargin<3
    tiff_flag = false;
end

if ~isempty(index) && istable(filepath)
    filepath = filepath(filepath.channel_num == index(1) &...
        filepath.z == index(2),:);
    if height(filepath) > 1
        filepath = filepath(filepath.y == index(3) &...
            filepath.x == index(4),:);
    end
    if isempty(filepath)
        img = [];
        return
    else
        filepath = filepath.file{1};
    end
end

if iscell(filepath)
    img = read_img_worker(filepath{1},tiff_flag);
else
    img = read_img_worker(filepath,tiff_flag);
end

if iscell(img)
    metadata = img{2};
    img = img{1};
end
end


function img = read_img_worker(filepath,tiff_flag)

[~,~,ext] =  fileparts(filepath);

switch ext
    case {'.tif','.tiff'}
        if ~tiff_flag
            info = imfinfo(filepath);
        end
        if tiff_flag || info.FileSize >1E7
            img = loadtiff(filepath);
        else
            img = imread(filepath);
        end
    case {'.nii','.mhd'}
        img{1} = niftiread(filepath);
        img{2} = niftiinfo(filepath);

    case {'.nrrd','.nhdr'}
        out = nhdr_nrrd_read(filepath,true);        
        img{1} = out.data;
        out = rmfield(out,'data');
        img{2} = out;
    otherwise
        error("Unrecognized file type.")
end

end