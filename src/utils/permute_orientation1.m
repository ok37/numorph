function img = permute_orientation1(img, or_in, or_out)
%--------------------------------------------------------------------------
% Permute sample anotomical orientation. Orientations are specified as 1x3 
% character vectors. Row, column, slice at idx=1 in the image should match 
% 1st, 2nd, 3rd characters in each vector in this order.
%
% Orientation key: 
% anterior(a)/posterior(p), 
% superior(s)/inferior(i), 
% left(l)/right(r)
%--------------------------------------------------------------------------
% Usage:
% img = permute_orientation(img, or_in, or_out)
%
% Example:
% img = permute_orientation(img,'rsa','ail');
% Permutes image from 'right','superior','anterior' at (y=1,x=1,z=1) to
% 'anterior','inferior','left'
%--------------------------------------------------------------------------
% Inputs:
% img: (numeric) Input 3D image. 
%
% or_in: (1x3 char) Input image orientation.
%
% or_out: (1x3 char) Output image orientation.
%--------------------------------------------------------------------------

% Check orientation characters
assert(length(or_in) == 3, "Input/output orientation should be 1x3 character arrays")
assert(length(or_out) == 3, "Input/output orientation should be 1x3 character arrays")

% Check if input/output are the same
if string(or_in) == string(or_out)
    return
end

% Scan characters to see which to permute
yxz_in = 1:3;
for i = 1:3
    if ismember({or_in(i)},{'a','p'})
        if contains(or_out,'p')
            idx = find(or_out == 'p');
            img = flip(img,idx);
        end
        yxz_in(i) = find(or_out == 'a' | or_out == 'p');

    elseif ismember({or_in(i)},{'s','i'})
        if contains(or_out,'s')
            idx = find(or_out == 's');
            img = flip(img,idx);
        end
        yxz_in(i) = find(or_out == 'i' | or_out == 's');
        
    elseif ismember({or_in(i)},{'l','r'})
        if contains(or_out,'l')
            idx = find(or_out == 'l');
            img = flip(img,idx);
        end
        yxz_in(i) = find(or_out == 'l' | or_out == 'r');

    else
        error("Unrecognized orientation character specified")
    end
end
assert(length(unique(yxz_in)) == 3,"All 3 axes not specified correctly")

% Permute
img = permute(img,yxz_in);

end
