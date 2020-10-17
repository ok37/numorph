function [I_new, B] = smooth_background_subtraction(I, apply_filter, r)
%--------------------------------------------------------------------------
% Background subtraction with smoothing for images with low signal to
% noise. More similar to Fiji's background subtraction module but still
% uses a flat morphological element.
%--------------------------------------------------------------------------

% Default arguments
if nargin<2
    apply_filter = 'false';
end

if nargin<3
    r=9;
end

% Check to see if uint16
img_class = class(I);
if ~isequal(img_class,'uint16')
    I = uint16(I);
end
    
% Define structuring element
se = strel('disk',r);

% Apply average filter to background
h = 1/3*ones(3,1);
h = h*h';
B1 = imfilter(I,h);

% Apply max filter
%B1 = imdilate(B1,ones(3));

% Apply morphological opening to subtract
B = imopen(B1,se);
I_new = I - B;

% Optional: filter final results
if isequal(apply_filter,'true')
    I_new = wiener2(I_new,[4 4]);
    I_new = medfilt2(I_new,[3,3]);
end

% Recast to original class if necessary
if ~isequal(I_new, img_class)
    I_new = cast(I_new,img_class);
end

end