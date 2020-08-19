function I = subtract_background(I,se,t_adj,y_adj,lowerThresh,upperThresh,gamma,filter)

%Check to see if single
img_class = class(I);
if ~isequal(img_class,'single')
    I = single(I);
end

%Filter image if specified
if isequal(filter,"true")
    I = imgaussfilt(I,0.2,'FilterSize',7,'FilterDomain','spatial');
end

% Apply any intensity adjustment
%I = apply_intensity_adjustment(I,t_adj,y_adj);

%Subtract background
%Morphological opening to get the background
nuc_radius

%Subtract. Here I chose to do only 95% of background to avoid 0 pixel
%values.
I = I - back;

%Use imadjust to remap low/high intensities and apply any gamma
I = imadjust(uint16(I),[0 upperThresh-lowerThresh],[],gamma);

if ~isequal(I, img_class)
    I = cast(I,img_class);
end

end
