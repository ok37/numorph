function I_adj = dog_adjust(I,nuc_radius,weight)
%--------------------------------------------------------------------------
% Enhacen blob objects using Difference of Gaussian filter
%--------------------------------------------------------------------------

% Use very small filter size
factor = 0.5;
s = factor*([nuc_radius*0.8 nuc_radius*1.2]);
filter_size = 2*ceil(s(2))+1;

%Check to see if single
img_class = class(I);
if ~isequal(img_class,'single')
    I = single(I);
end

% Instead of subtracting small from large and detecting edges, subtract
% large from small to detect blobs

I = I.^1.25;
I = imgaussfilt(I,2);
I2 = imgaussfilt(I,s(1),'FilterSize',filter_size);
I3 = imgaussfilt(I,s(2),'FilterSize',filter_size);
dog = I3-I2;
%dog = dog + min(dog(:));
%imshow(imadjust(uint16(dog)))
% Adjusted image is mean of original image and DoG image scaled by some
% weight
scale = (max(I(:))-min(I(:)))/(max(dog(:))-min(dog(:)));

if nargin <3
    weight = 1;
end

dog = dog*weight/scale;
%dog(dog<0) = 0;
I_adj = I-dog;

% Recast to original class if necessary
if ~isequal(I_adj, img_class)
    I_adj = cast(I_adj,img_class);
end

imshow(imadjust(uint16(I_adj),[0 0.02], [0 1],1.25))

end