function mask = generate_sampling_mask(I, mask_int_threshold,lowerThresh,upperThresh)

% Generate intensity threshold if not provided
%if isempty(mask_int_threshold)
%    mask_int_threshold(1) = (lowerThresh(1)+upperThresh(1))/2;
%    mask_int_threshold(2) = otsuthresh(im2uint16(mask(:)))*0.9;
%    mask_int_threshold(3) = mask_int_threshold1+mask_int_threshold2/2;
%end

%mask_int_threshold = 0.05;
    
% Generate mask
% Downsample to 20% resolution
I2 = imresize3(I,0.20,'linear'); 

if isempty(mask_int_threshold)
    mask_int_threshold = graythresh(I2)*0.75;
end

% Binarize mask and fill holes
mask = imbinarize(I2,mask_int_threshold);
mask = imfill(mask,26,'holes');
    
% Keep only largest compnent. Disconnected components will give errors
labels = bwconncomp(mask);
sizes = cellfun(@(s) length(s), labels.PixelIdxList);
[~, idx] = max(sizes);
for i = 1:length(sizes)
    if i ~= idx
        mask(labels.PixelIdxList{i}) = 0;
    end
end
    
% Resize back to original size
mask = imresize3(double(mask), size(I),'method', 'nearest');
signal = sum(mask(:))/numel(mask);

fprintf('\n Using a mask intensity threshold of %.4f \n', mask_int_threshold)
fprintf('\n Using %.4f percent of pixels \n', signal*100)

end
