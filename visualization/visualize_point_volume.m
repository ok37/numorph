function [p2, scat, df, I] = visualize_point_volume(df, I)
% Constants
rf = 0.25;  % Rescale factor
res1 = [1.208,1.208, 4];    % High resolution (centroids)
res2 = [10, 10, 10];        % Low resolution (mask)
res_view = [50,50,50];   % Resolution to view at
remove_l1 = "true";         % Remove layer 1

if nargin < 1
    %df = readmatrix(fullfile('evaluation','Results_svm','WT11L_centroids.csv'));
    df = readmatrix(fullfile('evaluation','Results_svm','TOP11L_centroids.csv'));

end

if nargin < 2
    %I = niftiread(fullfile('visualization','WT11L_C1_ToPro_resampled.nii'));
    I = niftiread(fullfile('visualization','TOP11L_C1_ToPro_resampled.nii'));
    %I = niftiread(fullfile('visualization','TOP16R_C1_ToPro_resampled.nii'));
    I = imresize3(I, rf);
end

%I = imresize3(I,round((res2./res_view).*size(I)));

df2 = df;
%df2(:,1:3) =df2(:,1:3).*(res1./res_view);
df2(:,1:3) =round(df2(:,1:3)*rf).*(res1./res2);

%I = permute(I,[3,1,2]);

%pad_amount = round(size(I)*0.2);
%I = padarray(I,pad_amount,'pre');
%pad_amount(1) = 0;
%I = padarray(I,pad_amount,'post');
%pad2 = pad_amount([2 1 3]);
%df2(:,1:3) = df2(:,1:3) + pad2;


[nrows, ncols, nslices] = size(I);
r1 = round(0.45*[nrows, ncols, nslices]);
r2 = round(0.55*[nrows, ncols, nslices]);

chunk = single(I(r1(1):r2(1), r1(2):r2(2), r1(3):r2(3)));
lowerThresh = (prctile(single(I(:)),2) + prctile(chunk(:),2))/(2*65535);

%
BW = imbinarize(I,lowerThresh);
se = strel('disk',6);
for i = 1:size(BW,3)
    %BWe = imerode(BW(:,:,i),se);
    %BW(:,:,i) = imreconstruct(BWe,BW(:,:,i));
    BW(:,:,i) = imopen(BW(:,:,i),se);
    BW(:,:,i) = imclose(BW(:,:,i),se);
    BW(:,:,i) = imfill(BW(:,:,i),'holes');
end

%BW = smooth3(BW);

%imagesc(BW(:,:,10))

%
cc = bwconncomp(BW,6);
cc_vol = cellfun(@(s) numel(s),cc.PixelIdxList);
max_cc = find(max(cc_vol));

BW(:) = 0;
BW(cc.PixelIdxList{max_cc}) = 1;

% Make the volume
s = isosurface(BW,0);
s = reducepatch(s,0.1);

s2 = smoothpatch(s,1,5,1);

s3 = bin_annotation_structures(df2(:,4),'layers');
rm_idx = ismember(s3,[0,1]);
df3 = df2(~rm_idx,:);

[~,idx] = sort(df3(:,1));
df3 = df3(idx,:);

% Plotting
[scat] = make_plot(df2,s2,BW,2,1);
\a = 1;


end

function keepAlpha(src,eventData,cEdgeColor)  
    scat = src.MarkerHandle;
    scat.EdgeColorData = EdgeColor;
    scat.FaceColorData = FaceColor;   
end

