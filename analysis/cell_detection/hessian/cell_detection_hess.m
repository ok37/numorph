% Cell detection
clear
%iname{1} = 'Training-Side0-[01x01].nii';
iname{1} = 'Training-Side60-[01x01].nii';
iname{2} = 'Training-Top0-[01x01].nii';
%iname{4} = 'Training-Top0-[02x02].nii';
%iname{5} = 'Training-Top60-[01x01].nii';

%lname{1} = 'OKu-Training-Side0-[01x01]-OKu2-ZH_7_8_lbl.nii';
lname{1} = 'OKu-Training-Side60-[01x01]-OKu-TL-SB_7_4_lbl.nii';
lname{2} = 'OKu-Training-Top0-[01x01]-OKu-EHF_5_11_lbl.nii';
%lname{4} = 'OKu-Training-Top0-[02x02]-OKu-EHF_8_31_lbl.nii';
%lname{5} = 'OKu-Training-Top60-[01x01]-OKu-SF-JM_7_19_lbl.nii';

for i = 1:length(lname)
   V = load_nii(lname{i});
    L = V.img;
    
    true_objects(i) = length(unique(L(:)))-1;
end






path1 = '/Users/Oleh/Documents/MATLAB/TCPipeline/Cell Detection';
diameter = 14;

for i = 1:length(iname)
resolution = [1 1 0.3];
threshold = -0.001; %-0.05,-0.01, -0.005

%Scales and thresholds
thresholds = threshold*(1:3);
factor = 0.2;
scales = linspace(0.8*diameter*factor, 1.3*diameter*factor, 2);

tic
%Load images
V = load_nii(char(iname{i}));
I = uint8(V.img);
I = imresize3(I,[124 124 38]);
I = imresize3(I,[200 200 60]);
%V = load_nii(lname);
%L = V.img;

%for i = 1:60
%   I(:,:,i) = imread(iname,i); 
%end

nDims = size(I,3);
radius = round((diameter/2-1)/2)*2+1;

%Preprocessing
se = strel('disk',diameter+1);
for ii = 1:size(I,3)
    I(:,:,ii) = imtophat(I(:,:,ii),se);
end
I = normalize_mins(I)* 255;
I = imadjustn(uint8(I),[],[],0.85);
I = single(I);

fprintf('Image pre-processing: %d seconds\n',toc)
I_max = max(I(:));
seeds = ones(size(I));

%Detect cells from Hessian eigenvalues
for scale = scales
    tic
    %Sigma for given scale
    sigmas = scale*resolution;
    
    % Apply 3D filter according to sigma
    I_filt = imGaussianFilter(I,round(sigmas*3),sigmas,'symmetric');
    I_filt = normalize_mins(I_filt) * single(I_max);

    [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = imHessian(I_filt, 1.2);
    [l1,l2,l3]=eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
    
    seeds = and(seeds, l1 < thresholds(1));
    seeds = and(seeds, l2 < thresholds(2));
    seeds = and(seeds, l3 < thresholds(3));
    
    fprintf('Image analysis: %d seconds\n',toc)
    %subplot(1,2,1); imagesc(I_filt(:,:,5)); title('first image');
    %subplot(1,2,2); imagesc(seeds(:,:,17)); title('second image');
    %subplot(1,3,1); imagesc(l1(:,:,5))
    %subplot(1,3,2); imagesc(l2(:,:,5))
    %subplot(1,3,3); imagesc(l3(:,:,5))
end

%Closing
tic
seeds = imclose(seeds, true(1, 1, 3));
seeds = imclose(seeds, true(1, 3, 1));
seeds = imclose(seeds, true(3, 1, 1));

%Opening
struc = false(3, 3, 3);
c = (3+1)/2;
struc(:, c, c) = true; struc(c, :, c) = true; struc(c, c, :) = true;
seeds = imopen(seeds, struc);

%Remove small seeds
%seeds = bwareaopen(seeds,7,26);

CC = bwconncomp(seeds,26);

ave_int = median(I(:));

for w = 1:CC.NumObjects
    ind = CC.PixelIdxList{w};
   pixels = I(ind);
   
   if max(pixels) < ave_int*2
       seeds(ind) = 0;
   end
end

CC = bwconncomp(seeds,26);
objects(i) = CC.NumObjects;
end
disp(true_objects)
disp(objects)


figure
imshowpair(seeds(:,:,1),imadjust(uint8(I(:,:,1))))
figure
imshowpair(seeds(:,:,2),imadjust(uint8(I(:,:,2))))
figure
imshowpair(seeds(:,:,3),imadjust(uint8(I(:,:,3))))

L = zeros(size(I));
index = randperm(CC.NumObjects);
for k = 1:CC.NumObjects
   ind = CC.PixelIdxList{k};
   L(ind) = index(k);
end

L2 = label2rgb3d(L,'jet');
L2 = imcomplement(L2);

for k = 1:size(L2,3)
   slice = cat(3,L2(:,:,k,1),  L2(:,:,k,2),  L2(:,:,k,3));
    imwrite(slice,sprintf('slice_%d.tif',k))
end






