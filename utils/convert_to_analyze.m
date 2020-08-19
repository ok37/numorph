%Allen Nissl reference
fname1 = 'fixed_nissl.tif';
%Downsampled nuclear stained image
fname2 = '180617_WT10R-ToPro-0-08x-0-017NA-6z_09-19-1.tif';
%Get info and number of images in the stack
info = imfinfo(fname1);
info2 = imfinfo(fname2);
num_images = numel(info);
num_images2 = numel(info2);

%Read Nissl files and convert to 8-bit
I_mov = zeros(info(1).Height, info(1).Width, num_images);
for k = 1:num_images
    I_mov(:,:,k) = imread(fname1, k);
end
I_mov = uint8(I_mov*255/65535);

%Read acquired files and convert to 8-bit
I_ref = uint16(zeros(info2(1).Height, info2(1).Width, num_images2));
for k = 1:num_images2
    I_ref(:,:,k) = imread(fname2, k);
end

%Adjust intensity for nuclear image
I_ref = imadjustn(I_ref,[],[],1.2);
I_ref = im2uint8(imadjustn(I_ref,stretchlim(I_ref(:),[0.001 .99])));

%Resize Nissl file to match nuclear image
ref_size = size(I_ref);
%I_mov = imresize3(I_mov,ref_size);
I_mov = imadjustn(I_mov,stretchlim(I_mov(:),[0.01 .975]),[],1.1);

%Specify volumes
fixedVolume = I_ref;
movingVolume = I_mov;

centerFixed = round(size(fixedVolume)/2);
centerMoving = round(size(movingVolume)/2);

Rfixed  = imref3d(size(fixedVolume),25,25,25);
Rmoving = imref3d(size(movingVolume),25,25,25);

[optimizer,metric] = imregconfig('multimodal');

tform = imregtform(movingVolume,Rmoving,fixedVolume,Rfixed, 'affine',....
    optimizer, metric,'PyramidLevels',4);
movingRegisteredVolume2 = imwarp(movingVolume,Rfixed,tform,'OutputView',Rfixed);

figure 
imshowpair(fixedVolume(:,:,centerFixed(3)), movingRegisteredVolume2(:,:,centerFixed(3)));


fnamec1= 'C1-annotation_25_left_colorR.tif';
fnamec2= 'C2-annotation_25_left_colorG.tif';
fnamec3= 'C3-annotation_25_left_colorB.tif';

for k = 1:num_images
    C1(:,:,k) = imread(fnamec1, k);
    C2(:,:,k) = imread(fnamec2, k);
    C3(:,:,k) = imread(fnamec3, k);
end

C1 = imresize3(C1,size(fixedVolume),'nearest'); 
C2 = imresize3(C2,size(fixedVolume),'nearest'); 
C3 = imresize3(C3,size(fixedVolume),'nearest'); 

C1 = imwarp(C1,Rfixed,tform,'OutputView',Rfixed,'interp','nearest');
C2 = imwarp(C2,Rfixed,tform,'OutputView',Rfixed,'interp','nearest');
C3 = imwarp(C3,Rfixed,tform,'OutputView',Rfixed,'interp','nearest');


writeanalyze(movingRegisteredVolume2,'mov',[25,25,25])
writeanalyze(fixedVolume,'ref',[25,25,25])

writeanalyze(C1,'C1',[25,25,25])
writeanalyze(C2,'C2',[25,25,25])
writeanalyze(C3,'C3',[25,25,25])
