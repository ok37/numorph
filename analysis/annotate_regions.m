%Function to open a registered image containing annotations for specific
%regions and creating a mask of these regions
function annotate_regions(path_to_atlas, path_to_id, save_flag)

A = nrrdread(strcat(pwd,'/supplementary files/annotation_25.nrrd'));
%I = loadtiff(strcat(paths,'/resampled/C3_ToPro__resampled.tif'));

id = csvread(strcat(paths,'/annotation/cortex.csv'));

%R = zeros(size(A));

%for i = 1:size(id,1)
%    R(A==id(i)) = 1;
%end

%options.color     = false;
%options.compress  = 'no';
%options.message   = false;
%options.append    = true;
%options.overwrite = true;

%I2 = single(R).*single(I);

%for i = 1:size(I2,3)
%    saveastiff(uint16(I2), 'cortex.tif', options);
%end



R = zeros(size(A));
G = zeros(size(A));
B = zeros(size(A));

for i = 1:size(id,1)
    R(A==id(i)) = id(i,2);
    G(A==id(i)) = id(i,3);
    B(A==id(i)) = id(i,4);
end

options.color     = true;
options.compress  = 'no';
options.message   = false;
options.append    = true;
options.overwrite = true;

for i = 1:size(R,3)
    I2(:,:,1) = R(:,:,i);
    I2(:,:,2) = G(:,:,i);
    I2(:,:,3) = B(:,:,i);

    saveastiff(uint8(I2), 'R_cortex.tif', options);
    %saveastiff(uint8(G), 'G_cortex.tif', options);
    %saveastiff(uint8(B), 'B_cortex.tif', options);

end
end