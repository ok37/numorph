function [r2,z_displacement] = zAlign(path_mov,path_ref,step,max_range,min_range)
%This function searches for an intial z displacement by taking sections
%within a moving image and performing 2D registration to a subvolume with
%the reference stack. The registration accuracy is then determined by
%structural similarity index and the z displacement with the highest
%average accuracy from all slices is used as the final z displacement

nfiles_mov = length(path_mov);    
nfiles_ref = length(path_ref);

%Check file lengths
if nfiles_mov ~= nfiles_ref
    error("Number of images do not match up")
end

%Pick reference images in range
z = [step:step:nfiles_ref];
image_ref = path_ref(z);


mid = round(length(z)/2);
ref_image = imread(image_ref(mid).name);
a = 1;
%Read moving images to register in range
for i = z(mid)-(max_range/2):z(mid)+(max_range/2)
        mov_filename = path_mov(i).name;
        mov_image = imread(mov_filename);
        
        tformEstimate = imregcorr(mov_image,ref_image,'translation'); 
        Rfixed = imref2d(size(ref_image));
        Mfixed = imwarp(mov_image,tformEstimate,'OutputView',Rfixed);
        r1(a) = ssim(Mfixed,ref_image)
        a = a+1;
end


for j = 1:length(r1)-(min_range)
    r_window(j) = mean(r1([j:j+min_range-1])); 
end

[~,ind] = max(r_window);
range = -max_range/2:1:max_range/2;
window = ind:1:ind+min_range-1;
range = range(window);

for i = z
    %Read reference image
    ref_filename = path_ref(i).name;
    ref_image = imread(ref_filename);
    
    a = 1;
    %Read moving images to register in range
    for j = i + range
        mov_filename = path_mov(j).name;
        mov_image = imread(mov_filename);
        
        tformEstimate = imregcorr(mov_image,ref_image,'translation'); 
        Rfixed = imref2d(size(ref_image));
        Mfixed = imwarp(mov_image,tformEstimate,'OutputView',Rfixed);
        r2(i/step,a) = ssim(Mfixed,ref_image)
        a = a+1;
    end
end

[~,z_displacement] = max(mean(r2));
z_displacement = range(z_displacement);
end