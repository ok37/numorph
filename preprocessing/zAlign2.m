function [z_displacement,ave_score] = zAlign2(path_mov,path_ref,z_positions,z_window,lowerThresh,z_initial)
%Set number of peaks and subixel value
peaks = 1;
usfac = 1;

%Read number of reference/moving images
nb_imgs_mov = height(path_mov);    
nb_imgs_ref = height(path_ref);

%Check file lengths
if nb_imgs_mov ~= nb_imgs_ref
    error("Number of images do not match up")
end

%Pick reference images in range
z = linspace(0,nb_imgs_ref,z_positions+2);
z = round(z(2:end-1));
path_ref = path_ref(z,:);

%Initial matrix of cross correlation values and perform registration unsing
%cross correlation
cc = zeros(z_positions,z_window*2+1);
max_signal = zeros(1,length(z));
a=1;
for i = 1:length(z)
    %Read reference image and specify Z range to look based on user-defined
    %window
    ref_img = imread(path_ref.file{i});
    z_range = z(i)-z_window+z_initial:1:z(i)+z_window+z_initial;
    signal = zeros(1,length(z_range));
    b=1;
    for j = z_range
        %Read moving image
        mov_img = imread(path_mov.file{j});
        
        %Calculate number of bright pixels
        signal(b)= numel(mov_img(mov_img>lowerThresh))/numel(mov_img);
        
        if signal(b)<0.01
            cc(a,b) = 0;
        else
        %Do simple crop in case image sizes don't match up
        if any(size(ref_img) ~= size(mov_img))
            [nrows_ref, ncols_ref] = size(ref_img);
            [nrows_mov, ncols_mov] = size(mov_img);
            
            %Just take smallest dimensions from either images
            y_dim = min(nrows_ref,nrows_mov);
            x_dim = min(ncols_ref,ncols_mov);

            ref_img = ref_img(1:y_dim, 1:x_dim);
            mov_img = mov_img(1:y_dim, 1:x_dim);
        end    
        
        %If image singal is low, apply filter to improve phase correlation
        %if signal(b)>0.001 && signal(b)<0.05
            %mov_img = imgaussfilt(mov_img,'FilterSize',5,'FilterDomain','spatial');
            %mov_img = imtophat(mov_img,strel('disk',20)); 
        %end
        
        %Calculate phase correlation
        [pc_img,ref_img2,~,type,cc(a,b)] = calculate_phase_correlation(mov_img,ref_img,peaks,usfac);
        
        %If pc didn't work, set cc to 0
        if type == 0
            cc(a,b) = 0;
        end
        %Display results
        %disp(cc)
        %if i == 5
        %    imshowpair(ref_img2*30,pc_img*250)
        %    adfg = 1;
        %end
        
        end
        b = b+1;
    end
    max_signal(i) = max(signal);
    a = a+1;
end

%Generate score by multiplying intensity value by cross-correlation. This
%will lower effect caused by noise in images with no features
[val, pos] = max(cc'.*max_signal);
idx = unique(pos);
score = zeros(1,length(idx));

%Sum scores for all unqiue positions
for i = 1:length(idx)
    score(i) = sum(val(pos==idx(i)));    
end

%Max score determines final z displacement
[~,k] = max(score);

r(1,:) = idx-(z_window+1)+z_initial;
r(2,:) = score;
disp(r)

%Adjust final z displacement based on initial z position
z_displacement = idx(k)-(z_window+1)+z_initial;
ave_score = mean(val);
end
