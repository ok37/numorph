function [z_displacement,q,low_flag] = z_align_stitch(path_mov,path_ref,overlap,z_positions,z_window,direction,signalThresh)
%--------------------------------------------------------------------------
% Determine pairwise z displacement between z-stack tiles by taking
% sections within a moving image and performing 2D phase correlation to the
% overlapped region with the reference stack. Average cross correlation for
% the selected slices in the stack are then used to determine the final
% z-displacement.
%--------------------------------------------------------------------------

% Defaults
peaks = 3;
usfac = 1;
min_overlap = 50;

% Adjust signalThresh to 16 bit
if signalThresh<1
    signalThresh = signalThresh*65535;
end

% Check file lengths
nfiles_mov = height(path_mov);    
nfiles_ref = height(path_ref);
if nfiles_mov ~= nfiles_ref
    error("Number of images do not match up")
end

% Adjust z_positions if specified number is too high
z_positions = min(z_positions,nb_imgs_mov-z_window*2);

% Pick reference images in range
z = round(linspace(z_window+1,nb_imgs_ref-z_window,z_positions));
path_ref = path_ref(z,:);

% Determine overlap region
ref_img = imread(path_ref.file{1});
[nrows, ncols] = size(ref_img);

if direction == 1 && ncols*overlap < min_overlap
    overlap = min(1,min_overlap/ncols);
elseif direction == 1 && nrows*overlap < min_overlap
    overlap = min(1,min_overlap/nrows);
end

if direction == 1
    %Horizontal image overlap regions
    overlap_min = {[1,nrows],[1,round(ncols*overlap)]};
    overlap_max = {[1,nrows],[ncols-overlap_min{2}(2)+1,ncols]};
else
    %Vertical image overlap regions
    overlap_min = {[1,round(nrows*overlap)],[1,ncols]};
    overlap_max = {[nrows-overlap_min{1}(2)+1,nrows],[1,ncols]};
end

%Perform registration using phase correlation
cc = zeros(z_positions,z_window*2+1);
f = zeros(1,length(z));
a=1;
for i = 1:length(z)
    % Define range for current z position
    z_sub_range = z(i)-z_window:1:z(i)+z_window;
    % Read ref image
    ref_img = imread(path_ref.file{i});
    % This is dumb - problem with imread loading 3 channel image
    ref_img = ref_img(overlap_max{1}(1):overlap_max{1}(2),overlap_max{2}(1):overlap_max{2}(2),:);
    % Check for multi-channel
    if size(ref_img,3)>1
       ref_img = ref_img(:,:,1); 
    end
    % Measure number of positive pixels
    f(i) = numel(ref_img(ref_img>signalThresh))/numel(ref_img);
    b=1;
    for j = z_sub_range
        % Read corresponding moving image
        mov_img = imread(path_mov.file{j});
        % This is dumb - problem with imread loading 3 channel image
        mov_img = mov_img(overlap_min{1}(1):overlap_min{1}(2),overlap_min{2}(1):overlap_min{2}(2),:);
        % Check for multi-channel
        if size(mov_img,3)>1
            mov_img = mov_img(:,:,1); 
        end
        % Get transform using phase correlation and measure cross
        % correlation
        [~,~,~,~,cc(a,b)] = calculate_phase_correlation(mov_img,ref_img,peaks,usfac);
        b = b+1;
    end
    a = a+1;
end

% Pick highest correlated
cc(isnan(cc))=0;

if sum(cc(:)) > 0
    [~,a1] = max(cc');
    f_num = length(f(f>0.005)); 
else
    f_num = 0;
end
disp(mean(cc))

%If images contain positive pixels, save best transform. Otherwise set as
%NaN and mark flag
if f_num>0
    a1(f<0.001)=NaN;
    a1(all(cc'==-1))=NaN;

    %Highest correlated is the one with 
    z1 = mode(a1);
    q1 = cc(:,z1);
    cc(:,z1) = -1;

    %Pick second highest correlated
    [~,a2] = max(cc');    
    z2 = mode(a2);
    q2 = cc(:,z2);

    %Take z_displacements of both
    z_displacement = z1 - (z_window+1);
    %disp(z_displacement)

    %q is the difference in quality of the 2 best matches
    q = mean(q1-q2);
    low_flag = 0;
else
    %Save tranform
    z_displacement = NaN;
    q = -1;
    low_flag = 1;
    %disp('low')
end
end