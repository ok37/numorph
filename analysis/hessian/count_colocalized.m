function df = count_colocalized(path_table,markers, thresholds, I_mask)

if ismepty(I_mask)
    fprintf('%s\t Ignoring image annotations.\n',datetime('now'));
end

df = readtable('/home/ok37/repos/3dunet-centroid/nuclei/centroids2.csv');
df = table2array(df);

%df2 = readtable('/home/ok37/repos/3dunet-centroid/nuclei/centroids_top.csv');
%df2 = table2array(df2);

%df = readtable('top_16_centroids.csv');
%df = table2array(df);
df = df+1;
df = df';

%df = load('WT_centroids.mat');
%df = df.centroids;
%df = load('TOP_centroids.mat');
%df = df.centroids;
%se = strel('disk',15);

%I_annotation = loadtiff(path_to_annotation);

I_annotation = I_mask;
num_images = height(path_table)/length(markers);
I_annotation = imresize3(I_annotation,[size(I_annotation,1),size(I_annotation,2),num_images],...
    'method','nearest');

img_list = unique(df(3,:));

img_list = img_list(:,img_list<1266);

disp(max(img_list))


v = zeros(length(markers),size(df,2));
df = vertcat(df,v);

tempI = imread(path_table.file{1});
num = 0;

for i = 2:length(markers)
   path_sub = path_table(path_table.markers == markers(i),:);
   threshold = thresholds(i-1);
      
   for j = img_list
       I = imread(path_sub.file{j});
       I_annotation_slice = imresize(I_annotation(:,:,j),size(I),'method','nearest');
       
       %if i == 3
        %I = imtophat(I,se);
        %y_adj = linspace(0.9,1.11,size(tempI,1));
        %I = uint16(double(I).*y_adj');
        %x_adj = linspace(0.9,1.11,size(tempI,2));
        %I = uint16(double(I).*x_adj);
       %else
       % I = imtophat(I,se);
       % y_adj = linspace(0.5,2,size(tempI,1));
       % I = uint16(double(I).*y_adj');
       % x_adj = linspace(0.85,1.17,size(tempI,2));
       % I = uint16(double(I).*x_adj);
       %end

       I = imgaussfilt(I,1);
       points_ind = df(3,:) == j;
       index = find(points_ind);
       

       for k = 1:length(index)
           value = I(df(1,index(k)),df(2,index(k)));
           if value > threshold
               df(3+i-1,index(k)) = 1;
		num = num+1;
           end

            structure_index = I_annotation_slice(df(1,index(k)),df(2,index(k)));
            df(end,index(k)) = structure_index;
       end
	fprintf('Counted %d cells at  slice: %d\n',num,j)
   end
end

end
