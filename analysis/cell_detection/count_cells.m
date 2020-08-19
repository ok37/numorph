% Check for stitched images
if exist(char(strcat(output_directory,'/stitched')),'dir') == 7
    path_stitch = dir(char(strcat(output_directory,'/stitched')));

    markercell = cellstr(repmat("",1,length(path_stitch)));
    [path_stitch.markers] = markercell{:};
    [path_stitch.z] = markercell{:};

    %Need to update these
    path_stitch = path_stitch(contains({path_stitch.name},'.tif'));
    
    name = strsplit([path_stitch.name],'_');
    z = str2double(name(3:3:end))';

    path_table_main = struct2table(path_stitch);

    path_table_main.markers(contains(path_table_main.name,"C1"))={'Ctip2'};
    path_table_main.markers(contains(path_table_main.name,"C2"))={'Cux1'};
    path_table_main.markers(contains(path_table_main.name,"C3"))={'ToPro'};
    
    path_table_main.z = z;
    path_table_main = sortrows(path_table_main,'z');
end



path_table_main = path_table_main(path_table_main.markers == "ToPro",:);
data = [0;0;0];

chunks = 1:5:height(path_table_main);

for index = 1:length(chunks)

if chunks(index) ~= chunks(end)
    z_pos = chunks(index):chunks(index+1)-1;
else
    z_pos = chunks(index):height(path_table_main);
end

fprintf(strcat(char(datetime('now')),'\t Segmenting nuclei\n'))
parfor i = z_pos
    name = char(strcat(output_directory,'/stitched/',path_table_main.name(i)));
    
    I=imread(name);
    blobs = count_cells_worker(I);
    blobs = vertcat(blobs,repmat(i,[1,length(blobs)]));
    
    data = horzcat(data,blobs);
end

data(:,1) = [];
point = zeros(3,1);
data = single(data);

fprintf(strcat(char(datetime('now')),'\t Finding Centroids\n'))
a = 1;
tic
for i = 1:length(data)
   z = data(3,i);
   
   if z ~= 0
    x = data(1,i);
    y = data(2,i);

    x_min = x-2;
    x_max = x+2;
   
    y_min = y-2;
    y_max = y+2;
   
    z_min = z-4;
    z_max = z+4;

    ind = data(1,:)>x_min & data(1,:)<x_max & data(2,:)>y_min &...
        data(2,:)<y_max & data(3,:)>z_min & data(3,:)<z_max;
    
    xyz = data(:,ind);
    xyz(:,xyz(3,:)==0)=[];
   
   if size(xyz,2)>3 
   %for j = 1:length(xyz)-1
   %    dist = pdist([xyz(1:2,1),xyz(1:2,j+1)],'euclidean');
   %    if dist<2
   %        xyz2 = horzcat(xyz2,xyz(:,j+1));
   %    end
   %end
   point = horzcat(point,round(mean(xyz,2)));
   end
   data(3,ind) = 0;

    disp(i)
   end

end


end





point(:,1) = [];

%%

a = 1;



function blobs = count_cells_worker(I)

se = strel('disk',15);
I_filt = imtophat(I,se);

%Generate binary
B = imbinarize(I_filt,0.002);

B = bwpropfilt(B,'Area',[15,200]);
B = imfill(B,'holes');

D = -bwdist(~B);
D(~B) = -Inf;
L = watershed(D);
L(~B) = 0;

BW= imbinarize(L, 0);
BW = bwareaopen(BW,5);

c = regionprops(BW,'Centroid','Area');

blobs = reshape([c.Centroid],[2,size(c)]);
%blobs = vertcat(blobs,[c.Area]);
end