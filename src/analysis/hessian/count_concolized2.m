% Check for stitched images
output_directory = "/home/ok37/Stitching/WT1";
load('centroids.mat')

options.color     = true;
options.compress  = 'no';
options.message   = false;
options.append    = true;
options.overwrite = true;


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
centroid = centroids{1};

for i = 2:length(centroids)
    centroid = horzcat(centroid, centroids{i});
end

path_table_1 = path_table_main(path_table_main.markers == "Ctip2",:);
path_table_2 = path_table_main(path_table_main.markers == "Cux1",:);
path_table_3 = path_table_main(path_table_main.markers == "ToPro",:);

M = loadtiff(char(strcat(output_directory,'/annotation/annotation_file.tiff')));
id = csvread(char(strcat(output_directory,'/annotation/cortex.csv')));

M = uint32(M);
M = imresize3(M,[size(M,1),size(M,2),height(path_table_2)]);

path_table_2(634,:) = [];
path_table_3(634,:) = [];
M(:,:,634) = [];

v = path_table_1.z;

counts = zeros(1,5);

a = 1;
for i = 635:864%1:length(v)
   i = v(i); 
   disp(i)
   cen_z = centroid(:,centroid(3,:) == i);
   
   
   
   name1 = char(strcat(output_directory,'/stitched/',path_table_1.name(i)));
   name2 = char(strcat(output_directory,'/stitched/',path_table_2.name(i)));
   
   R = imread(name1);

   R2 = zeros(size(R));
   G2 = zeros(size(R));
   B2 = zeros(size(R));


 if ~isempty(cen_z)
   G = imread(name2);

   M2 = imresize(M(:,:,i),[size(R,1),size(R,2)]);
   
   
   for j = 1:size(cen_z,2)
      
      y = cen_z(1,j);
      x = cen_z(2,j);
      
      R_int = R(x,y);
      G_int = G(x,y);
      
      counts(a,1) = x;
      counts(a,2) = y;
      counts(a,3) = i;
      counts(a,4) = 0;
      counts(a,5) = M2(x,y);
      if R_int >900    
        counts(a,4) = 1;
        G2(x,y) = 1;
        B2(x,y) = 1;
      end

      if G_int >15 && counts(a,4) == 0
          counts(a,4) = 2;
          R2(x,y) = 1;
          G2(x,y) = 1;
      elseif G_int > 20 
          counts(a,4) = 3;
          R2(x,y) = 1;
      end
            a = a+1;
   end
   end
   
   %I_save = cat(3, R2*255, G2*255, B2*255);
   %I_save = imresize(I_save,0.25);
   
   %name_save = char(strcat(output_directory,'/counts/Annotated_',string(i),'.png'));
   %saveastiff(uint8(I_save),name_save,options)
   %imwrite(I_save,name_save)
   
   
end


csvwrite(char(strcat(output_directory,'/cell_counts4.csv')),counts)






