function evaluate_classification(config)

path_centroids = fullfile(config.output_directory, sprintf('%s_centroids.csv',config.sample_name));
centroids = readmatrix(path_centroids);

itable = readtable(fullfile(config.output_directory,sprintf('%s_patch_info.csv',config.sample_name)));
ctable = readtable(fullfile(config.output_directory,'classifier',sprintf('%s_classifications.csv',config.sample_name)));

cen_idx = table2array(itable(:,1));
ct = ctable.Type;

cen_sub = centroids(cen_idx,:);

gm = cen_sub(:,8);
gm(gm==0) = 1;
gm(gm>3) = 4;

svm = cen_sub(:,9);
svm(svm==0) = 1;
svm(svm>3) = 4;

ct(ct>3) = 4;

low_thresh = 0.5;
markers = 1:3;
k_idx = zeros(size(centroids,1),1);
for i = 1:length(markers)-1
   idx = i + 5;
   thresh = prctile(centroids(:,idx),low_thresh*100);
   k_idx = k_idx | centroids(:,idx)>thresh;
end
k_sum = sum(k_idx);
nk_sum = sum(~k_idx);

gm = sum(ct==gm)/1000;
svm = sum(ct==svm)/1000;

gm1 = ((gm*k_sum)+nk_sum)/(length(k_idx));
svm1 = ((svm*k_sum)+nk_sum)/(length(k_idx));

disp(gm1)
disp(svm1)



end