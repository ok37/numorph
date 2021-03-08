function centroids = reannotate_centroids(centroids,config)
% Load mask and reannotate centroids 

I_mask = load_var(config,strcat(config.sample_id,'_mask'),'I_mask');
z_pos = unique(centroids(:,3));
s = config.resolution{1}/config.resample_resolution;

% Resize mask
[mrows,mcols,mslices] = size(I_mask);

% Adjust for python base 0 indexing
centroids(:,1:3) = centroids(:,1:3);

yxz = round(centroids(:,1:3).*s(1:3));
yxz(yxz == 0) = 1;
yxz(:,1) = min(yxz(:,1),mrows);
yxz(:,2) = min(yxz(:,2),mcols);
yxz(:,3) = min(yxz(:,3),mslices);

yxz = sub2ind([mrows,mcols,mslices],yxz(:,1),yxz(:,2),yxz(:,3));
a = I_mask(yxz);
centroids(:,4) = a;

end






