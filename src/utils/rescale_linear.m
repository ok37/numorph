function v_out = rescale_linear(v,img_in,res_in,res_out)
% Rescale linear indexes from 1 resolution to another

[y,x,z] = ind2sub(size(img_in),v);
res_adj = res_in./res_out;
s = [y,x,z].*res_adj;

new_dims = round(size(img_in).*res_adj);
v_out = sub2ind(new_dims,s(1),s(2),s(3));

end