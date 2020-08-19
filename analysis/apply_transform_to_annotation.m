function trans_annot = apply_transform_to_annotation(reg_params,home_path)
% Apply transform to annotation
annotation_path = fullfile(home_path,'supplementary_files','annotation_25-left.nii');

annot_img = niftiread(annotation_path);
annot_img = imrotate(annot_img,90);
annot_img = flipdim(annot_img,1);

for i = 1:length(reg_params.TransformParameters)
    reg_params.TransformParameters{i}.FinalBSplineInterpolationOrder = 0;
end

trans_annot = transformix(annot_img,reg_params,[1 1 1], []);
end