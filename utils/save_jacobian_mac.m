function J = save_jacobian(reg_params_path, outputDir, sample_n)
% Calculate determinant of the jacobian from calculated registration
% parameters to calculate regions of volumetric compression/expansion

% Placeholder
% Just for Mac
path1 = getenv('PATH');
path1 = [path1, ':/Users/Oleh/Programs/elastix-5.0.0-mac/bin'];
setenv('PATH',path1)
path2 = '/Users/Oleh/Programs/elastix-5.0.0-mac/lib';
setenv('LD_LIBRARY_PATH',path2)

load(reg_params_path,'reg_params')

J = transformix([],reg_params,[1,1,1],[],'jac');

file_name = fullfile(outputDir,...
    sprintf('%s_spatialJacobian.nii',sample_n));

niftiwrite(J, file_name)

end