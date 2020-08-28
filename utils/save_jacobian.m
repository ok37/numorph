function J = save_jacobian(config, direction)
%--------------------------------------------------------------------------
% Calculate determinant of the jacobian from calculated registration
% parameters to determine regions of volumetric compression/expansion
%--------------------------------------------------------------------------

if nargin<2
    direction = 'atlas_to_img';
end

% Placeholder
% For mac, sometimes need to set PATH and LIBRARY variables explicitly each
% time
%path1 = getenv('PATH');
%path1 = [path1, ':/Users/Oleh/Programs/elastix-5.0.0-mac/bin'];
%setenv('PATH',path1)
%path2 = '/Users/Oleh/Programs/elastix-5.0.0-mac/lib';
%setenv('LD_LIBRARY_PATH',path2)

load(fullfile(config.output_directory,'variables','reg_params.mat'))

J = transformix([],reg_params.(direction),[1,1,1],[],'jac');

save_path = fullfile(config.output_directory,'resampled');
file_name = fullfile(save_path,...
    sprintf('%s_spatialJacobian.nii',config.sample));

niftiwrite(J, file_name)

end