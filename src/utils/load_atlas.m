function img = load_atlas(atlas_type,hemisphere,res_out,or_out)
% Load atlas or annotations and resize to specific resolution

if nargin<1
    atlas_type = 'nissl';
end

if nargin<2
    hemisphere = 'left';
end

if nargin<3
    res_out = 25;
end

if nargin<4
    or_out = 'ail';
end

home_path = fileparts(which('NM_config'));
switch atlas_type
    case 'nissl'
        img = read_img(fullfile(home_path,'data','atlas','ara_nissl_25.nii'));
        img = get_hemisphere(img,hemisphere,'atlas');
        img = permute_orientation(img,'ail',or_out);
        img = imresize3(img,25/res_out);   
        
    case 'average'
        img = read_img(fullfile(home_path,'data','atlas','average_template_25.nii'));
        img = get_hemisphere(img,hemisphere,'atlas');
        img = permute_orientation(img,'ail',or_out);
        img = imresize3(img,25/res_out);  
        
    case 'annotation'
        img = load(fullfile(home_path,'data','annotation_data','annotationData.mat'),'annotationVolume');
        img = img.annotationVolume;
        img = imresize3(img,10/res_out,'Method','nearest');  
        img = get_hemisphere(img,hemisphere,'annotation');
        img = permute_orientation(img,'sra',or_out);
        
    otherwise
        error("Unrecognized atlas type. Options are 'nissl','average', or 'annotation'")
end

end


function img = get_hemisphere(img,hemisphere,atlas_type)

if isequal(atlas_type,'atlas')
    if isequal(hemisphere,'right')
        img = flip(img,3);
    elseif isequal(hemisphere,'both')
        img2 = flip(img,3);
        img = cat(3,img,img2);
    end
else
    if isequal(hemisphere,'left')
        img = flip(img,2);
    elseif isequal(hemisphere,'both')
        img2 = flip(img,2);
        img = cat(2,img,img2);
    end
end

end