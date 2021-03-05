function ud = generate_classifications(inputs, patch_path, info_path)
%--------------------------------------------------------------------------
% Generate manual classifications for image patches to train SVM
% classifier. This creates a minimal app for users to quickly view image
% patches and assign classifier values.
%--------------------------------------------------------------------------
% Usage:
% ud = generate_classifications(inputs, patch_path, info_path)
% 
%--------------------------------------------------------------------------
% Inputs:
% inputs: (config structure or string) Sample id or config structure from 
% NM_analyze. If config structure, all other information loaded from this.
% If sample id, searches current folder for sample patches and info.
%
% patch_path: (string) Path to images patches. (optional)
%
% info_path: (string) Path to patch info csv file to display during 
% annotation. (optional)
%
%--------------------------------------------------------------------------
% Outputs:
% ud: Image patch object.
%
%--------------------------------------------------------------------------

% Get path to centroid from config structure
if isstruct(inputs)
    img_path = fullfile(config.output_directory,'classifier',...
        sprintf('%s_patches.tif',config.sample_id));
    img = loadtiff(img_path);
    itable = fullfile(config.output_directory,'classifier',...
        sprintf('%s_patch_info.csv',config.sample_id));
    info = readmatrix(itable);
    sample = config.sample_id;
end

if ischar(inputs) || isstring(inputs)
    sample = inputs;
        
    % Load patch
    if nargin>1
        img = loadtiff(patch_path);
    else
        patch_name = sprintf('%s_patches.tif',sample);
        patch_path = fullfile(pwd,patch_name);
        if isfile(patch_path)
            img = loadtiff(patch_path);
        else
            error("Could not locate %s in current directory",patch_name)
        end
    end
    
    % Load info
    if nargin>2
        info = readmatrix(info_path);
    else
        info_name = sprintf('%s_patch_info.csv',sample);
        info_path = fullfile(pwd,info_name);
        if isfile(info_path)
            info = readmatrix(info_path);
        else
            error("Could not locate %s in current directory",info_name)
        end
    end
end

if size(img,4) == 3
    img = permute(img,[1,2,4,3]);
end

ud.z = 1;
ud.c = 'c';
ud.s = 'all';
ud.sample = sample;
ud.class = zeros(1,size(img,4));
layers = bin_annotation_structures(info(:,5),'layers');

[nrows,~] = size(img);
point = zeros(nrows,nrows);
point(round(nrows/2),round(nrows/2)) = 1;
point = imdilate(point,[0 1 0; 1 1 1; 0 1 0]);

for i = 1:size(img,4)
    a = squeeze(img(:,:,:,i));
    vs.slice{i} = a;
    vs.slice{i} = imoverlay(vs.slice{i},point,[1,1,1]);
    vs.layer(i) = layers(i);
    vs.structure{i} = "cortex";
end

clear img
f = figure;
ud = display_slice(f, vs, ud);
truesize(f,[800,800])

set(f, 'UserData', ud);
set(f, 'KeyPressFcn', @(f,k)hotkeyFcn(f, k, vs)); % Hot keys
display_help_menu

clear ud
end

function hotkeyFcn(f, keydata, vs)

ud = get(f, 'UserData');

switch keydata.Key
    case 'rightarrow'
        ud = advance_slice(f, vs, ud, 'r');
    case 'leftarrow'
        ud = advance_slice(f, vs, ud, 'l');
    case 'c'
        ud.c = 'c';
        ud = display_slice(f, vs, ud);
    case 'r'
        if ud.c == 'r'
            ud.c = 'c';
        else
            ud.c = 'r';
        end
        ud = display_slice(f, vs, ud);
    case 'g'
        if ud.c == 'g'
            ud.c = 'c';
        else
            ud.c = 'g';
        end
        ud = display_slice(f, vs, ud);
    case 'b'
        if ud.c == 'b'
            ud.c = 'c';
        else
            ud.c = 'b';
        end
        ud = display_slice(f, vs, ud);
    case 'f'
        if ~isequal(ud.s,'all')
            return
        end        
        idx = find(ud.class == 0);
        if ~isempty(idx)
            ud.z = idx(1);
            ud = display_slice(f, vs, ud);
        else
            disp("All patches annotated!")
        end
    case 's'
        x = input('Enter cell-type index or ''all'' to release subset: ','s');
        ud.s = x;
        if ~isequal(x,'all')
            ud.s = str2double(x);
        end
        ud = advance_slice(f, vs, ud, 'r');
        figure(f)
    case 'w'
        ttable = array2table(ud.class','VariableNames',{'Type'});
        fname = sprintf('%s_classifications.csv',ud.sample);
        writetable(ttable,fname)
    case 'l'
        fname = sprintf('%s_classifications.csv',ud.sample);
        if isfile(fname)
            ttable = readtable(fname);
            ud.class = ttable{:,1}';
            ud = display_slice(f, vs, ud);
        else
            fprintf("Could not locate previous classification table in current "+...
                "directory\n")
        end
    case 'h'
        display_help_menu
    case '1'
        ud.class(ud.z) = 1;
    case '2'
        ud.class(ud.z) = 2;
    case '3'
        ud.class(ud.z) = 3;
    case '4'
        ud.class(ud.z) = 4;
    case '5'
        ud.class(ud.z) = 5;
    case '6'
        ud.class(ud.z) = 6;
    case '7'
        ud.class(ud.z) = 7;
    case '8'
        ud.class(ud.z) = 8;
    case '9'
        ud.class(ud.z) = 9;
end

if ismember(keydata.Key,['1','2','3','4','5','6','7','8','9'])
    ud = advance_slice(f, vs, ud, 'r');
elseif ismember(keydata.Key,['c','r','g','b'])
    ud = display_slice(f, vs, ud);
end

set(f, 'UserData', ud);

end


function ud = display_slice(f, vs, ud)

img = vs.slice{ud.z};
switch ud.c
    case 'r'
    img(:,:,2:3) = 0;
    case 'g'
    img(:,:,[1 3]) = 0;
    case 'b'
    img(:,:,1:2) = 0;
end

ud.img = img;

if isempty(f.Children)
    imshow(ud.img);
else
    imshow(ud.img, 'Parent',f.Children)
end

text(2,4,string(ud.z),'FontSize',30, 'Color',[1,1,1])
text(4,size(vs.slice{1},1)-6,string(vs.layer(ud.z)),'FontSize',30, 'Color',[1,1,1])

text(size(vs.slice{1},1)-6,4,string(ud.class(ud.z)),...
    'FontSize',30, 'Color',[1,1,1])
end


function [ud, vs] = advance_slice(f, vs, ud, direction)

if isequal(ud.s,'all')
    idx = ones(1,length(ud.class),'logical');
elseif ~any(ud.class==ud.s)
    fprintf("No cells with the selected type detected \n")
    return
else
    idx = ud.class==ud.s;
end
    
if isequal(direction,'r')
    if ud.z == length(vs.slice)
        ud.z = 1;
    else
        ud.z = ud.z+1;
        while ~idx(ud.z)
            ud.z = ud.z+1;
            if ud.z == length(vs.slice)+1
                ud.z = 1;            
            end
        end
    end
else
    
    if ud.z == 1
        ud.z = length(vs);
    else
        ud.z = ud.z-1;
        while ~idx(ud.z)
            ud.z = ud.z-1;
            if ud.z < 1
                ud.z = length(vs.slice);
            end
        end
    end
end

ud = display_slice(f, vs, ud);
end

function display_help_menu
msg = "Manual annotation key strokes:\n"+...
    "  l/r arrow: navigate patch number\n"+...
    "  1-9:       assign class annotation\n"+...
    "  r:         display only red (1st) channel\n"+...
    "  g:         display only green (2nd) channel\n"+...
    "  b          display only blue (3rd) channel\n"+...
    "  c:         display composite (all channels)\n"+...
    "  h:         display help menu\n"+...
    "  s:         subset patches of a given annotation\n"+...
    "  f:         find nearest unannotated patch\n"+...
    "  w:         save results to csv file\n"+...
    "  l:         load previous results from csv file\n";
fprintf(msg)
end



