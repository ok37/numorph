function visualize_flat_cortex(input, key, flatview)
%--------------------------------------------------------------------------
% Visualize flattened cortex at 10um resolution. 
%--------------------------------------------------------------------------
% Usage:
% visualize_flat_cortex(input, key, flatview)
%
%--------------------------------------------------------------------------
% Inputs:
% input:    1. Flatmap structure containg voxelized cortex data for heatmap
%              representation. 
%           2. Table containing fold change and p values by cortical 
%              annotation index.
%           3. 2D image at the same size as 10um flatmap (1360x1360)
%
% key: (1x3 integer) First index indicates cell class index. Second index
% indicates statistic (1=fold change, 2=p-value, 3=p.adj). Third index
% indicates thresholding signficant regions (1=true, 0=false).
% 
% flatview: ('harris' or 'seventeen') Bin cortical structures. 
%
%--------------------------------------------------------------------------

%rescale_flag = true;
n_colors = 201;
limits = [-100,100];
adj = (limits(2)-limits(1))/(n_colors);

if nargin<2
    key = [1,1,1];
elseif length(key)<3
    key = [key,0];
end

if nargin<3
    flatview = 'harris';
end

if key(2) == 1
    limits = [-100,100];
elseif key(2) > 1
    if key(3) == 0
        limits = [0,5];
    else
        limits = [1.301,5];
    end
end

%figure;imagesc(results2)

% Load flatview cortex
fc = load(fullfile(fileparts(which('NM_config')),'data','annotation_data','flatviewCortex.mat'),flatview);
if isequal(flatview,'harris')
    fc = fc.harris;
    ftable = readtable(fullfile(fileparts(which('NM_config')),'annotations','custom_annotations','harris_cortical_groupings.xls'));
else
    fc = fc.seventeen;
    ftable = readtable(fullfile(fileparts(which('NM_config')),'annotations','custom_annotations','cortex_17regions.xls'));
end

% Read cortical groupings
indexes = ftable.index;
    
if isstruct(input)
    img = squeeze(input(key(1)).data(:,:,key(2)));
    
    if key(3) == 1
        sig_bw = squeeze(input(key(1)).data(:,:,3));
        sig_bw = imgaussfilt(sig_bw,15);
        sig_bw = -log10(sig_bw)<1.301;
    else
        sig_bw = [];
    end
    
    if key(2) == 1
        img = img*100-100;
    else
        img = -log10(img);
    end
    img = imgaussfilt(img,15);
    img(sig_bw) = 0;
    t_disp = input(key(1)).marker;
    
elseif istable(input)
    input = input(ismember(input.index,indexes),:);
    v = input.(beta_col);
    if ~isempty(p_col)
        p = input.(p_col);
        %v(p>0.05) = 1;
        sig = p<0.05;
    end
    i_index = input.index;
    [~, i] = sort(i_index);
    v = v(i);
    if ~isempty(p_col)
        sig = sig(i);
    end
    
    % Round vector to nearest hundreth
    v = round(v,2);
    v(v>2) = 2;
    img = zeros(size(fc.flatCtx));
    for i = 1:length(indexes)
        img(fc.flatCtx==indexes(i)) = v(i);
    end
elseif isnumeric(input)
    img = input;
    img(isnan(img)) = 0;
    
    
    
    img = imgaussfilt(img,10,'FilterSize',13);    
    p = [];
    t_disp = "Flatview Image";
    
end

% Get mask and boundaries
boundaries = fc.boundaries;
bw = fc.flatCtx==0 & ~boundaries;
%lowerThresh = min(1,min(v));
%upperThresh = max(1,max(img(:)));

if ishandle(1)
    close(gcf)
end
ax1 = axes;
imagesc(img);

%%img(bw) = limits(2)+adj;
%%img(fc.boundaries) = limits(2)+adj*2;

% Load colors
colors = brewermap(n_colors,'*RdBu');

%x1 = linspace(limits(1),limits(2),n_colors);
%idx = x1>=lowerThresh & x1<=upperThresh;
%colors = colors(idx,:);

%colors = [colors;[1,1,1];[0,0,0]];
colormap(colors)
caxis([limits(1) limits(2)+adj*2])
%h = colorbar;
%h.Limits = [limits(1),limits(2)];
%h.Ticks = limits(1):diff(limits)/4:limits(2);
set(ax1,'XTick',[], 'YTick', [])
daspect([1,1,1])

view(2)
ax2 = axes;
a2 = image(bw);
a2.AlphaData = bw;
ax2.Colormap = [[0,0,0];[1,1,1]];
ax3 = axes;
a3 = image(fc.boundaries);
a3.AlphaData = fc.boundaries;
ax3.Colormap = [[1,1,1];[0,0,0]];
linkaxes([ax1,ax2,ax3])
ax1.Visible = 'off';
ax1.XTick = [];
ax1.YTick = [];
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
ax3.Visible = 'off';
ax3.XTick = [];
ax3.YTick = [];

P = get(ax1,'Position');
XLIM = get(ax1,'XLim');
YLIM = get(ax1,'YLim');
PA = get(ax1,'PlotBoxAspectRatio');

set(ax2,'Position',P,'XLim',XLIM,'YLim',YLIM,'PlotBoxAspectRatio',PA)
set(ax3,'Position',P,'XLim',XLIM,'YLim',YLIM,'PlotBoxAspectRatio',PA)

% Add structure names
x = fc.x;
y = fc.y;
for i = 1:length(indexes)
    if isstruct(input) || isempty(p) || ~sig(i)
        text(x(i),y(i),fc.acronym{i},'FontName','Arial','FontSize',11)
    else
        text(x(i),y(i),strcat(fc.acronym{i},'*'),...
            'FontName','Arial','FontSize',11,'FontWeight','bold')
    end
end
text(0,-30,t_disp,'FontName','Arial','FontSize',11,'FontWeight','bold')

set(gcf,'CurrentAxes',ax1)
h = colorbar;
h.Limits = [limits(1),limits(2)];
h.Ticks = limits(1):diff(limits)/4:limits(2);
h.Position(1) = 0.85;


end
