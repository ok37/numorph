function visualize_flat_cortex(input, beta_col, p_col, flatview)

%rescale_flag = true;
n_colors = 201;
limits = [0,2];

adj = (limits(2)-limits(1))/(n_colors);

if nargin<2
    beta_col = 'value';
end

if nargin<3
    p_col = [];
    p = [];
end

if nargin<4
    if size(input,1)<18
        flatview = 'seventeen';
    else
        flatview = 'harris';
    end
end

% Load flatview cortex
load(fullfile(fileparts(which('NM_config')),'data','flatviewCortex.mat'),'fc')
if isequal(flatview,'harris')
    fc = fc.harris;
    ftable = readtable(fullfile(fileparts(which('NM_config')),'annotations','harris_cortical_groupings.csv'));
else
    fc = fc.seventeen;
    ftable = readtable(fullfile(fileparts(which('NM_config')),'annotations','cortex_17regions.xls'));
end

% Read cortical groupings
indexes = ftable.index;
    
if istable(input)
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
end

% Get mask and boundaries
boundaries = fc.boundaries;
bw = fc.flatCtx==0 & ~boundaries;

% Round vector to nearest hundreth
v = round(v,2);
v(v>2) = 2;
img = zeros(size(fc.flatCtx));
for i = 1:length(indexes)
    img(fc.flatCtx==indexes(i)) = v(i);
end

lowerThresh = min(1,min(v));
upperThresh = max(1,max(img(:)));


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
title(beta_col)

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
    if isempty(p) || ~sig(i)
        text(x(i),y(i),fc.acronym{i},'FontName','Arial','FontSize',11)
    else
        text(x(i),y(i),strcat(fc.acronym{i},'*'),...
            'FontName','Arial','FontSize',11,'FontWeight','bold')
    end
end

set(gcf,'CurrentAxes',ax1)
h = colorbar;
h.Limits = [limits(1),limits(2)];
h.Ticks = limits(1):diff(limits)/4:limits(2);
h.Position(1) = 0.85;

end
