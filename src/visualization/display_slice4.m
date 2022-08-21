function [p1,p2] = display_slice4(vs,key)
% User defined keys
% z: z_position
% marker: which marker
% category: counts, volume, density
% stat: which stat
% key: [z_positition, marker, stat];

if key{1}(3) == 1 && key{1}(2) ~= 0
    key{1}(2) = 1;
end
if key{2}(3) == 1 && key{1}(2) ~= 0
    key{2}(2) = 1;
end
dp = set_display_parameters(vs, key);

% Initialize figure
f = figure;
f.Position = [5,600,855,600]; %AR = 1.425

% Plot the left hand side
[p1,cmin,cmax] = display_plot(dp,vs,f);
dp.cmin(1) = cmin; dp.cmax(1) = cmax;
p1 = sanitize_plot(p1,1);

% Plot the right hand side
dp.pos = 2;
[p2,cmin,cmax] = display_plot(dp,vs,f);
dp.cmin(2) = cmin; dp.cmax(2) = cmax;
p2 = sanitize_plot(p2,2);

% Link the axes
linkaxes([p1 p2],'xy')

% Set colorbar
dp.pos = 1;
p1 = add_colorbar(p1,vs,dp);
dp.pos = 2;
p2 = add_colorbar(p2,vs,dp);

% Set text
% Add structure info top left
dp.stext = annotation('textbox', [0 0.95 0.4 0.05], ...
    'EdgeColor', 'none', 'Color', 'k', 'FontSize',12, 'FontWeight','bold');

% Add position info top right
dp.ptext = annotation('textbox', [0.6 0.95 0.4 0.05], ...
    'EdgeColor', 'none', 'Color', 'k','FontSize',12,'HorizontalAlignment','right');

% Add p1 stats to the top left 
dp.p1text = annotation('textbox', [0.0 0.9 0.4 0.05], ...
    'EdgeColor', 'none', 'Color', 'k','FontSize',12,'HorizontalAlignment','left');

% Add p2 stats to the top right
dp.p2text = annotation('textbox', [0.6 0.9 0.4 0.05], ...
    'EdgeColor', 'none', 'Color', 'k','FontSize',12,'HorizontalAlignment','right');

set(dp.ptext, 'String', sprintf('Slice:%d/%d \t %.2f AP', dp.z,size(vs.av,3),vs.ap(dp.z)));

end

function [p,cmin,cmax] = display_plot(dp,vs,f)
% Initial plot display

if isequal(dp.init,'true')
    p = subplot(1,2,dp.pos);
else
    child = get(f,'Children');
    if dp.pos == 1
        p = child(4);
    else
        p = child(2);
    end
end

[slice,colors,title] = retrieve_slice(dp,vs);
h1 = imagesc(single(slice));

% Set caxis mappings
if dp.marker(dp.pos) == 0
    cmin = 0; cmax = length(colors)+1;
    h1.CDataMapping = 'direct';
    caxis([cmin cmax])
else
    [~,~,~,cmin,cmax] = retrieve_data(vs,dp);
    h1.CDataMapping = 'scaled';    
    caxis([cmin cmax])
end

% Add title
p.Title.String = title;
p.Title.FontSize = 14;

% Overlay boundaries
hold on

% Add alpha to annotation volume
a = ~ismember(single(vs.av(:,:,dp.z)),single(vs.indexes));
h3 = image(a);
h1.CDataMapping = 'scaled'; 

if dp.marker(dp.pos) == 0
    set(h3,'AlphaData',a*dp.alpha_val)
else
    set(h3,'AlphaData',a)
end

h3.Visible = 'off';
%if dp.marker(dp.pos) ~= 0
%    h3.Visible  = 'on';
%elseif ~isequal(dp.alpha,'true')    
%    h3.Visible  = 'off';
%end

boundaries = vs.boundaries(:,:,dp.z)*-1;
h2 = image(boundaries);
boundaries = boundaries<0;
set(h2,'AlphaData',boundaries)
hold off

colormap(p,[[0.25,0.25,0.25];colors;])
end

function [title,name,values,cmin,cmax] = retrieve_data(vs,dp)

% Get axis
a = dp.pos;

if isequal(dp.category(a),1)
    % Do volume
    values = vs.markers(dp.marker(a)).Volume.values(:,dp.stat(a));
    title = vs.markers(dp.marker(a)).Volume.title(dp.stat(a));
    name = vs.markers(dp.marker(a)).Volume.name;
    cmin = vs.markers(dp.marker(a)).Volume.cmin(:,dp.stat(a));
    cmax = vs.markers(dp.marker(a)).Volume.cmax(:,dp.stat(a));
elseif isequal(dp.category(a),2)
    % Do counts
    values = vs.markers(dp.marker(a)).counts.values(:,dp.stat(a));
    title = vs.markers(dp.marker(a)).counts.title(dp.stat(a));
    name = vs.markers(dp.marker(a)).counts.name;
    cmin = vs.markers(dp.marker(a)).counts.cmin(:,dp.stat(a));
    cmax = vs.markers(dp.marker(a)).counts.cmax(:,dp.stat(a));
        
elseif isequal(dp.category(a),3)
    % Do density
    values = vs.markers(dp.marker(a)).density.values(:,dp.stat(a));
    title = vs.markers(dp.marker(a)).density.title(dp.stat(a));   
    name = vs.markers(dp.marker(a)).density.name;
    cmin = vs.markers(dp.marker(a)).density.cmin(:,dp.stat(a));
    cmax = vs.markers(dp.marker(a)).density.cmax(:,dp.stat(a));
end

end



function [slice,colors,title] = retrieve_slice(dp,vs)

% Get which axis
a = dp.pos;

% Displaying annotation volume
if dp.marker(a) == 0
    slice = vs.av(:,:,dp.z);
    colors = vs.av_colors;
    
    s1 = 1:length(vs.info);
    s1 = s1(ismember(s1,vs.indexes));
    s1 = s1(ismember(s1,vs.av_idx{dp.z}));
    
    acronyms = retrieve_acronyms(vs,s1);
    dp.acronyms = acronyms;
    colors = [1,1,1;colors(s1,:)];
    
    slice(~ismember(slice,s1)) = 1;
    
    [~, ~,i3] = unique(slice);
    s1 = [1,s1];
    s1 = 1:length(s1);
    slice2 = s1(i3);
    slice = reshape(slice2, size(slice));    
        
    %c1 = vs.av_colors(s1,:);
    %img = vs.av(:,:,dp.z);
    %img(~ismember(img,s1)) = 0;
    
    
    %colormap(c1)
    
    title = "Annotation Volume";
else
    % Displaying some statistics
    slice = double(vs.av(:,:,dp.z));
    idxs = vs.av_idx{dp.z};
    if isequal(dp.category(a),1)
        % Do volume
        stats = vs.markers(1).Volume.values(:,dp.stat(a));
        title = vs.markers(1).Volume.title(dp.stat(a));
        
    elseif isequal(dp.category(a),2)
        % Do counts
        stats = vs.markers(dp.marker(a)).counts.values(:,dp.stat(a));
        title = vs.markers(dp.marker(a)).counts.title(dp.stat(a));
        
    elseif isequal(dp.category(a),3)
        %  Do density
        stats = vs.markers(dp.marker(a)).density.values(:,dp.stat(a));
        title = vs.markers(dp.marker(a)).density.title(dp.stat(a));
        
    else
        error("Invalid stats category selected")
    end

    % Recolor slice according to stats
    for i = 1:length(idxs)
        slice(slice == idxs(i)) = stats(idxs(i));
    end
    
    % Get colormap
    colors = vs.markers(dp.marker(a)).colors{dp.stat(a)};
end

end

function p = add_colorbar(p,vs,dp)

if isequal(dp.init,'false')
    if dp.marker(dp.pos) == 0
        update_annotation_colorbar(p,vs,dp);
    end
    
elseif dp.pos == 1
% Set colorbar
c = colorbar(p);
c.Limits = [dp.cmin(1) dp.cmax(1)];

% Set colorbar positions
c.Position(1) = 0.0550;
c.Position(3) = 0.015;

    % If annotation volume, add tick labels with acronyms for major areas
    if dp.marker(1) == 0
        c.Ticks = vs.av_cmap_pos;
        %c.TickLabels = vs.av_cmap_labels;

        c.TickLabels = [];
        c.Limits = [1 dp.cmax(1)];
        labels = {vs.info.acronym};

        
        idx = vs.av_idx{dp.z};
        idx = idx(vs.bin_idx{dp.z});

        colors = vs.av_colors(idx,:);
        
        [~,b] = unique(colors,'rows');
        c.Ticks = sort(b+1);
        c.TickLabels = {vs.info(sort(idx(b)')).acronym};
        
       %%%%%
        %s1 = 1:length(vs.info);
        %s1 = s1(ismember(s1,vs.indexes));
        %s1 = s1(ismember(s1,vs.av_idx{dp.z}));
        %c1 = vs.av_colors(s1,:);
        %img = vs.av(:,:,dp.z);
        %img(~ismember(img,s1)) = 0;
        %colormap(c1)    
        %s2 = 1:length(s1);
        %colormap s2
        %c.Limits = [min(s2) max(s2)];

        %idxs = ismember(vs.indexes,vs.av_idx{dp.z});
        %idxs = vs.av_idx{dp.z}(idxs);
        %c.Ticks = 1:length(idxs);
        %c.Limits = [min(idxs) max(idxs)];  
        %c.Ticks = round((idxs(1:end-1)+idxs(2:end))/2);
        %c.Ticks = unique(idxs)-1;
    end 
end

if dp.pos == 2
    % Set colorbar
    c = colorbar(p);

    % Set colorbar positions
    c.Position(1) = 0.93;
    c.Position(3) = 0.015;
    c.Limits = [dp.cmin(2) dp.cmax(2)];

    % If annotation volume, add tick labels with acronyms for major areas
    if dp.marker(2) == 0
        c.Ticks = vs.av_cmap_pos;
        c.TickLabels = vs.av_cmap_labels;
    end
end

end

function update_annotation_colorbar(p,vs,dp)


idx = vs.av_idx{dp.z};
idx = idx(vs.bin_idx{dp.z});

colors = vs.av_colors(idx,:);

a = 1;

end




function updateSlice(f, evt, vs)

dp = get(f, 'UserData');
pre_dp_pos = dp.pos;


% scroll through slices
dp.z = dp.z - evt.VerticalScrollCount;

if dp.z>size(vs.av,3); dp.z = 1; end %wrap around
if dp.z<1; dp.z = size(vs.av,3); end %wrap around

vs.z = dp.z;

% Update image p1
dp.pos = 1;
[slice,~,~] = retrieve_slice(dp,vs);
f.Children(4).Children(3).CData = single(slice);
add_colorbar(f.Children(1),vs,dp);

% Update alpha p1
a = ~ismember(single(vs.av(:,:,dp.z)),vs.indexes);
f.Children(4).Children(2).CData = a;

if dp.marker(dp.pos) ~= 0
    f.Children(4).Children(2).AlphaData = a;
else
    f.Children(4).Children(2).AlphaData = a*dp.alpha_val;
end

% Overlay boundaries p1
boundaries = vs.boundaries(:,:,dp.z)*-1;
f.Children(4).Children(1).CData = boundaries;
boundaries = boundaries<0;
f.Children(4).Children(1).AlphaData = boundaries;

% Update image p2
dp.pos = 2;
[slice,~,~] = retrieve_slice(dp,vs);
f.Children(2).Children(3).CData = single(slice);

% Update alpha p2
a = ~ismember(single(vs.av(:,:,dp.z)),vs.indexes);
f.Children(2).Children(2).CData = a;
if dp.marker(dp.pos) ~= 0
    f.Children(2).Children(2).AlphaData = a;
else
    f.Children(2).Children(2).AlphaData = a*dp.alpha_val;
end

% Overlay boundaries p2
boundaries = vs.boundaries(:,:,dp.z)*-1;
f.Children(2).Children(1).CData = boundaries;
boundaries= boundaries<0;
f.Children(2).Children(1).AlphaData = boundaries;


% Update AP coordinates
dp.ptext.String = sprintf('Slice:%d/%d \t %.2f AP',dp.z,size(vs.av,3),vs.ap(dp.z));

dp.pos = pre_dp_pos;

set(f, 'UserData', dp);

end

function p = sanitize_plot(p, idx)

set(p, 'Units', 'normalized');
set(p,'xtick',[])
set(p,'ytick',[])

set(p,'XColor','none')
set(p,'YColor','none')

daspect([1 1 1])

if idx == 1
    set(p, 'Position', [0.075, 0.025, 0.425, 0.85]);
else
    set(p, 'XDir','reverse')
    set(p, 'Position', [0.5, 0.025, 0.425, 0.85]);
end

end

function acronyms = retrieve_acronyms(vs,indexes)

acronyms = string({vs.info.acronym});
acronyms = acronyms(indexes);

end

function dp = set_display_parameters(vs, key)

% User data structure
dp.z = key{1}(1);
dp.marker = [key{1}(2),key{2}(2)];
dp.category = [key{1}(3),key{2}(3)];
dp.stat = [key{1}(4),key{2}(4)];
dp.pos = 1;
dp.alpha = 'false';
dp.alpha_val = 1;
dp.init = 'true';
dp.imgL = zeros(size(vs.av(:,:,1)));
dp.imgR = zeros(size(vs.av(:,:,1)));

end


