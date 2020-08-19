function visualize_flat_cortex()

load(fullfile('supplementary_data','flatviewCortex.mat'))
ftable = readtable(fullfile('annotations','harris_cortical_groupings.csv'));
indexes = ftable.index;

df = readtable(fullfile('evaluation','TCe_harris_cortical_groupings_stats.csv'));
indexes = df.index;


%v = aggregate_by_group(df,{'WT','TOP'},{'Counts'});

%wt = df.WT_Mean_Ctip2_Counts + df.WT_Mean_Cux1_Counts;
%top1 = df.TOP_Mean_Ctip2_Counts + df.TOP_Mean_Cux1_Counts;
wt = df.WT_Mean_Volume;
top1 = df.TOP_Mean_Volume;

x = [929,735,902,520,555,785,542,710,791,631,392,390,250,417,325,395,294,490,757,450,605,392,742,407,277,1040,1120,1098,1193,758,1050,909,557,228,499,760,870,1027,744,607,205,115,164];
y = [296,535,476,662,810,771,522,671,860,705,617,425,599,845,860,938,762,953,971,1056,1114,1162,1044,1013,1048,557,638,404,408,145,265,164,261,487,224,1195,1098,1060,913,938,882,773,814];

%df = df(ismember(df.index,indexes),:);
%stats = df(:,10:end);
%wt_idx = cellfun(@(s) contains(s,'WT'),stats.Properties.VariableNames);
%top_idx = ~wt_idx;

%stats = table2array(stats);
%wt_mean = mean(stats(:,wt_idx),2);
%top1_mean = mean(stats(:,top_idx),2);

prct = 100*(top1-wt)./wt;
%prct = log2(abs(top1-wt));

img = zeros(size(fc.flatCtx));
for i = 1:length(indexes)
    img(fc.flatCtx==indexes(i)) = prct(i);
end


img(img>0) = -1;
img(fc.boundaries) = 1;
f = imagesc(img);
colormap([brewermap(100,'*Blues');[1,1,1];[0,0,0]])
caxis([-100 1])

for i = 1:length(indexes)
    text(x(i),y(i),fc.acronym{i},'FontName','Arial')
end

h = colorbar;
h.Limits = [-100,0];
h.Ticks = [-100,0];
set(gca,'XTick',[], 'YTick', [])
daspect([1,1,1])

end