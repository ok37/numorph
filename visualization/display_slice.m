function [p1,p2] = display_slice2(Left,Right,z)

f = figure;
%screen_size = get(0,'ScreenSize');
%screen_size = screen_size(3:4)./[2560 1440];
%set(f,'Name','Atlas Viewer','Position', [1000*screen_size(1) 250*screen_size(2) 1100*screen_size(1) 850*screen_size(2)])
movegui(f,'onscreen')
p1 = subplot(1,2,1);
imagesc(Left(:,:,z))
p2 = subplot(1,2,2);
imagesc(Right(:,:,z))
set(p2, 'XDir','reverse')
set(p1, 'Units', 'normalized');
set(p2, 'Units', 'normalized');
set(p1,'xtick',[])
set(p2,'ytick',[])
set(p1,'ytick',[])
set(p2,'xtick',[])
set(p1, 'Position', [0.05, 0.05, 0.45, 0.8]);
set(p2, 'Position', [0.5, 0.05, 0.45, 0.8]);

linkaxes([p1 p2],'xy')


end