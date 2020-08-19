function Iout = apply_guided_filter(I)
% Apply image guided filtering. This preserves edges slightly better than
% gaussian filtering

neighborhood = [4,4];
amount = 1E8;

Iout = imguidedfilter(I,'DegreeofSmoothing',1E8,'NeighborhoodSize',[4 4]);

end