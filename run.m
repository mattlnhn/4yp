%RUN

filename = "";

mat = zeros(3, 3, 3);
mat(1, :, :) = [; ;];
mat(2, :, :) = [; ;];
mat(3, :, :) = [; ;];

geom = cell(5, 1);
geom{1} = ;
geom{2} = ;
geom{3} = ;
geom{4} = ;
geom{5} = ;

dt = ;

params = ;

conv = ;

h = inverse(filename, geom, mat, dt, params, conv);