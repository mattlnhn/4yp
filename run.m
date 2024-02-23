%RUN

filename = "";

%% materials
mat = zeros(3, 3, 3);
% [a_k b_k c_k;
% a_rho b_rho c_rho;
% a_cp b_cp c_cp]
% where f(T) = a + bT + cT^2
h25     = [9.9905357 .0205437 -.000003;
            9070 0 0;
            396.5228931 .2075422 .0000134];
inco718 = [9.5164989 .0216787 -.0000039;
            8190 0 0;
            363.8195515 .1233661 .0000527];
mat(1, :, :) = inco718;
mat(2, :, :) = h25;
mat(3, :, :) = inco718;

%% geometry
geom = cell(5, 1);
% lengths of each section [l1, l2, l3]
% for instrumented sections, length from interface to final thermocouple
% for central disc, height of disc.
geom{1} = [33.5e-3, 10e-3, 33.5e-3];
% dx
geom{2} = ;
% nodes in each section = geom{1}./geom{2}
geom{3} = geom{1}./geom{2};
% total nodes = sum(geom{3})
geom{4} = sum(geom{3});
% thermocouple nodes #IMPORTANT: THINK ABOUT COORDINATE SYSTEM
geom{5} = [xcoord2node([], geom{2})];

%% parameters for transient solvers
% time step
dt = ;
% [r, epsilon, hi]
% hi = initial value of h
params = ;
% [RTOLh, RTOLe, pmax]
% convergence criteria
conv = ;

%% execute
h = inverse(filename, geom, mat, dt, params, conv);

function n = xcoord2node(x, dx)
    n = ((x./dx) + 1)./2;
end
