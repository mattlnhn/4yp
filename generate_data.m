%GENERATE_DATA
%% init
clear

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
geom{2} = 5e-5;
% nodes in each section = geom{1}./geom{2}
geom{3} = geom{1}./geom{2};
% total nodes = sum(geom{3})
geom{4} = sum(geom{3});
% thermocouple nodes #IMPORTANT: THINK ABOUT COORDINATE SYSTEM
geom{5} = [xcoord2node([10 20 30 32 45 47 57 67]*1e-3, geom{2})];

%% friendly names
dx = geom{2};
n1 = geom{3}(1);
n2 = geom{3}(2);
n3 = geom{3}(3);
n = geom{4};

%% parameters for transient solver
Tmin = 25;
Tmax = 150;
% properties at Tmin
k1_Tmin = [1 Tmin Tmin^2]*reshape(mat(1, 1, :), [], 1, 1);
rho1_Tmin = [1 Tmin Tmin^2]*reshape(mat(1, 2, :), [], 1, 1);
cp1_Tmin = [1 Tmin Tmin^2]*reshape(mat(1, 3, :), [], 1, 1);
k2_Tmin = [1 Tmin Tmin^2]*reshape(mat(2, 1, :), [], 1, 1);
rho2_Tmin = [1 Tmin Tmin^2]*reshape(mat(2, 2, :), [], 1, 1);
cp2_Tmin = [1 Tmin Tmin^2]*reshape(mat(2, 3, :), [], 1, 1);
alpha1_Tmin = k1_Tmin/(rho1_Tmin*cp1_Tmin);
alpha2_Tmin = k2_Tmin/(rho2_Tmin*cp2_Tmin);

k1_Tmax = [1 Tmax Tmax^2]*reshape(mat(1, 1, :), [], 1, 1);
rho1_Tmax = [1 Tmax Tmax^2]*reshape(mat(1, 2, :), [], 1, 1);
cp1_Tmax = [1 Tmax Tmax^2]*reshape(mat(1, 3, :), [], 1, 1);
k2_Tmax = [1 Tmax Tmax^2]*reshape(mat(2, 1, :), [], 1, 1);
rho2_Tmax = [1 Tmax Tmax^2]*reshape(mat(2, 2, :), [], 1, 1);
cp2_Tmax = [1 Tmax Tmax^2]*reshape(mat(2, 3, :), [], 1, 1);
alpha1_Tmax = k1_Tmax/(rho1_Tmax*cp1_Tmax);
alpha2_Tmax = k2_Tmax/(rho2_Tmax*cp2_Tmax);

% time step
dt = .5*dx^2./[alpha1_Tmin alpha2_Tmin alpha1_Tmax alpha2_Tmax];
dt = min(dt);

%% build basic A matrix
% relatively costly, only needs doing once
A = diag(ones(n-1, 1), -1) + diag(-2*ones(n, 1), 0) + ...
    diag(ones(n-1, 1), 1);
A(1, 1) = -3;
A(n, n) = -3;

%% establish initial temperature distribution
% Ti = Tmin*ones(n, 1);
% Tb = [150 Tmin];
% E = 1e99*ones(n, 1);
% tolE = 1e-6;
% 
% % clock
% elapsed = 0;
% 
% while any E above tolerance continue iterating
% while any(E>tolE)
%     elapsed = elapsed + dt;
%     fprintf("%.6f s\n", elapsed)
%     Tii = direct(Tb, Ti, 1e12, dt, geom, mat, A);
%     E = (Tii-Ti)./Ti;
%     Ti = Tii;
% end
% save("Ti.mat", "Ti")

load("Ti.mat")

%% function for h
T = 120; % seconds, period for simulation
s = ceil(T/dt);
h = interp1([0 T/2 T], [10000 1000 10000], linspace(0, T, s), "nearest");

%% generate and save noisy temperature data
% iterations
p = 0;
pmax = s;

% save properties
f_save = 20; % Hz, save frequency
T_save = zeros(n, T*f_save);

% clock
elapsed = 0; % time
c = 0; % counter

while p < pmax
    % increment counter
    p = p + 1;
    % increment clock
    elapsed = elapsed + dt;
    
    % update temperature
    Tii = direct(Tb, Ti, h(p), dt, geom, mat, A);
    Ti = Tii;

    % if time to save
    if elapsed >= (1/f_save)
        % increment clock counter
        c = c + 1;
        % reset clock
        elapsed = 0;
        % add noise drawn from zero mean gaussian where 3*sigma = 1.5 degC,
        % save
        T_save(:, c) = Tii + normrnd(0, 0.5, [n 1]);
    end
end

save("T-"+string(datetime("now","Format","yyyyMMdd-HHmmss"))+".mat",...
    "T_save")

%FUNCTIONS
function n = xcoord2node(x, dx)
    n = ((x./dx) + 1)./2;
end

