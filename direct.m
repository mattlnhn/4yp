function Tii = direct(Tb, Ti, h, dt, geom, mat, A)
%DIRECT
% Tb        boundary temperatures
% Ti        current temperature distribution
% h         thermal contact conductance
% dt        time step
% geom      {1} lengths of each section [l1, l2, l3]
%           {2} dx
%           {3} nodes in each section = geometry{1}./geometry{2}
%           {4} total nodes = sum(geometry{3})
% mat       3x3x3 matrix
%           dim 1: material 1, 2, 3
%           dim 2: k, rho, cp
%           dim 3: a, b, c (where f(T) = a + bT + cT^2)
% A         precalculated tri-diagonal A

    %% friendly variable names
    dx = geom{2};
    n1 = geom{3}(1);
    n2 = geom{3}(2);
    n3 = geom{3}(3);
    n = geom{4};

    %% temperature dependent properties
    % mean temperatures
    Tm = [sum(Ti(1:n1))/n1 sum(Ti(n1+1:n1+n2))/n2 sum(Ti(n1+n2+1:n))/n3]; 

    k1 = mat(1, 1, 1) + mat(1, 1, 2)*Tm(1) + mat(1, 1, 3)*Tm(1)^2;
    rho1 = mat(1, 2, 1) + mat(1, 2, 2)*Tm(1) + mat(1, 2, 3)*Tm(1)^2;
    cp1 = mat(1, 3, 1) + mat(1, 3, 2)*Tm(1) + mat(1, 3, 3)*Tm(1)^2;
    alpha1 = k1/(rho1*cp1);
    tau1 = alpha1*dt/dx^2;

    k2 = mat(2, 1, 1) + mat(2, 1, 2)*Tm(2) + mat(2, 1, 3)*Tm(2)^2;
    rho2 = mat(2, 2, 1) + mat(2, 2, 2)*Tm(2) + mat(2, 2, 3)*Tm(2)^2;
    cp2 = mat(2, 3, 1) + mat(2, 3, 2)*Tm(2) + mat(2, 3, 3)*Tm(3)^2;
    alpha2 = k2/(rho2*cp2);
    tau2 = alpha2*dt/dx^2;

    k3 = mat(3, 1, 1) + mat(3, 1, 2)*Tm(3) + mat(3, 1, 3)*Tm(3)^2;
    rho3 = mat(3, 2, 1) + mat(3, 2, 2)*Tm(3) + mat(3, 2, 3)*Tm(3)^2;
    cp3 = mat(3, 3, 1) + mat(3, 3, 2)*Tm(3) + mat(3, 3, 3)*Tm(3)^2;
    alpha3 = k3/(rho3*cp3);
    tau3 = alpha3*dt/dx^2;

    %% build tau matrix
    tau = zeros(n);
    tau(1:n+1:n*n1) = tau1;
    tau(n*n1 + n1 +1:n+1:n*(n1+n2)) = tau2;
    tau(n*(n1+n2) + (n1+n2) + 1:n+1:n*n) = tau3;

    %% build b vector
    % temperature boundary conditions
    b = zeros(n, 1);
    b(1) = 2*Tb(1);
    b(end) = 2*Tb(2);

    %% adjust A matrix
    % see paper for maths
    phi1a = 2*k1/(h*dx);
    phi2a = 2*k2/(h*dx);
    xi1a = 2*phi1a/(1-(1+phi1a)*(1+phi2a));
    xi2a = 2*phi2a/(1-(1+phi1a)*(1+phi2a));
    A(n1, n1) = xi2a-1;
    A(n1, n1+1) = -xi2a;
    A(n1+1, n1) = xi1a;
    A(n1+1, n1+1) = xi1a-1;

    phi1b = 2*k2/(h*dx);
    phi2b = 2*k3/(h*dx);
    xi1b = 2*phi1b/(1-(1+phi1b)*(1+phi2b));
    xi2b = 2*phi2b/(1-(1+phi1b)*(1+phi2b));
    A(n1+n2, n1+n2) = xi2b-1;
    A(n1+n2, n1+n2+1) = -xi2b;
    A(n1+n2+1, n1+n2) = -xi1b;
    A(n1+n2+1, n1+n2+1) = xi1b-1;

    %% update
    Tii = tau*(A*Ti + b) + Ti;

end

