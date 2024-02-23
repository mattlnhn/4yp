function h_out = inverse(filename, geom, mat, dt, params, conv)
%INVERSE
% geom      {1} lengths of each section [l1, l2, l3]
%           {2} dx
%           {3} nodes in each section = geom{1}./geom{2}
%           {4} total nodes = sum(geom{3})
%           {5} thermocouple nodes #IMPORTANT: THINK ABOUT COORDINATES
% mat       3x3x3 matrix
%           dim 1: material 1, 2, 3
%           dim 2: k, rho, cp
%           dim 3: a, b, c (where f(T) = a + bT + cT^2)
% dt        time step
% params    inverse solver parameters [r, epsilon, hi]
% conv      convergence criteria [RTOLh, RTOLe, pmax]
    %% friendly variable names
    n = geom{4};
    r = params(1);
    epsilon = params(2);
    hi = params(3);
    RTOLh = conv(1);
    RTOLe = conv(2);
    pmax = conv(3);

    %% build basic A matrix
    % relatively costly, only needs doing once
    A = diag(ones(n-1, 1), -1) + diag(-2*ones(n, 1), 0) + ...
        diag(ones(n-1, 1), 1);
    A(1, 1) = -3;
    A(n, n) = -3;

    %% import data
    dat = readtable(filename);
    M = size(dat, 2); % total no. of time steps
    h_out = zeros(M-r+1);
    h_out(1) = hi;

    %% initial temperature distribution
    Ti = interp1([1; geom{5}'; n], ...
        [dat.T0(1); dat.T1(1); dat.T2(1); dat.T3(1); dat.T4(1); ...
        dat.T5(1); dat.T6(1); dat.T7(1); dat.T8(1); dat.T9(1)]);

    %% loop
    % index mm = m-1 i.e. previous time step.
    % annoying but reduces confusion
    for mm = 1:M-r
        done = 0; % convergence flag
        p = 0; % iteration counter
        E = 1e99; % initial error

        t = dat.time(mm:mm+r); % times corresponding to m-1, m,...,m+r-1
        s = ceil((t - t(1))/dt); % step #s where Y updates

        % matrix of temp. vectors for each step
        % dim 1: temperature
        % dim 2: time
        T = zeros(N, s(end)+1);

        T(:, 1) = Ti; % initial temp
        Tdh = T; % matrix of temp. vectors for h update

        h = h_out(mm);

        % measured temperature
        % dim 1: temperature
        % dim 2: time
        Y = [dat.T1(mm:mm+r)'; dat.T2(mm:mm+r)'; dat.T3(mm:mm+r)';
            dat.T4(mm:mm+r)'; dat.T5(mm:mm+r)'; dat.T6(mm:mm+r)';
            dat.T7(mm:mm+r)'; dat.T8(mm:mm+r)'];
        
        % boundary temperatures
        Tb = [dat.T0(mm:mm+r)'; dat.T9(mm:mm+r)'];
        % linear interpolation
        Tb = interp1(1:r, Tb, linspace(1, r, s(end)));

        while done == 0
            p = p + 1;
            Eprev = E; % previous error

            % index ii = time steps after the initial temp.
            for ii = 1:s(end)
                T(:, ii+1) = direct(Tb(ii), T(:, ii), h, dt, geom, mat, A);
                Tdh(:, ii+1) = direct(Tb(ii), T(:, ii), h, dt, geom, mat, A);
            end

            % temperatures at sensor locations and measurment times
            Ts = T(geom{5}, s(2:end)+1);
            Tdhs = Tdh(geom{5}, s(2:end)+1);
    
            % sensitivity coefficients
            X = (Tdhs - Ts)./(epsilon*h);
            % update
            dh = sum((Y-Tdhs).*X, 'all')/sum(X.^2, 'all');
            h = h + dh;
    
            % error
            E = sum((Y-Ts).^2, 'all');
    
            % convergence criteria
            % relative change in h
            if abs(dh/h) < RTOLh
                done = 1;
            end
            % relative change in error
            if abs((E-Eprev)/E) < RTOLe
                done = 1;
            end
            % timeout
            if p >= pmax
                done = 1;
            end

        end

        h_out(mm+1) = h;

        for n = 1:s(2)
            Ti = direct(Tb(n), Ti, h, dt, geom, mat, A);
        end
    end
end
