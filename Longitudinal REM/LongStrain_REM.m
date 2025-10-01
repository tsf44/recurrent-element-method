function [t,yt,yxt,x_elem,phi_elem,omega] = LongStrain_REM(M,fs,options)
% --- DESCRIPTION ------------------------------------------------------- %
% LongStrain_REM() is a recurrent element model (REM) for modeling 
% longitudinal strain of a beam of normalized length L, where each element 
% has its own effective stiffness, 'Y', and effective mass, 'rho'. The 
% response of the beam is simulated over time as a solution to the linear, 
% lossless wave equation. Each element is modeled using a relative length 
% of 0 to Lx, where M * Lx = L, though each element is thought of as having
% unit length.
%
% The position at element joints are equal over time, such that the
% position and at the x=0 boundary of the m-th element is equal to the
% position and at the x=Lx boundary of the (m-1)-th element. The strain has
% a similar boundary condition constraint, with the exception that a change 
% in effective stiffness will cause a small discontinuity in the strain. 
% Said change in effective stiffness can be caused by having a 
% multi-material structure along the length of the beam, or due to damage 
% (e.g., cracking) in the structure.
%
% The natural frequencies of each element are calculated using embedded
% transcendental equations that look for the root value to said equation
% corresponding the n-th mode.
%
% The model is effective at modeling longitudinal strain for a beam having 
% uniform material properties. It's results are nearly identical with the
% theoretical model and has been tested up to the third mode shape, but 
% additional modes can theoretically be obtained using this code. The REM's
% performance is also independent of the number of elements used to define
% the REM, with the exception that it will take additional computation time
% to calculate the mode shapes for M-number of elements.
%
% For beams having non-uniform material properties, it is assumed that the
% REM will not provide exact mode shapes (given that the REM uses only the
% previous element to influence the shape of the current element), but it
% is expected that the REM may be a simple, but useful, approximation. More
% accurate structure models can be made using FEM or BEM techniques, but at
% the expense of higher computational cost.
% ----------------------------------------------------------------------- %
% AUTHOR: (c) July 2023 Trent Furlong
% ----------------------------------------------------------------------- %
    arguments
        M   % Number of elements and sensors
        fs  % Sampling frequency, [Hz]
        options.n = 1;
        options.L = 1;
        options.dx = 0.005;
        options.Y = 100*ones(1,M);
        options.SA = 1*ones(1,M);
        options.rho = 1*ones(1,M);
        options.T = 1;
        options.sigma = 0;
        options.animateTF = false;
        options.BCcase = "Free-Free";
        options.StrainTF = false;
        options.normalize = true;
        options.SLocs = (1/M) * (1:M).' - (1/(2*M));
    end

% --- Parse 'options' argument ------------------------------------------ %
n = options.n;                  % nth mode shape
L = options.L;                  % Length of the beam, [normalized length]
dx = options.dx;                % Spatial resolution, [normalized length]
Y = options.Y;                  % Young's modulus of each element, [Pa]
SA = options.SA;                % Surface area of each element, [m^2]
rho = options.rho;              % Density of each elements, [kg/m^3]
T = options.T;                  % Duration, [s]
sigma = options.sigma;          % White noise standard deviation, [WU]
animateTF = options.animateTF;  % Animate vibration response, {logical}
BCcase = options.BCcase;        % Boundary conditions, {string}
StrainTF = options.StrainTF;    % Express results in strain, {logical}
normalize = options.normalize;  % Normalize results, {logical}
SLocs = options.SLocs;          % Sensor locations relative to L

% --- Ensure dimensionality --------------------------------------------- %
if numel(Y) == 1
    Y = Y * ones(1,M);
end
if numel(SA) == 1
    SA = SA * ones(1,M);
end
if numel(rho) == 1
    rho = rho * ones(1,M);
end

% --- Spatial information ----------------------------------------------- %
x = (0:dx:L-dx).';
x_elem = reshape(x,[],M);
Lx = L/M;
X = (0:dx:Lx).';

% --- Time information -------------------------------------------------- %
dt = 1/fs;
t = (0:dt:T-dt).';

% --- Calculate the natural frequencies of the individual elements ------ %
c = sqrt(Y ./ rho);     % Quasi-longitudinal wave speed, [m/s]
switch BCcase
    case "Fixed-Fixed"
        k = n * pi;             % Wave number, [rad/m]
        P0 = 0;                 % Value of P for m = 1
        Q0 = 1;                 % Value of Q for m = 1
    case "Free-Free"
        k = n * pi;             % Wave number, [rad/m]
        P0 = 1;                 % Value of P for m = 1
        Q0 = 0;                 % Value of Q for m = 1
    case "Fixed-Free"
        k = (2*n - 1) * pi / 2; % Wave number, [rad/m]
        P0 = 0;                 % Value of P for m = 1
        Q0 = 1;                 % Value of Q for m = 1
end
k = k * ones(size(c));  % Make into an array equal in size to c
omega = c .* k;         % Angular frequency, [rad/s]
ES = Y .* SA;           % Young's modulus times surface area, [Pa*m^2]

% --- Animate the vibration response (T/F) ------------------------------ %
if animateTF == true
    figure;
end

% --- Define sensor locations, 'SLocs', w/in 'x' ------------------------ %
if iscolumn(SLocs)
    SLocs = SLocs.'; % Ensure 'SLocs' is a row vector
end
% Find x-values where x == SLocs
xLocs = round(x,4) == round(SLocs,4);
% Condense logical results into a single row vector
xLocs = sum(xLocs,2);
% Reshape 'xLocs' to match shape of 'x_elem' and convert to logical
xLocs = logical(reshape(xLocs, [], M));
% Repeat logical for each time step
xLocs = repmat(xLocs,[1 1 numel(t)]);

% --- Define and preallocate loop variables ----------------------------- %
yxt = zeros([size(x_elem,1), M, length(t)]);
phi_elem = zeros(size(x_elem));

% --- Calculate the response over time for each element ----------------- %
for ii = 1:length(t)
    % Reset parameters at each time step
    wx = sigma * randn(size(x_elem)); % White noise
    tau = 1;
    epsilon = 1;
    psi = 1;
    P = P0;
    Q = Q0;
    % --- Calculate the response of the m-th element -------------------- %
    for mm = 1:M
        % --- FIRST ELEMENT --------------------------------------------- %
        if mm == 1
            % Define the mode shape
            phi_m = ModeShapes(X(1:end-1),k(mm),P,Q,tau,epsilon,psi);
            % Define the spatial derivatives of the mode shape
            strain_m = Phi_1st_Deriv(X(1:end-1),k(mm),P,Q,tau,epsilon,psi);
        
        % --- ELEMENTS 2:M ---------------------------------------------- %
        else
            % --- Update parameters P,Q --------------------------------- %
            % Get last value of 'phi_m' and its first derivative 
            p = phi_m(end);
            q = strain_m(end);
            % Define P,Q based on derived definitions
            P = p;
            Q = (-1 / (1j * k(mm-1))) * q;

            % --- Update parameter values of 'tau','epsilon', and 'psi' ---
            if ii ~= 1
                tau = sin(omega(mm-1) * t(ii)) ...
                    / sin(omega(mm) * t(ii));
                epsilon = ES(mm-1) / ES(mm);
            end
            psi = k(mm-1) / k(mm);

            % Define the mode shape
            phi_m = ModeShapes(X,k(mm),P,Q,tau,epsilon,psi);
            % Define the spatial derivatives of the mode shape
            strain_m = Phi_1st_Deriv(X,k(mm),P,Q,tau,epsilon,psi);
        end

        % --- Populate element response variable ------------------------ %
        % Return either the longitudinal displacement or strain
        if length(phi_m) == length(X)
            if StrainTF == true
                % Strain
                yxt(:,mm,ii) = sin(omega(mm)*t(ii)) * strain_m(2:end) + wx(:,mm);
            else
                % Displacement
                yxt(:,mm,ii) = sin(omega(mm)*t(ii)) * phi_m(2:end) + wx(:,mm);
            end
        else
            if StrainTF == true
                % Strain
                yxt(:,mm,ii) = sin(omega(mm)*t(ii)) * strain_m + wx(:,mm);
            else
                % Displacement
                yxt(:,mm,ii) = sin(omega(mm)*t(ii)) * phi_m + wx(:,mm);
            end
        end

        % --- Populate mode shape --------------------------------------- %
        if ii == 1
            if mm == 1
                if StrainTF == true
                    phi_elem(:,mm) = strain_m;
                else
                    phi_elem(:,mm) = phi_m;
                end
            else
                if StrainTF == true
                    phi_elem(:,mm) = strain_m(2:end);
                else
                    phi_elem(:,mm) = phi_m(2:end);
                end
            end
        end

    end
end
% --- Normalize the mode shape by the maximum amplitude ----------------- %
% Find maximum value over space and time
if normalize == true
    maxVal = max(max(max(abs(yxt))));
    phi_elem = phi_elem / maxVal;
else
    maxVal = 1;
end
switch BCcase
    case "Fixed-Fixed"
        yxt = 1j * (yxt / maxVal);
    case "Free-Free"
        yxt = yxt / maxVal;
    case "Fixed-Free"
        yxt = 1j * (yxt / maxVal);
end
% [~,maxIdx] = max(yt_true(end/2,1:100));

% --- Animate the vibration response over time (T/F) -------------------- %
if animateTF == true
    for ii = 1:length(t)
        plot(x_elem,yxt(:,:,ii),'LineWidth',2)
        ax = gca;
        ax.FontSize = 18;
        xlabel('Length [Normalized]')
        if StrainTF == true
            ylabel('Strain [Normalized]')
        else
            ylabel('Displacement [Normalized]')
        end
        grid minor
        ylim([-1.15 1.15])
        xlim([0 1])
        drawnow
    end
end

% --- Define M sensor responses at x=Lx/2 of each element --------------- %
% Displacement of 'm' elements over time
% yt = squeeze(yxt(floor(end/2),:,:)).';
% --- Define sensor responses at x == SLocs ----------------------------- %
yt = yxt(xLocs);    % Returns matrix of size [SLocs * numel(t), 1]
% Reshape to size [numel(SLocs), numel(t)]
yt = reshape(yt,[numel(SLocs),numel(t)]);  
% Transpose; columns are sensor data
yt = yt.';          
% ----------------------------------------------------------------------- %
end
