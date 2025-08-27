function [t, yt, yxt, x_elem, phi_elem, omega] = ...
    FlexuralBending_REM(M, fs, options)
% --- DESCRIPTION ---------------------------------------- %
% FlexuralBending_REM() is a recurrent element model (REM) 
% of a flexural bending beam of normalized length L, where 
% each element has its own flexural rigidity, EI, and mass 
% per unit length, mu. The response of the beam is simulated 
% over time using the Euler-Bernoulli beam equation. Each 
% element is modeled using a relative length of 0 to Lx, 
% where M * Lx = L,though each element is thought of as 
% having unit length.
%
% The position and slope at element joints are equal over 
% time, such that  position and slope at the x=0 boundary of 
% the m-th element are equal to position and slope at the 
% x=Lx boundary of the (m-1)-th element. The curvature and 
% shear force have similar boundary conditions, with the 
% exception that a change in flexural rigidity will cause a 
% small discontinuity in the curvature and shear force. Said 
% change in flexural rigidity can be caused by having a 
% multi-material structure along the length of the beam, or 
% due to damage (e.g., cracking) in the structure.
%
% The natural frequencies of each element are calculated 
% using embedded transcendental equations that look for the 
% root value to said equation corresponding the n-th mode.
%
% The model is effective at modeling a flexural bending beam 
% having uniform material properties. It's results are 
% nearly identical with the theoretical model and has been 
% tested up to the third mode shape, but additional modes 
% can theoretically be obtained using this code. The REM's 
% performance is also independent of the number of elements 
% used to define the REM, with the exception that it will 
% take additional computation time to calculate the mode 
% shapes for M-number of elements.
% -------------------------------------------------------- %
% AUTHOR: (c) July 2023 Trent Furlong
% -------------------------------------------------------- %
    arguments
        M   % Number of elements
        fs  % Sampling frequency, [Hz]
        options.n = 1;
        options.L = 1;
        options.dx = 0.005;
        options.EI = 100*ones(1,M);
        options.mu = 1*ones(1,M);
        options.T = 1;
        options.sigma = 0;
        options.animateTF = false;
        options.BCcase = "Clamped-Free"
        options.StrainTF = false;
        options.normalize = true;
        options.SLocs = (1/M) * (1:M) - (1/(2*M));
    end

% --- Parse 'options' argument --------------------------- %
n = options.n;    % n-th mode shape
L = options.L;    % Length of the beam, [normalized length]
dx = options.dx;  % Spatial resolution, [normalized length]
EI = options.EI;  % Flexural rigidity array, [Pa*m^4]
mu = options.mu;  % Linear mass density array, [kg/m]
T = options.T;    % Duration, [s]
sigma = options.sigma; % Noise standard deviation, [WU]
animateTF = options.animateTF; % Animate vibration response
BCcase = options.BCcase;       % Boundary conditions
StrainTF = options.StrainTF;   % Express results in strain
normalize = options.normalize; % Normalize results
SLocs = options.SLocs; % Sensor locations relative to L

% --- Ensure dimensionality of 'EI' and 'mu' match 'M' --- %
if numel(EI) == 1
    EI = EI * ones(1, M);
end
if numel(mu) == 1
    mu = mu * ones(1, M);
end

% --- Spatial information -------------------------------- %
x = (0:dx:L-dx).';          % Position along beam length
x_elem = reshape(x, [], M); % Element position along beam 
Lx = L/M;                   % 'Length' of each element
X = (0:dx:Lx).';            % Relative element position
X1 = X(1:end-1);            % Position of first element
% NOTE: For the variable X, the point at Lx is needed to get 
% the values of P, Q, R, S

% --- Time information ----------------------------------- %
dt = 1/fs;                  % Time step, [s]
t = (0:dt:T-dt).';          % Time array of duration T, [s]

% --- Calculate 'omega' and 'beta' for each element ------ %
% -> omega : Angular frequency array, [rad/s]
% ->  beta : Root values of transcendental eq., [unitless]
[omega, beta] = CalcBetaOmega(n, EI, mu, BCcase);

% --- Define the inital values of P,Q,R,S ---------------- %
[P0, Q0, R0, S0] = InitializePQRS(BCcase, n, beta(1));

% --- Animate the vibration response (T/F) --------------- %
if animateTF == true
    figure;
end

% --- Define sensor locations, 'SLocs', w/in 'x' --------- %
if iscolumn(SLocs)
    SLocs = SLocs.'; % Ensure 'SLocs' is a row vector
end
% Find x-values where x == SLocs
xLocs = round(x, 6) == round(SLocs, 6);
% Condense logical results into a single row vector
xLocs = sum(xLocs, 2);
% Reshape 'xLocs' to match 'x_elem' and convert to logical
xLocs = logical(reshape(xLocs, [], M));
% Repeat logical for each time step
xLocs = repmat(xLocs, [1, 1, numel(t)]);

% --- Define and preallocate loop variables -------------- %
% Response over time and space
yxt = zeros([size(x_elem, 1), M, length(t)]);
% Mode shape response of each element
phi_elem = zeros(size(x_elem));

% --- Calculate response over time for each element ------ %
for ii = 1:length(t)
    % Reset parameters at each time step
    wx = sigma*randn(size(x_elem)); % White noise (spatial)
    tau = 1;
    epsilon = 1;
    psi = 1;
    P = P0;
    Q = Q0;
    R = R0;
    S = S0;
    % --- Calculate the response of the m-th element ----- %
    for mm = 1:M
        % --- FIRST ELEMENT ------------------------------ %
        if mm == 1
            % Define current value of 'beta'
            beta_m = beta(mm);
            % Define the mode shape, 'phi_m'
            phi_m = ModeShape(X1, beta_m, ...
                P, Q, R, S, tau, epsilon, psi);
            % Define the spatial derivatives of 'phi_m'
            dphi_dx_m = Phi_1st_Deriv(X1, beta_m, ...
                P, Q, R, S, tau, epsilon, psi);
            d2phi_dx2_m = Phi_2nd_Deriv(X1, beta_m, ...
                P, Q, R, S, tau, epsilon, psi);
            d3phi_dx3_m = Phi_3rd_Deriv(X1, beta_m, ...
                P, Q, R, S, tau, epsilon, psi);
        
        % --- ELEMENTS 2:M ------------------------------- %
        else
            % --- Update parameters P,Q,R,S -------------- %
            % Get last value of 'phi_m' and its derivatives 
            p = phi_m(end);
            q = dphi_dx_m(end);
            r = d2phi_dx2_m(end);
            s = d3phi_dx3_m(end);
            % Define P,Q,R,S based on derived definitions
            P = p;
            Q = q / beta(mm-1);
            R = r / (beta(mm-1).^2);
            S = s / (beta(mm-1).^3);

            % Update values of 'tau','epsilon', and 'psi'
            if ii ~= 1
                tau = sin(omega(mm-1) * t(ii)) ...
                    / sin(omega(mm) * t(ii));
                epsilon = EI(mm-1) / EI(mm);
            end
            psi = beta(mm-1) / beta(mm);

            % Define current value of 'beta'
            beta_m = beta(mm);
            % Define the mode shape, 'phi_m'
            phi_m = ModeShape(X, beta_m, ...
                P, Q, R, S, tau, epsilon, psi);
            % Define the spatial derivatives of 'phi_m'
            dphi_dx_m = Phi_1st_Deriv(X, beta_m, ...
                P, Q, R, S, tau, epsilon, psi);
            d2phi_dx2_m = Phi_2nd_Deriv(X, beta_m, ...
                P, Q, R, S, tau, epsilon, psi);
            d3phi_dx3_m = Phi_3rd_Deriv(X, beta_m, ...
                P, Q, R, S, tau, epsilon, psi);
        end

        % --- Populate element response variable --------- %
        if length(phi_m) == length(X) % mm == 2:M
            % Get indices 2:end to not double count the 
            % x=0 point from the mth element and the x=Lx 
            % point from the (m-1)th element.
            if StrainTF == true
                yxt(:, mm, ii) = sin(omega(mm) * t(ii)) ...
                    * dphi_dx_m(2:end) + wx(:, mm);
            else
                yxt(:, mm, ii) = sin(omega(mm) * t(ii)) ...
                    * phi_m(2:end) + wx(:, mm);
            end
        else % mm == 1
            if StrainTF == true
                yxt(:, mm, ii) = sin(omega(mm) * t(ii)) ...
                    * dphi_dx_m + wx(:, mm);
            else
                yxt(:, mm, ii) = sin(omega(mm) * t(ii)) ...
                    * phi_m + wx(:, mm);
            end
        end

        % --- Populate mode shape ------------------------ %
        if ii == 1
            if mm == 1
                if StrainTF == true
                    phi_elem(:, mm) = dphi_dx_m;
                else
                    phi_elem(:, mm) = phi_m;
                end
            % Get indices 2:end to not double count the x=0 
            % point from the mth element and the x=Lx point 
            % from the (m-1)th element.
            else
                if StrainTF == true
                    phi_elem(:, mm) = dphi_dx_m(2:end);
                else
                    phi_elem(:, mm) = phi_m(2:end);
                end
            end
        end
        
    end
end

% --- Normalize the mode shape by the max amplitude ------- %
if normalize == true
    % Find the maximum absolute value of the displacement 
    % over time
    maxVal = max(max(max(abs(yxt))));
else
    maxVal = 1;
end
yxt = yxt / maxVal;             % Normalized modal response
phi_elem = phi_elem / maxVal;   % Normalized mode shapes

% --- Animate the vibration response over time (T/F) ----- %
if animateTF == true
    for ii = 1:length(t)
        plot(x_elem,yxt(:, :, ii), 'LineWidth', 2)
        ax = gca;
        ax.FontSize = 18;
        xlabel('Length [Normalized]')
        ylabel('Displacement [Normalized]')
        grid minor
        ylim([-1.15 1.15])
        xlim([0 1])
        drawnow
    end
end

% --- Define sensor responses at x == SLocs -------------- %
yt = yxt(xLocs); % Matrix of size [(SLocs * numel(t)), 1]
% Reshape to size [numel(SLocs), numel(t)]
yt = reshape(yt, [numel(SLocs), numel(t)]);  
% Transpose; columns are sensor data
yt = yt.';          
% -------------------------------------------------------- %
% -------------------------------------------------------- %
end
