function [t, yt, yxt, x_elem, phi_elem, omega] = ...
    FlexBend_REM_General(X, x_elem, fs, options)
% --- DESCRIPTION ---------------------------------------- %
% FlexBend_REM_NonUniform() is a recurrent element model 
% (REM) of a flexuralbending beam of normalized length L,  
% where each element has its ownflexural rigidity, EI, and 
% mass per unit length, mu. The response of the beam is 
% simulated over time using the Euler-Bernoulli beam 
% equation. Each element is modeled using a relative length 
% of 0 to Lx, where M * Lx = L,though each element is 
% thought of as having unit length.
%
% The position and slope at element joints are equal over   
% time, such thatposition and slope at the x=0 boundary of 
% the m-th element are equal to position and slope at the 
% x=Lx boundary of the (m-1)-th element. The curvature and 
% shear force have similar boundary conditions, with the 
% exception that a change in flexural rigidity will cause a 
% small discontinuity in the curvature and shear force. Said 
% change in flexural rigidity can be caused by having a 
% multi-material structure along thelength of the beam, or 
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
% AUTHOR: (c) Oct. 2023 Trent Furlong
% -------------------------------------------------------- %
    arguments
        X       % Relative element length per element {cell}
        x_elem  % Element postion along x=(0:dx:L-dx) {cell}
        fs      % Sampling frequency, [Hz]
        options.n = 1;
        options.EI = 100*ones(1,numel(X));
        options.mu = 1*ones(1,numel(X));
        options.T = 1;
        options.sigma = 0;
        options.animateTF = false;
        options.BCcase = "Clamped-Free"
        options.StrainTF = false;
        options.normalize = true;
    end

% --- Parse 'options' argument --------------------------- %
M = numel(X);    % Number of elements
n = options.n;   % n-th mode shape
EI = options.EI; % Flexural rigidity array,  [Pa*m^4]
mu = options.mu; % Linear mass density array, [kg/m]
T = options.T;   % Duration, [s]
sigma = options.sigma; % Noise standard deviation, [WU]
animateTF = options.animateTF;  % Animate vibration response
BCcase = options.BCcase;        % Boundary conditions
StrainTF = options.StrainTF;    % Express results in strain
normalize = options.normalize;  % Normalize results

% --- Ensure dimensionality of 'EI' and 'mu' match 'M' --- %
if numel(EI) == 1
    EI = EI * ones(1, M);
end
if numel(mu) == 1
    mu = mu * ones(1, M);
end

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

% --- Define and preallocate loop variables -------------- %
% Response over time and space
yxt = cell([numel(t), numel(X)]);
% Mode shape response of each element
phi_elem = cell(size(X));

% --- Calculate response over time for each element ------ %
for ii = 1:length(t)
    % Reset parameters at each time step
    wx = sigma * zeros(size(X)); % White noise (spatial)
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
            phi_m = ModeShape(X{mm}, beta_m, ...
                P, Q, R, S, tau, epsilon, psi);
                
            % Define the spatial derivatives of 'phi_m'
            dphi_dx_m = Phi_1st_Deriv(X{mm}, beta_m, ...
                P, Q, R, S, tau, epsilon, psi);
            d2phi_dx2_m = Phi_2nd_Deriv(X{mm}, beta_m, ...
                P, Q, R, S, tau, epsilon, psi);
            d3phi_dx3_m = Phi_3rd_Deriv(X{mm}, beta_m, ...
                P, Q, R, S, tau, epsilon, psi);
        
        % --- ELEMENTS 2:M ------------------------------- %
        else
            % --- Update parameters P,Q,R,S -------------- %
            % Get last value of 'phi_m' and each of its 
            % derivatives 
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
            phi_m = ModeShape(X{mm}, beta_m, ...
                P, Q, R, S, tau, epsilon, psi);
            
            % Define the spatial derivatives of 'phi_m'
            dphi_dx_m = Phi_1st_Deriv(X{mm} ,beta_m, ...
                P, Q, R, S, tau, epsilon, psi);
            d2phi_dx2_m = Phi_2nd_Deriv(X{mm}, beta_m, ...
                P, Q, R, S, tau, epsilon, psi);
            d3phi_dx3_m = Phi_3rd_Deriv(X{mm}, beta_m, ...
                P, Q, R, S, tau, epsilon, psi);
        end

        % --- Populate element response variable --------- %
        if StrainTF == true
            yxt{ii, mm} = sin(omega(mm)*t(ii)) ...
                * dphi_dx_m(1:end-1) + wx(:,mm);
        else
            yxt{ii, mm} = sin(omega(mm)*t(ii)) ...
                * phi_m(1:end-1) + wx(:,mm);
        end

        % --- Populate mode shape ------------------------ %
        if ii == 1
            if StrainTF == true
                phi_elem{mm} = dphi_dx_m(1:end-1);
            else
                phi_elem{mm} = phi_m(1:end-1);
            end
        end
        
    end
end

% --- Normalize the mode shape by the max amplitude ------ %
if normalize == true
    % Find the maximum absolute value of the displacement 
    % over time
    temp_cell = cellfun(@abs, yxt, 'UniformOutput', false); 
    maxVal = max(max(cellfun(@max, temp_cell)));
else
    maxVal = 1;
end
yxt = cellfun(@(x) x / maxVal, yxt, 'UniformOutput', false);
phi_elem = cellfun(@(x) x / maxVal, phi_elem, ...
    'UniformOutput', false);

% --- Animate the vibration response over time (T/F) ----- %
if animateTF == true
    for ii = 1:length(t)
        cla
        hold on
        cellfun(@plot,x_elem,yxt(ii,:))
        ax = gca;
        ax.LineWidth = 2;
        ax.FontSize = 18;
        xlabel('Length [Normalized]')
        ylabel('Displacement [Normalized]')
        ylim([-1.15 1.15])
        xlim([0 1])
        hold off;
        drawnow
    end
end

yt = NaN; % (CURRENTLY UNUSED IN THIS VERSION OF CODE)
% -------------------------------------------------------- %
end
