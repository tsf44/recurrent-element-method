function [P0, Q0, R0, S0] = InitializePQRS(BCcase, n, beta)
    % Values of P, Q, R, S depend on the boundary condition, 
    % and certain parameters are dependent on the value of 
    % 'n' in order to normalize the flexural displacement.
    switch BCcase % Boundary Conditions of the Beam
        case "Clamped-Free"
            P0 = 0;
            Q0 = 0;
            % Value of R0 dependent on the mode number
            switch n
                case 1
                    R0 = 1;
                case 2
                    R0 = 1.05579;
                case 3
                    R0 = 1.58916;
            end
            S0 = -R0 * ((cosh(beta) + cos(beta))...
                    ./ (sinh(beta) + sin(beta)));
        case "Free-Free"
            % Value of P0 dependent on the mode number
            switch n
                case 1
                    P0 = 1;
                case 2
                    P0 = 1.05579;
                case 3
                    P0 = 1.58916;
            end
            Q0 = -P0 * ((cosh(beta) - cos(beta))...
                    ./ (sinh(beta) - sin(beta)));
            R0 = 0;
            S0 = 0;
        case "Pinned-Pinned"
            P0 = 0;
            % Value of Q0 dependent on the mode number
            switch n
                case 1
                    Q0 = 1;
                case 2
                    Q0 = 1.05579;
                case 3
                    Q0 = 1.58916;
            end
            R0 = 0;
            S0 = -Q0;     
    end
end
