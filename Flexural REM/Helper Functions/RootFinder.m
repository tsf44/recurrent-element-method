function beta = RootFinder(m_ratio, n)
    % Initialize value of 'beta'
    beta = zeros(size(m_ratio));
    % Find root for each element
    for m_idx = 1:numel(m_ratio)
        % Transcendental equation
        z = @(x) -m_ratio(m_idx) + (1/x) ...
            * ((1 + cosh(x)*cos(x)) ...
            / (cosh(x)*sin(x) - sinh(x)*cos(x)));
        % Initial value, x0, dependent on the mode number, n
        switch n
            case 1
                x0 = pi/2;
            case 2
                x0 = 3*pi/2;
            case 3
                x0 = 5*pi/2;
            case 4
                x0 = 7*pi/2;
        end
        % Find the nth root to the transcendental equation
        beta(m_idx) = fzero(z, x0);
    end
end
