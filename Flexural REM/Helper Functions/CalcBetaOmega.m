function [omega, beta] = CalcBetaOmega(n, EI, mu, BCcase)
    % Define 'm_ratio'
    m_ratio = Calculate_M_Ratio(mu);
    % Define value of 'beta'
    switch BCcase
        case "Clamped-Free"
            beta = RootFinder(m_ratio, n);
        case "Free-Free"
            beta = RootFinder(m_ratio, n+1);
        case "Pinned-Pinned"
            beta = (n*pi) * ones(size(EI));
    end
    % Calculate the oscillating frequency of the element
    omega = (beta).^2 .* sqrt(EI ./ mu); % [rad/s]
end
