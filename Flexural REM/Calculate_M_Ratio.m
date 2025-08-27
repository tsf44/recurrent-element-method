function m_ratio = Calculate_M_Ratio(mu)
    % Initialize 'm_ratio'
    m_ratio = zeros(size(mu));
    for idx = 1:numel(mu)
        % Calculate difference between current and past 'mu'
        if idx == 1
            % Initialize variables for first element
            m_last = 0;
            delta_mu = 0;
        else
            % Difference in 'mu'
            delta_mu = mu(idx) - mu(idx-1);
        end
        
        % Define m_ratio based on the value of delta_mu
        if delta_mu == 0
            % Set 'm_ratio' to last value of m_ratio
            m_ratio(idx) = m_last;
        else
            % Update value of 'm_ratio' if change in 'mu' 
            % from previous element
            m_ratio(idx) = delta_mu / mu(idx-1);
            % Update value of 'm_last' to the current value 
            % of 'm_ratio'
            m_last = m_ratio(idx);
        end
    end
end
