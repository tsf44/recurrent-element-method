function dphi_dx = Phi_1st_Deriv(x,K,Pm_1,Qm_1,tau,epsilon,psi)
    % Define input argument
    KX = K * x;
    % Define the mode shape, phi(x)
    dphi_dx = - K*tau .* (Pm_1*sin(KX) + 1j*epsilon*psi*Qm_1*cos(KX));
end
