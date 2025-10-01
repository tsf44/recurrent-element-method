function phi = ModeShapes(x,K,Pm_1,Qm_1,tau,epsilon,psi)
    % Define input argument
    KX = K * x;
    % Define the mode shape, phi(x)
    phi = tau .* (Pm_1*cos(KX) - 1j*epsilon*psi*Qm_1*sin(KX));
end
