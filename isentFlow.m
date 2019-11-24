function [M,p_p0,rho_rho0,T_T0,V_a0,mp] = isentFlow(Ap,gam,sign,Mini)
    %This function calculates the isentropic flow properties across a
    %nozzle with area ratio Ap (A/Amin). 
    N = length(Ap);
    for i = 1:N
        M(i) = fsolve(@(M)mach(M,Ap(i),gam),Mini);
        Mini=M(i)+0.05*sign;
%         Mini=M(i);
    end
    p_p0 = (1+(gam-1)/2*M.^2).^(-gam/(gam-1));
    rho_rho0 = (1+(gam-1)/2*M.^2).^(-1/(gam-1));
    T_T0 = (1+(gam-1)/2*M.^2).^(-1);
    V_a0 = M.*sqrt(T_T0);
    mp = rho_rho0.*V_a0.*Ap;
end