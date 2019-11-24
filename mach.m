function F = mach(M,Ap,gam)
 F = (Ap)^2 - 1/M^2*(2/(gam+1)*(1+(gam-1)/2*M^2))^((gam+1)/(gam-1));
end