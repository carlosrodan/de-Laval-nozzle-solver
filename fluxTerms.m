function [F1,F2,F3,J2] = fluxTerms(U1,U2,U3,gam,d_dxp_lnAp)
    F1 = U2;
    F2 = U2.^2./U1+(gam-1)/gam*(U3-gam/2*U2.^2./U1);
    F3 = gam*U2.*U3./U1-gam*(gam-1)/2*U2.^3./U1.^2;
    J2 = (gam-1)/gam*(U3-gam/2*U2.^2./U1).*d_dxp_lnAp;
end