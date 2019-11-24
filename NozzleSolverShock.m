clear all;close all;clc;

N =61; %Number of nodes
% Nt = 580 for Artifical Viscosity
Nt = 600; %Number of time steps
gam = 1.4;
R = 287.058; %J/kg.k
C = 1; %Courant number
Cx = 0.3; %0.01< Cx <0.3

xp = linspace(0,3,N);
incrXp = 3/N;
Ap = 1+2.2*(xp-1.5).^2; %A'= A
d_dxp_lnAp = 4.4*(xp-1.5)./(1+2.2*(xp-1.5).^2);

%% Boundary conditions
%Inlet
rho1bc = 1;
T1bc = 1;
%Outlet
ppbc = 0.6784;

%% Analytical results
%--Position of shock wave---
pe_p01 = ppbc;
peAe_pe0Amin = pe_p01*Ap(end);
Me = fsolve(@(M)1/M*(2/(gam+1))^((gam+1)/(2*(gam-1)))*(1+(gam-1)/2*M^2)^(-1/2)-peAe_pe0Amin,0.1);
pe_p0e = (1+(gam-1)/2*Me^2)^(-3.5);
p02_p01 = pe_p01/pe_p0e;
M1 = fsolve(@(M) p02_p01-(((gam+1)*M^2)/(2+(gam-1)*M^2))^(gam/(gam-1))*((gam+1)/(2*gam*M^2-gam+1))^(1/(gam-1)),2);
Ap_shock = sqrt(1/M1^2*(2/(gam+1)*(1+(gam-1)/2*M1^2))^((gam+1)/(gam-1)));
%Look for index from the second half of the nozzle
[out,idx] = sort(abs(Ap-Ap_shock));
pos = idx(2);
%Position of the shock. x/L:
xp_shock = xp(pos);

%New Amin (throat area) Amin2.
Ae_Amine = 1/Me*(2/(gam+1)*(1+(gam-1)/2*Me^2))^((gam+1)/(2*(gam-1))); %Amin = A*, A*_e = A*_2
%Ae_Amine = Ap(end)*Amin1/Amin2
Amin1_Amin2 = Ae_Amine/Ap(end);
Ap2 = Ap*Amin1_Amin2;

%Flow properties before the shock wave (isentropic relations)
sign = 1; %Mach icreases.
Mini = 0.5;
[M_1,p1_p0,rho1_rho0,T1_T0,V1_a0,mp1] = isentFlow(Ap(1:pos),gam,sign,Mini);

%Flow properties right after the shock wave
M2 = sqrt((2+(gam-1)*M1^2)/(2*gam*M1^2-gam+1));
p2_p1 = (2*gam*M1^2-gam+1)/(gam+1);
T2_T1 = (2*gam*M1^2-gam+1)*(2+(gam-1)*M1^2)/((gam+1)^2*M1^2);
rho2_rho1 = (gam+1)*M1^2/(2+(gam-1)*M1^2);
V2_V1 = M2/M1.*sqrt(T2_T1);

p2_p0 = p1_p0(pos)*p2_p1;
T2_T0 = T1_T0(pos)*T2_T1;
T02_T01 = 1;
rho2_rho0 = rho1_rho0(pos)*rho2_rho1;
V2_a0 = V1_a0(pos)*V2_V1;
mp2 = rho2_rho1*V2_V1*rho1_rho0(pos)*V1_a0(pos)*Ap_shock;

%Flow properties after the shock wave (isentropic relations)
sign = -1; %Mach decreases.
Mini = M2;
[M_2,p_p02,~,T_T02,V_a02,mp_2] = isentFlow(Ap2(pos+1:end),gam,sign,Mini);
%Correction p/p01 = p/p02* p02/p01 (T01 = T02)
p_p01_2 = p_p02*p02_p01;
%The density is calculated from the corrected pressure ratio
rho_rho01_2 = p_p01_2.^(1/gam);

M_an = [M_1(1:pos-1),M2,M_2];
p_p0_an = [p1_p0(1:pos-1),p2_p0,p_p01_2];
T_T0_an = [T1_T0(1:pos-1),T2_T0,T_T02];
rho_rho0_an = [rho1_rho0(1:pos-1),rho2_rho0,rho_rho01_2];
V_a0_an = [V1_a0(1:pos-1),V2_a0,V_a02];
mp_an = [mp1(1:pos-1),mp2,mp_2];

% figure(1);
% plot(xp,M_an,xp,p_p0_an,xp,T_T0_an,xp,rho_rho0_an,xp,V_a0_an,xp,mp_an);
% legend('Mach','p','T','rho','V','mp');
% stop;

%% Anderson's numerical solution
load('Anderson.mat');
xp_And = Anderson_Shock(:,2);
rhop_And = Anderson_Shock(:,4);
Vp_And = Anderson_Shock(:,5);
Tp_And = Anderson_Shock(:,6);
pp_And = Anderson_Shock(:,7);
M_And = Anderson_Shock(:,8);
mp_And = Anderson_Shock(:,9);

%% Inital conditons
%Find where xp < 0.5
pos05 = max(find(xp<0.5));
%Find where xp < 1.5
pos15 = max(find(xp<1.5));
%Find where xp < 2.1
pos21 = max(find(xp<2.1));
% -- Primtve varables (non-dmensonal) --
% 0<xp<0.5
rho1 = 1*ones(1,pos05);
T1 = 1*ones(1,pos05);
% 0.5<xp<1.5
rho2 = 1 - 0.366*(xp(pos05+1:pos15)-0.5);
T2 = 1 - 0.167*(xp(pos05+1:pos15)-0.5);
% 1.5<xp<2.1
rho3 = 0.634-0.702*(xp(pos15+1:pos21)-1.5);
T3 = 0.833-0.4908*(xp(pos15+1:pos21)-1.5);
% 2.1<xp<3
rho4 = 0.5892+0.10228*(xp(pos21+1:end)-2.1);
T4 = 0.93968+0.0622*(xp(pos21+1:end)-2.1);

rhop = [rho1,rho2,rho3,rho4];
Tp = [T1,T2,T3,T4];
pp = rhop.*Tp;

% -- Dependent variables --
U1 = rhop.*Ap;
Vp = 0.59./U1;
U2 = 0.59*ones(1,N);
ep = Tp;
U3 = rhop.*(ep/(gam-1)+gam/2*Vp.^2).*Ap;
%The initial conditions are correct!. Same as in Table 7.9. P.346 Anderson

%% Preallocation
d_dtp_U1_t = zeros(1,N-1);
d_dtp_U2_t = zeros(1,N-1);
d_dtp_U3_t = zeros(1,N-1);
d_dtp_U1_t_dt = zeros(1,N-1);
d_dtp_U2_t_dt = zeros(1,N-1);
d_dtp_U3_t_dt = zeros(1,N-1);
%% Time marching loop
tp(1) = 0;
for t = 1:Nt
    %%  Time step
    ap = sqrt(Tp);
    incr_tp = min(C*incrXp./(ap+Vp));
    tp(t+1) = tp(t)+incr_tp;
    % -- Flux terms (pure form) --
    [F1,F2,F3,J2] = fluxTerms(U1,U2,U3,gam,d_dxp_lnAp);
    %% Predictor step (Forward difference)
    for i = 2:N-1
        d_dtp_U1_t(i) = - (F1(i+1)-F1(i))/incrXp;
        d_dtp_U2_t(i) = - (F2(i+1)-F2(i))/incrXp+J2(i);
        d_dtp_U3_t(i) = - (F3(i+1)-F3(i))/incrXp;
        % -- Artificial viscosites --
        S_num = Cx*abs(pp(i+1)-2*pp(i)+pp(i-1))/(pp(i+1)+2*pp(i)+pp(i-1));
        S1(i) = S_num*(U1(i+1)-2*U1(i)+U1(i-1));
        S2(i) = S_num*(U2(i+1)-2*U2(i)+U2(i-1));
        S3(i) = S_num*(U3(i+1)-2*U3(i)+U3(i-1));
%         S1 = 0;S2 = 0;S3 = 0;
    end
    %Predicted values
    U1_pred = U1(1:N-1) + d_dtp_U1_t * incr_tp + S1;
    U2_pred = U2(1:N-1) + d_dtp_U2_t * incr_tp + S2;
    U3_pred = U3(1:N-1) + d_dtp_U3_t * incr_tp + S3;
    [F1_pred,F2_pred,F3_pred,J2_pred] = fluxTerms(U1_pred,U2_pred,U3_pred,gam,d_dxp_lnAp(1:N-1));
    %--Apply Boundary conditions at the inlet (i=1)--
    U1_pred(1) = rho1bc*Ap(1);
    U2_pred(1) = 2*U2_pred(2)-U2_pred(3);
    U3_pred(1) = U1_pred(1)*(T1bc/(gam-1)+gam/2*Vp(1)^2);
    %--Linear extrapolation for values at outlet (i=N)--
    U1_pred(N) = 2*U1_pred(N-1)-U1_pred(N-2);
    U2_pred(N) = 2*U2_pred(N-1)-U2_pred(N-2);
    U3_pred(N) = 2*U3_pred(N-1)-U3_pred(N-2);
   % -- pp_pred ---
    rhop_pred = U1_pred./Ap;
    Vp_pred = U2_pred./U1_pred; 
    Tp_pred = (gam-1)*(U3_pred./U1_pred-gam/2*Vp_pred.^2);
    pp_pred = rhop_pred.*Tp_pred;
    %% Corrector step (Rearward difference)
    for i = 2:N-1
       d_dtp_U1_t_dt(i) = - (F1_pred(i)-F1_pred(i-1))/incrXp;
       d_dtp_U2_t_dt(i) = - (F2_pred(i)-F2_pred(i-1))/incrXp + J2_pred(i);
       d_dtp_U3_t_dt(i) = - (F3_pred(i)-F3_pred(i-1))/incrXp;
       S_num = Cx*abs(pp_pred(i+1)-2*pp_pred(i)+pp_pred(i-1))/(pp_pred(i+1)+2*pp_pred(i)+pp_pred(i-1));
       S1(i) = S_num*(U1_pred(i+1)-2*U1_pred(i)+U1_pred(i-1));
       S2(i) = S_num*(U2_pred(i+1)-2*U2_pred(i)+U2_pred(i-1));
       S3(i) = S_num*(U3_pred(i+1)-2*U3_pred(i)+U3_pred(i-1));
% S1 = 0;S2 = 0;S3 = 0;
    end
    %% Averaged derivatives
    d_dt_U1_av = 0.5*(d_dtp_U1_t+d_dtp_U1_t_dt);
    d_dt_U2_av = 0.5*(d_dtp_U2_t+d_dtp_U2_t_dt);
    d_dt_U3_av = 0.5*(d_dtp_U3_t+d_dtp_U3_t_dt);
    %% Final corrected values (from i=2 to i=N-1)
    incrU1 = d_dt_U1_av * incr_tp + S1;
    incrU2 = d_dt_U2_av * incr_tp + S2;
    incrU3 = d_dt_U3_av * incr_tp + S3;
    Error1(t+1) = mean(abs(incrU1));
    Error2(t+1) = mean(abs(incrU2));
    Error3(t+1) = mean(abs(incrU3));
    U1 = U1(1:N-1) + incrU1;
    U2 = U2(1:N-1) + incrU2;
    U3 = U3(1:N-1) + incrU3;
    %--Apply Boundary conditions at the inlet (i=1)--
    U1(1) = rho1bc*Ap(1);
    U2(1) = 2*U2(2)-U2(3);
    U3(1) = U1(1)*(T1bc/(gam-1)+gam/2*Vp(1)^2);
    %--Linear extrapolation for values at outlet (i=N)--
    U1(N) = 2*U1(N-1)-U1(N-2);
    U2(N) = 2*U2(N-1)-U2(N-2);
    Vp(N) = U2(N)/U1(N);
    %The value of U3(N) is determined by the specified exit pressure ratio.
    U3(N) = ppbc*Ap(N)/(gam-1)+gam/2*U2(N)*Vp(N);
    %--------------------------------------------
    %% Final primitive variables
    rhop = U1./Ap; %density ratio
    Vp = U2./U1;   %Velocity ratio
    Tp = (gam-1)*(U3./U1-gam/2*Vp.^2); %Temperature ratio
    ep = Tp; %Energy ratio
    pp = rhop.*Tp; %Pressure ratio
    M = (Vp./ap); %Mach number
    mp = rhop.*Ap.*Vp; %Mass flow rate ratio

    %% Graphs
    scrsz = get(0,'ScreenSize');
    figure(1);
    %Area ratio
    subplot(3,3,1);
    plot(xp,Ap,'Linewidth',1);
    xlabel('x/L'); ylabel('A/A_*'); grid minor;
    title('Area ratio');
    %Pressure ratio
    subplot(3,3,2);
    plot(xp,pp,xp,p_p0_an,'Linewidth',1);
%     plot(xp,pp,xp,p_p0_an,xp_And,pp_And,'o','MarkerSize',3,'Linewidth',1);
    xlabel('x/L'); ylabel('p/p_0'); grid minor;
    title('Pressure ratio');
    %Density ratio
    subplot(3,3,3);
    plot(xp,rhop,xp,rho_rho0_an,'Linewidth',1);
%     plot(xp,rhop,xp,rho_rho0_an,xp_And,rhop_And,'o','MarkerSize',3,'Linewidth',1);
    xlabel('x/L'); ylabel('\rho/\rho_0'); grid minor;
    title('Density ratio');
    %Temperature ratio
    subplot(3,3,4);
    plot(xp,Tp,xp,T_T0_an,'Linewidth',1);
%     plot(xp,Tp,xp,T_T0_an,xp_And,Tp_And,'o','MarkerSize',3,'Linewidth',1);
    xlabel('x/L'); ylabel('T/T_0'); grid minor;
    title('Temperature ratio');
    %Velocity ratio
    subplot(3,3,5);
    plot(xp,Vp,xp,V_a0_an,'Linewidth',1);
%     plot(xp,Vp,xp,V_a0_an,xp_And,Vp_And,'o','MarkerSize',3,'Linewidth',1);
    xlabel('x/L'); ylabel('V/a_0'); grid minor;
    title('Velocity ratio'); 
    %Mach number
    subplot(3,3,6);
    plot(xp,M,xp,M_an,'Linewidth',1);
%     plot(xp,M,xp,M_an,xp_And,M_And,'o','MarkerSize',3,'Linewidth',1);
    xlabel('x/L'); ylabel('M'); grid minor;
    title('Mach number');
    %Mass flow rate ratio
    subplot(3,3,7);
    plot(xp,mp,xp,mp_an,'Linewidth',1);
%     plot(xp,mp,xp,mp_an,xp_And,mp_And,'o','MarkerSize',3,'Linewidth',1);
    xlabel('x/L'); ylabel('$\dot{m}$', 'Interpreter','latex'); grid minor;
    title('Mass flow rate ratio'); ylim([0.5,1]);

    %Legend
    subplot(3,3,8);
    plot(1:10, nan(1,10), 1:10, nan(1,10),1:10, nan(1,10));
    axis off;
    [~, hobj, ~, ~] = legend({'Numerical','Analytical',"Anderson's numerical solution"},'Fontsize',9,'Location','north');
    hl = findobj(hobj,'type','line');
    set(hl,'LineWidth',1);
    
    %Residuals
    subplot(3,3,9);
    plot(tp,Error1,tp,Error2,tp,Error3);
    xlabel('time (tp)'); ylabel('Residuals'); grid minor;
    leg = legend({'incrU1','incrU2','incrU3'});
    leg.ItemTokenSize = [10,30];
    set(leg,'Fontsize',9,'Location', 'northoutside','Orientation','horizontal');

    %Display current iteration
    iter = t
    
end

%% Error calculation vs. analytical (Residual Sum of Squares) 
rho_err = sum(abs(rhop-rho_rho0_an))
p_err = sum(abs(pp-p_p0_an))
T_err = sum(abs(Tp-T_T0_an))
V_err = sum(abs(Vp-V_a0_an))
M_err = sum(abs(M-M_an))
mp_err = sum(abs(mp-mp_an))
    
    
%% Save final results
% save('NumericalData_ArtfVisc06784.mat','rhop','Vp','Tp','pp','M','mp','p_p0_an','rho_rho0_an','T_T0_an','V_a0_an','M_an','mp_an','rho_err','p_err','T_err','V_err','M_err','mp_err');
% save('NumericalData_NoArtfVisc.mat','rhop','Vp','Tp','pp','M','mp','rho_err','p_err','T_err','V_err','M_err','mp_err');

% %% Save all figures
% figure(2);
% %Area ratio
% plot(xp,Ap,'Linewidth',1);
% xlabel('x/L'); ylabel('A/A_*'); grid minor;
% title('Area ratio');
% print(gcf,'Ap.png','-dpng','-r300');
% 
% %Pressure ratio
% figure(3);
% plot(xp,pp,xp,p_p0_an,xp_And,pp_And,'o','MarkerSize',3,'Linewidth',1);
% xlabel('x/L'); ylabel('p/p_0'); grid minor;
% title('Pressure ratio');
% legend('Numerical','Analytical',"Anderson's numerical solution",'Location','best');
% print(gcf,'pp.png','-dpng','-r300');
% 
% %Density ratio
% figure(4);
% plot(xp,rhop,xp,rho_rho0_an,xp_And,rhop_And,'o','MarkerSize',3,'Linewidth',1);
% xlabel('x/L'); ylabel('\rho/\rho_0'); grid minor;
% title('Density ratio');
% legend('Numerical','Analytical',"Anderson's numerical solution",'Location','best');
% print(gcf,'rhop.png','-dpng','-r300');
% 
% %Temperature ratio
% figure(5);
% plot(xp,Tp,xp,T_T0_an,xp_And,Tp_And,'o','MarkerSize',3,'Linewidth',1);
% xlabel('x/L'); ylabel('T/T_0'); grid minor;
% title('Temperature ratio');
% legend('Numerical','Analytical',"Anderson's numerical solution",'Location','best');
% print(gcf,'Tp.png','-dpng','-r300');
% 
% %Velocity ratio
% figure(6);
% plot(xp,Vp,xp,V_a0_an,xp_And,Vp_And,'o','MarkerSize',3,'Linewidth',1);
% xlabel('x/L'); ylabel('V/a_0'); grid minor;
% title('Velocity ratio'); 
% legend('Numerical','Analytical',"Anderson's numerical solution",'Location','best');
% print(gcf,'Vp.png','-dpng','-r300');
% 
% %Mach number
% figure(7);
% plot(xp,M,xp,M_an,xp_And,M_And,'o','MarkerSize',3,'Linewidth',1);
% xlabel('x/L'); ylabel('M'); grid minor;
% title('Mach number');
% legend('Numerical','Analytical',"Anderson's numerical solution",'Location','best');
% print(gcf,'Mach.png','-dpng','-r300');
% 
% %Mass flow rate ratio
% figure(8);
% plot(xp,mp,xp,mp_an,xp_And,mp_And,'o','MarkerSize',3,'Linewidth',1);
% xlabel('x/L'); ylabel('$\dot{m}$', 'Interpreter','latex'); grid minor;
% title('Mass flow rate ratio'); ylim([0.5,0.75]);
% legend('Numerical','Analytical',"Anderson's numerical solution",'Location','best');
% print(gcf,'mp.png','-dpng','-r300');
% 
% %Residuals
% figure(9);
% plot(tp,Error1,tp,Error2,tp,Error3);
% xlabel('time (tp)'); ylabel('Residuals'); grid minor;
% legend('Error1','Error2','Error3','Location','Best');
% print(gcf,'residuals.png','-dpng','-r300');
