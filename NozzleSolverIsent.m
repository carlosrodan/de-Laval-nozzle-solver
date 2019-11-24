clear all;close all;clc;

N =62; %Number of nodes
Nt = 180; %Number of time steps
gam = 1.4;
R = 287.058; %J/kg.k
C = 1; %Courant number
xp = linspace(0,3,N);
incrXp = 3/N;
Ap = 1+2.2*(xp-1.5).^2; %A'= A
d_dxp_lnAp = 4.4*(xp-1.5)./(1+2.2*(xp-1.5).^2);
Apmin = min(Ap);

%% Analitical solution
sign = 1; %Mach icreases.
Mini = 0.5;
[M_an,p_p0_an,rho_rho0_an,T_T0_an,V_a0_an,mp_an] = isentFlow(Ap,gam,sign,Mini);

%% Anderson's numerical solution
load('Anderson.mat');
xp_And = Anderson_Isent(:,1);
rhop_And = Anderson_Isent(:,3);
Vp_And = Anderson_Isent(:,4);
Tp_And = Anderson_Isent(:,5);
pp_And = Anderson_Isent(:,6);
M_And = Anderson_Isent(:,7);
mp_And = Anderson_Isent(:,8);

%% Boundary conditions at the inflow
rho1bc = 1;
T1bc = 1;

%% Inital conditons
%Find where xp < 0.5
pos05 = max(find(xp<0.5));
%Find where xp < 1.5
pos15 = max(find(xp<1.5));
% -- Primtve varables (non-dmensonal) --
% 0<xp<0.5
rho1 = 1*ones(1,pos05);
T1 = 1*ones(1,pos05);
% 0.5<xp<1.5
rho2 = 1 - 0.366*(xp(pos05+1:pos15)-0.5);
T2 = 1 - 0.167*(xp(pos05+1:pos15)-0.5);
% 1.5<xp<3
rho3 = 0.634 - 0.3879*(xp(pos15+1:end)-1.5);
T3 = 0.833 - 0.3507*(xp(pos15+1:end)-1.5);

rhop = [rho1,rho2,rho3];
Tp = [T1,T2,T3];

% -- Dependent variables --
U1 = rhop.*Ap;
Vp = 0.59./U1;
U2 = 0.59*ones(1,N);
ep = Tp;
U3 = rhop.*(ep/(gam-1)+gam/2*Vp.^2).*Ap;
%The initial conditions are correct!. Same as in Table 7.9. P.346 Anderson

%% Time marching loop
tp(1) = 0;
for t = 1:Nt
    %%  Time step
    ap = sqrt(Tp);
    incr_tp = min(C*incrXp./(ap+Vp));
    %time
    tp(t+1) = tp(t)+incr_tp;
    %% Flux terms (pure form)
    [F1,F2,F3,J2] = fluxTerms(U1,U2,U3,gam,d_dxp_lnAp);
    %% Predictor step (Forward difference)
    for i = 2:N-1
        d_dtp_U1_t(i) = - (F1(i+1)-F1(i))/incrXp;
        d_dtp_U2_t(i) = - (F2(i+1)-F2(i))/incrXp+J2(i);
        d_dtp_U3_t(i) = - (F3(i+1)-F3(i))/incrXp;
    end
    %Predicted values
    U1_pred = U1(1:N-1) + d_dtp_U1_t * incr_tp;
    U2_pred = U2(1:N-1) + d_dtp_U2_t * incr_tp;
    U3_pred = U3(1:N-1) + d_dtp_U3_t * incr_tp;
    %% Flux terms (pure form)
    [F1_pred,F2_pred,F3_pred,J2_pred] = fluxTerms(U1_pred,U2_pred,U3_pred,gam,d_dxp_lnAp(2:N));
    %% Corrector step (Rearward difference)
    for i = 2:N-1
       d_dtp_U1_t_dt(i) = - (F1_pred(i)-F1_pred(i-1))/incrXp;
       d_dtp_U2_t_dt(i) = - (F2_pred(i)-F2_pred(i-1))/incrXp + J2_pred(i);
       d_dtp_U3_t_dt(i) = - (F3_pred(i)-F3_pred(i-1))/incrXp;
    end
    %% Averaged derivatives
    d_dt_U1_av = 0.5*(d_dtp_U1_t+d_dtp_U1_t_dt);
    d_dt_U2_av = 0.5*(d_dtp_U2_t+d_dtp_U2_t_dt);
    d_dt_U3_av = 0.5*(d_dtp_U3_t+d_dtp_U3_t_dt);
    
    incrU1 = d_dt_U1_av * incr_tp;
    incrU2 = d_dt_U2_av * incr_tp;
    incrU3 = d_dt_U3_av * incr_tp;
    Error1(t+1) = mean(abs(incrU1));
    Error2(t+1) = mean(abs(incrU2));
    Error3(t+1) = mean(abs(incrU3));
    %% Final corrected values (from i=2 to i=N-1)
    U1 = U1(1:N-1) + incrU1;
    U2 = U2(1:N-1) + incrU2;
    U3 = U3(1:N-1) + incrU3;
    %--Apply Boundary conditions at the inlet (i=1)--
    U1(1) = rho1bc*Ap(1);
    U2(1) = 2*U2(2)-U2(3);
    Vp(1) = U2(1)/U1(1);
    U3(1) = U1(1)*(T1bc/(gam-1)+gam/2*Vp(1)^2);
    %--Linear extrapolation for values at outlet (i=N)--
    U1(N) = 2*U1(N-1)-U1(N-2);
    U2(N) = 2*U2(N-1)-U2(N-2);
    U3(N) = 2*U3(N-1)-U3(N-2);
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
    plot(xp,pp,xp,p_p0_an,xp_And,pp_And,'o','MarkerSize',3,'Linewidth',1);
    xlabel('x/L'); ylabel('p/p_0'); grid minor;
    title('Pressure ratio');
    %Density ratio
    subplot(3,3,3);
    plot(xp,rhop,xp,rho_rho0_an,xp_And,rhop_And,'o','MarkerSize',3,'Linewidth',1);
    xlabel('x/L'); ylabel('\rho/\rho_0'); grid minor;
    title('Density ratio');
    %Temperature ratio
    subplot(3,3,4);
    plot(xp,Tp,xp,T_T0_an,xp_And,Tp_And,'o','MarkerSize',3,'Linewidth',1);
    xlabel('x/L'); ylabel('T/T_0'); grid minor;
    title('Temperature ratio');
    %Velocity ratio
    subplot(3,3,5);
    plot(xp,Vp,xp,V_a0_an,xp_And,Vp_And,'o','MarkerSize',3,'Linewidth',1);
    xlabel('x/L'); ylabel('V/a_0'); grid minor;
    title('Velocity ratio'); 
    %Mach number
    subplot(3,3,6);
    plot(xp,M,xp,M_an,xp_And,M_And,'o','MarkerSize',3,'Linewidth',1);
    xlabel('x/L'); ylabel('M'); grid minor;
    title('Mach number');
    %Mass flow rate ratio
    subplot(3,3,7);
    plot(xp,mp,xp,mp_an,xp_And,mp_And,'o','MarkerSize',3,'Linewidth',1);
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
    set(leg,'Fontsize',9,'Location', 'northeast','Orientation','horizontal');
    pause(0.001);

    
end
