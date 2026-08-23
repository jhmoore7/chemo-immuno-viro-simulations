function dXdt = DissModel(t,var, param)
%% define the parameters
alpha_1        = param(1);  % tumor cell growth rate 
alpha_2   = param(2);  % HSPC growth rate 
alpha_H  = param(3);  % M1 growth rate
alpha_3  = param(4);  % M2 growth rate 
alpha_4    = param(5);  % M2-induced tumor growth rate
alpha_KU        = param(6);  % tumor cell carrying capacity
alpha_I  = param(7);  % Loss of tumor cells due to M1
gamma_1      = param(8);  % HSPC half-saturation constant(growth)
gamma_2      = param(9);  % tumor death rate
gamma_3      = param(10); % HSPC death rate
mu_1      = param(11); % M2 death rate
mu_H    = param(12);
mu_2    = param(13);
mu_I    = param(14);
mu_KI   = param(15);
mu_3    = param(16);
mu_KU   = param(17);
mu_D    = param(18);
mu_V    = param(19);
q       = param(20);
delta_1 = param(21);
delta_2 = param(22);
delta_3 = param(23);
delta_4 = param(24);
Beta_1  = param(25);
Beta_2  = param(26);
Beta_3  = param(27);
Beta_4  = param(28);
psi     = param(29);
phi_V   = param(30);
phi_D   = param(31);

%% define the variables
H  = var(1);     % Tumor cells
I  = var(2);     % HSPCs from bone marrow
K_I = var(3);     % M1 macrophages
K_U = var(4);     % M2 macrophages
D   = var(5);
V   = var(6);
%% introduce the equations
dH   = alpha_1*H*(1-alpha_H*H)-gamma_1*K_U*H-mu_1*D*H-mu_H*H;
dI   = q+psi+alpha_2*I*(1-alpha_I*I)-gamma_2*K_U*I-mu_2*D*I-mu_I*I;
dK_I = delta_1*K_U*V/(delta_2+K_U)-delta_3*K_I*I-delta_4*K_I*D-mu_KI*K_I;
dK_U = alpha_3*K_U*(1-alpha_KU*(K_I+K_U))-delta_1*K_U*V/(delta_2+K_U)-gamma_3*K_U*I-mu_3*K_U*D-mu_KU*K_U;
dD   = phi_D-Beta_1*I*D-Beta_2*H*D-Beta_3*(K_U+K_I)*D-mu_D*D;
dV   = phi_V+Beta_4*(delta_3*I+delta_4*D)*K_I-mu_V*V;
dXdt = [dH;dI;dK_I;dK_U;dD;dV];
end