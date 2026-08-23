
%% define the parameters
alpha_1 = 0.01;   
alpha_2 = 0.024;   
alpha_H = 10^(-8);  
alpha_3 = 0.693;   
alpha_4 = 0.01;  
alpha_KU= (9.7*10^9)^(-1);  
alpha_I = (5*10^7)^(-1);  
gamma_1 = 10^(-8);  
gamma_2 = 0.00001;  
gamma_3 = 10^(-5); 
mu_1    = 0.1; 
mu_H    = 0.0005;
mu_2    = 0.2;
mu_I    = 0.024;
mu_KI   = 3.66;
mu_3    = 0.3;
mu_KU   = 0.01;
mu_D    = 1.5;  
mu_V    = 1;
q       = 10^5;
delta_1 = 0.054;
delta_2 = 200;
delta_3 = 0.027;
delta_4 = 0.0178;
Beta_1  = 10^(-10);
Beta_2  = 10^(-10);
Beta_3  = 10^(-10);
Beta_4  = 50.8;
psi     = 0;
phi_V   = 0;
phi_D   = 0;

param = [alpha_1, alpha_2, alpha_H, alpha_3, alpha_4,...
    alpha_KU, alpha_I, gamma_1, gamma_2, gamma_3,...
    mu_1, mu_H, mu_2, mu_I, mu_KI, mu_3, mu_KU,...
    mu_D, mu_V, q, delta_1, delta_2, delta_3, delta_4, ...
    Beta_1, Beta_2, Beta_3, Beta_4, psi, phi_D, phi_V] ;

init = [10^6, 10^3, 0, 10^5, 0, 0];
tspan = [1:0.1:7];
ChemoDose = 1.6;
opt = odeset('RelTol',1e-12,'AbsTol',1e-14);
[time1,y1] = ode45(@(t,x)DissModel(t,x,param),tspan, init,opt);

init = y1(end,:);
init(5) = y1(end,5)+ChemoDose; 
tspan = [7:0.1:14];
[time2,y2] = ode45(@(t,x)DissModel(t,x,param),tspan, init,opt);

init = y2(end,:);
init(5) = y2(end,5)+ChemoDose; 
tspan = [14:0.1:21];
[time3,y3] = ode45(@(t,x)DissModel(t,x,param),tspan, init,opt);

init = y3(end,:);
init(5) = y3(end,5)+ChemoDose;
tspan = [21:0.1:30];
[time4,y4] = ode45(@(t,x)DissModel(t,x,param),tspan, init,opt);

time = cat(1,time1(1:end-1,:),time2(1:end-1,:),time3(1:end-1,:),time4(1:end-1,:));
y = cat(1,y1(1:end-1,:),y2(1:end-1,:),y3(1:end-1,:),y4(1:end-1,:));

% Plot the results
figure;
plot(time, y(:,1),"Color",[0.18,0.55,0.34],'LineWidth',2);
hold on
plot(time, y(:,2),"Color",[0,0,1],'LineWidth',2);
plot(time, y(:,4),"Color",[1,0,0],'LineWidth',2);
xlabel('Time (Days)');
legend('Healthy Cells','Immune Cells','Cancer Cells');
title('Three Treatments of Chemotherapy');
set(gcf,'Color',[1,1,1])
set(gca,'LineWidth',2,'FontSize',20)
