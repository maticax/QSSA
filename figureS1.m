%% viral dynamics
clear all
clc
close all
set(0,'DefaultAxesFontSize',16)
 set(0,'DefaultAxesFontWeight','bold');
set(0, 'DefaultLineLineWidth', 3);
% viral dynamics 
% T'(t) = -bTV ; 
% I'(t) = bTV-delta*I ; 
% V'(t) = pI-cV ;

% qssa model 
% T'(t) =- b*p/c *T*I;
% I'(t) = p*b/c*T*I-delta*I;
% v(t) = p/c *I; 

% ill qssa model 
% T'(t) = - b*TV;
% v'(t) = (p*b)/c*TV-delat*v;
% I(t) = c/p *v;

% parameter set
V=300;
c_vec = [10 20 30 40 50 60 ];
%c_vec = [10 20 30 40 50 60];
delta =2.1;
%p_vec = [100 1000 5000 11000 15000 20000];
for i=1:length(c_vec)
params_full=[3.15*1e-7/V 2.1 11000 c_vec(i)]; % b , delta, p, c 
%params_full=[0.01 0.002 0.0002 0.002]; % b , delta, p, c 

b = params_full(1);
delta = params_full(2);
p = params_full(3);
c = params_full(4);
pp=p/c;
pprime(i)=pp;
%%
% time span
tspan = linspace(0,10,2000); 
dt = tspan(2)-tspan(1);
% Initial values
Init = [100000*V 0 10000*V]; 
% function
fun = @(t,y,params) [-params(1)*y(1)*y(3);params(1)*y(1)*y(3)-params(2)*y(2);...
    params(3)*y(2)-params(4)*y(3)];
% solving ODE
[t,viral_full] = ode45(@(t,y) fun(t,y,params_full),tspan,Init);
% relative error
% plotting 
figure(1)
subplot(3,2,i)
plot(t,viral_full(:,1)+viral_full(:,2),'r','LineWidth',3);ylabel('Total'); hold on
xlabel('Time')
ylabel('Total Cells')
title(num2str(delta/c))
figure(2)
subplot(3,2,i)
plot(t,viral_full(:,2),'r','LineWidth',3); hold on
xlabel('Time')
ylabel('Infected Cells')
ylim([1e-3 inf])
title(num2str(delta/c))
figure(3)
subplot(3,2,i)
plot(t,viral_full(:,3),'r','LineWidth',3);hold on
xlabel('Time')
ylabel('Virus')
title(num2str(delta/c))
ylim([1e-3 inf])
%% QSSA 
params_qssa = [b*pp delta];
beta =b ;
        epsilon = 0.001;

Initials = Init;

Init_qssa = [Initials(1) beta*Initials(1)*Initials(3)*exp((beta*Initials(3)*exp(-1)-c)/c)/c]; 

% function
fun_qssa = @(t,y,params) [-params(1)*y(1)*y(2);params(1)*y(1)*y(2)-params(2)*y(2) ];
% solving ODE
[t,vi_qssa] =ode45(@(t,y) fun_qssa(t,y,params_qssa),tspan,Init_qssa);
v=pp*vi_qssa(:,2);
% relative error
idx_vec = find(abs(tspan-1/c)<dt); 
idx = idx_vec(1);

m1= sqrt(squeeze(sum((viral_full(idx:end,1)+viral_full(idx:end,2)-vi_qssa(idx:end,1)...
    -vi_qssa(idx:end,2)).^2,1)));
m2= sqrt(squeeze(sum((viral_full(idx:end,1)+viral_full(idx:end,2)).^2,1)));
m3= sqrt(squeeze(sum((viral_full(idx:end,2)-vi_qssa(idx:end,2)).^2,1)));
m4= sqrt(squeeze(sum((viral_full(idx:end,2)).^2,1)));
m5= sqrt(squeeze(sum((viral_full(idx:end,3)-v(idx:end)).^2,1)));
m6= sqrt(squeeze(sum((viral_full(idx:end,3)).^2,1)));

rel.target_q = m1./m2;
rel.infect_q = m3./m4;
rel.viral_q = m5./m6;

rel.target_q_vec(i) = m1./m2;
rel.infect_q_vec(i) = m3./m4;
rel.viral_q_vec(i) = m5./m6;
% plotting 

figure (1)
subplot(3,2,i)
plot(t,vi_qssa(:,1)+vi_qssa(:,2),'g--','LineWidth',3);ylabel('Total'); hold on;
figure(2)
subplot(3,2,i)
plot(t,vi_qssa(:,2),'g--','LineWidth',3); hold on;
figure(3)
subplot(3,2,i)
plot(t,p/c*vi_qssa(:,2),'g--','LineWidth',3);hold on;

%% QSSA_il 
params_qssa_il = [b p*b/c delta];
Init_qssa_il = [Init(1) Init(3)]; 
% function
fun_qssa_il = @(t,y,params) [-params(1)*y(1)*y(2);params(2)*y(1)*y(2)-params(3)*y(2) ];
% solving ODE
[t,vi_qssa_il] =ode45(@(t,y) fun_qssa_il(t,y,params_qssa_il),tspan,Init_qssa_il);
I=c/p*vi_qssa_il(:,2);
% relative error
m11= sqrt(squeeze(sum((viral_full(idx:end,1)+viral_full(idx:end,2)...
    -vi_qssa_il(idx:end,1)-I(idx:end)).^2,1)));
m22= sqrt(squeeze(sum((viral_full(idx:end,1)+viral_full(idx:end,2)).^2,1)));
m33= sqrt(squeeze(sum((viral_full(idx:end,3)-vi_qssa_il(idx:end,2)).^2,1)));
m44= sqrt(squeeze(sum((viral_full(idx:end,3)).^2,1)));
m55= sqrt(squeeze(sum((viral_full(idx:end,2)-I(idx:end)).^2,1)));
m66= sqrt(squeeze(sum((viral_full(idx:end,2)).^2,1)));

rel.target_qi = m11./m22;
rel.viral_qi = m33./m44;
rel.infect_qi = m55./m66;

rel.target_qi_vec(i) = m11./m22;
rel.viral_qi_vec(i) = m33./m44;
rel.infect_qi_vec(i) = m55./m66;
% plotting 
figure (1)
subplot(3,2,i)
plot(t,vi_qssa_il(:,1)+I,'b:','LineWidth',3);hold on
legend('full',['QSSA, Rel.err= ',num2str(rel.target_q)] ,['QSSA_{il}, Rel.err= ',num2str(rel.target_qi)],'FontSize',8)
ylabel('Total')
figure(2)
subplot(3,2,i)
plot(t,c/p*vi_qssa_il(:,2),'b:','LineWidth',3)
legend('full',['QSSA, Rel.err= ',num2str(rel.infect_q)],['QSSA_{il}, Rel.err= ',num2str(rel.infect_qi)],'FontSize',8)

figure(3)
subplot(3,2,i)
plot(t,vi_qssa_il(:,2),'b:','LineWidth',3); hold on;

legend('full',['QSSA, Rel.err= ',num2str(rel.viral_q)],['QSSA_{il}, Rel.err= ',num2str(rel.viral_qi)],'FontSize',8)
end
figure(4)
rel.target_q_vec(i) = m1./m2;
rel.infect_q_vec(i) = m3./m4;
rel.viral_q_vec(i) = m5./m6;
rel.target_qi_vec(i) = m11./m22;
rel.viral_qi_vec(i) = m33./m44;
rel.infect_qi_vec(i) = m55./m66;
subplot(1,3,1)
aaa = delta./c_vec;
loglog(aaa,rel.target_q_vec,'ro-.',aaa,rel.target_qi_vec,'bx:');
xlabel('C_v')
ylabel('Rel.err (Total cells)')
legend('QSSA','QSSA_{il}')
subplot(1,3,2)
loglog(aaa,rel.infect_q_vec,'ro-.',aaa,rel.infect_qi_vec,'bx:');
xlabel('C_v')
ylabel('Rel.err (Infected cells)')
legend('QSSA','QSSA_{il}')
subplot(1,3,3)
loglog(aaa,rel.viral_q_vec,'ro-.',aaa,rel.viral_qi_vec,'bx:');
xlabel('C_v')
ylabel('Rel.err (Virus)')
legend('QSSA','QSSA_{il}')

%%
figure(5)
semilogy(t,viral_full(:,3),'r','LineWidth',3);hold on
semilogy(t,vi_qssa_il(:,2),'b:','LineWidth',3); 
legend('Viral model', 'QSSA_{il}');
xlabel('Time');ylabel('Viral load'); 
title('Invild QSSA (\delta/c = 0.035)')
xlim([0 5])

figure(6)
semilogy(t,viral_full(:,3),'r','LineWidth',3);hold on
semilogy(t,vi_qssa_il(:,2),'b:','LineWidth',3); 
semilogy(t,p/c*vi_qssa(:,2),'g--','LineWidth',3);

legend('Viral model', 'QSSA_{il}', 'QSSA');
xlabel('Time');ylabel('Viral load'); 
title('Corrected QSSA (\delta/c = 0.035)')
xlim([0 5])


figure(7)
loglog(aaa,rel.viral_q_vec,'ro-.',aaa,rel.viral_qi_vec,'bx:');
xlabel('\delta/c')
ylabel('Rel.err (Virus)')
legend('QSSA','QSSA_{il}')