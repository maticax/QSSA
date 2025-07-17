%% Parameter estimation 
% take some time!!!!!!!!!!!!!!!
clc
clear all
close all
set(0,'DefaultAxesFontSize',16)
 set(0,'DefaultAxesFontWeight','bold');
set(0, 'DefaultLineLineWidth', 3);
% viral dynamics 


% time scale & Initials

tspan = linspace(0,10,20000); 
t=tspan;
dt = tspan(2)-tspan(1);
T0=100000; %copies/ml
I0=0;
v0=10000;
V = 1; %5kg*70ml/kg=350ml
Initials=[T0*V I0*V v0*V];
beta = 3.15e-07/V;
p = 11000;
c = 60; 
k=6; % eclipse
delta = 2.1;



T=zeros(length(t),1);
I=zeros(length(t),1);
I2= zeros(length(t),1);
v=zeros(length(t),1);
T(1)=Initials(1);
I(1)=Initials(2);
v(1)=Initials(3);
%%
for i=1:length(t)-1
   
T(i+1)=T(i)-beta*T(i)*v(i)*dt;
I(i+1)=I(i)+(beta*T(i)*v(i)-delta*I(i))*dt;
v(i+1)=v(i)+(p*I(i)-c*v(i))*dt;
T1 = [T(i+1)>0];
I1 = [I(i+1)>0];
I3 = [I2(i+1)>0];
v1 = [v(i+1)>0];
T(i+1) = T(i+1)*T1;
I(i+1) = I(i+1)*I1;
v(i+1) = v(i+1)*v1;
end
ITOT =I;
CTOT = T+I;
vTOT = v;
fulldata = [t' T+I t' v]; 
rng(6)
m=10;
partial_data = (datasample(fulldata(1:15000,:),6));
%%


figure(1)
semilogy(partial_data(:,1),partial_data(:,2),'o',partial_data(:,3),partial_data(:,4),'x')
xlabel('Postinfection (day)')
legend('T_{data}+I_{data}','v_{data}')



%% initial guess and simulation
T=zeros(length(t),1);
I=zeros(length(t),1);
I2= zeros(length(t),1);
v=zeros(length(t),1);
params1(1) = 3e-7;% beta 3.15e-7; %v(0)*V
params1(2) = k;%6; %k 

params1(3) = 1.9 ;%2.1; % delta 
params1(4) =10000;%11000; % p
params1(5) = c;%500; %p'=p/c
%%
Initials=[T0*V,I0*V,v0*V];
T(1)=Initials(1);
I(1)=Initials(2);
v(1)=Initials(3);
% full model
for i=1:length(t)-1
T(i+1)=T(i)-params1(1)*T(i)*v(i)*dt;
I(i+1)=I(i)+(params1(1)*T(i)*v(i)-params1(3)*I(i))*dt;
%I2(i+1)=I2(i)+(params1(2)*I(i)-params1(3)*I2(i))*dt;
v(i+1)=v(i)+(params1(4)*I(i)-c*v(i))*dt;
T1 = [T(i+1)>0];I1 = [I(i+1)>0];I3 = [I2(i+1)>0];v1 = [v(i+1)>0];
T(i+1) = T(i+1)*T1;I(i+1) = I(i+1)*I1;I2(i+1) = I2(i+1)*I3;v(i+1) = v(i+1)*v1;
end
% qssa model 
Tr=zeros(length(t),1);Ir=zeros(length(t),1);Ir2 = zeros(length(t),1);vr=zeros(length(t),1);
Tr(1)=Initials(1);


Ir(1) = params1(1)*Initials(1)*Initials(3)*exp((params1(1)*Initials(3)*exp(-1)-params1(1)*Initials(3)-params1(3))/c)/c;
for i=1:length(t)-1    
Tr(i+1)=Tr(i)-(params1(1)*params1(4)/params1(5).*Tr(i).*Ir(i))*dt ;
Ir(i+1)=Ir(i)+(params1(1)*params1(4)/params1(5).*Tr(i).*Ir(i)-params1(3)*Ir(i))*dt;
%Ir2(i+1)=Ir2(i)+(params1(2)*Ir(i)-params1(3).*Ir2(i))*dt;
Tr1 = [Tr(i+1)>0];Ir1 = [Ir(i+1)>0];Ir3 = [Ir2(i+1)>0];
Tr(i+1) = Tr(i+1).*Tr1;Ir(i+1) = Ir(i+1).*Ir1;Ir2(i+1) = Ir2(i+1).*Ir3;

end
vr=params1(4)/params1(5)*Ir;
%% qssa_il model
% qssa model 
Tr_il=zeros(length(t),1);Ir_il=zeros(length(t),1);Ir2_il = zeros(length(t),1);vr_il=zeros(length(t),1);
Tr_il(1)=Initials(1);
vr_il(1)=Initials(3);

for i=1:length(t)-1    
Tr_il(i+1)=Tr_il(i)-(params1(1).*Tr_il(i).*vr_il(i))*dt ;
%Ir_il(i+1)=Ir_il(i)+(params1(1).*Tr_il(i).*vr_il(i)-params1(2)*Ir_il(i))*dt;
vr_il(i+1)=vr_il(i)+(params1(1)*params1(4)/params1(5).*Tr_il(i).*vr_il(i)-params1(3).*vr_il(i))*dt;
Tr1 = [Tr_il(i+1)>0];Ir1 = [Ir_il(i+1)>0];vr3 = [vr_il(i+1)>0];
Tr_il(i+1) = Tr_il(i+1).*Tr1;Ir_il(i+1) = Ir_il(i+1).*Ir1;vr_il(i+1) = vr_il(i+1).*vr3;

end
Ir2_il=vr_il*params1(5)/params1(4);
%%
tspan = t;
idx_vec = find(abs(tspan-1/c)<dt); 
idx = idx_vec(1);

m1= sqrt(squeeze(sum((T(idx:end)+I(idx:end)+I2(idx:end)-CTOT(idx:end)).^2,1)));
m2= sqrt(squeeze(sum((CTOT(idx:end)).^2,1)));
m3= sqrt(squeeze(sum((I(idx:end)+I2(idx:end)-ITOT(idx:end)).^2,1)));
m4= sqrt(squeeze(sum((ITOT(idx:end)).^2,1)));
m5= sqrt(squeeze(sum((v(idx:end)-vTOT(idx:end)).^2,1)));
m6= sqrt(squeeze(sum((vTOT(idx:end)).^2,1)));

m11= sqrt(squeeze(sum((CTOT(idx:end)-Tr(idx:end)-Ir(idx:end)-Ir2(idx:end)).^2,1)));
m22= sqrt(squeeze(sum((CTOT(idx:end)).^2,1)));
m33= sqrt(squeeze(sum((ITOT(idx:end)-Ir(idx:end)-Ir2(idx:end)).^2,1)));
m44= sqrt(squeeze(sum((ITOT(idx:end)).^2,1)));
m55= sqrt(squeeze(sum((vTOT(idx:end)-vr(idx:end)).^2,1)));
m66= sqrt(squeeze(sum((vTOT(idx:end)).^2,1)));

m111= sqrt(squeeze(sum((CTOT(idx:end)-Tr_il(idx:end)-Ir_il(idx:end)-Ir2_il(idx:end)).^2,1)));
m222= sqrt(squeeze(sum((CTOT(idx:end)).^2,1)));
m333= sqrt(squeeze(sum((ITOT(idx:end)-Ir_il(idx:end)-Ir2_il(idx:end)).^2,1)));
m444= sqrt(squeeze(sum((ITOT(idx:end)).^2,1)));
m555= sqrt(squeeze(sum((vTOT(idx:end)-vr_il(idx:end)).^2,1)));
m666= sqrt(squeeze(sum((vTOT(idx:end)).^2,1)));

rel.target_f = round(m1./m2,2);
rel.infect_f = round(m3./m4,2);
rel.viral_f = round(m5./m6,2);

rel.target_q = round(m11./m22,2);
rel.infect_q = round(m33./m44,2);
rel.viral_q = round(m55./m66,2);

rel.target_q_il = round(m111./m222,2);
rel.infect_q_il = round(m333./m444,2);
rel.viral_q_il = round(m555./m666,2);

%%
% model check
%  figure(2) 
% 
% semilogy(t,CTOT,t,T+I+I2,'r-.',t,Tr+Ir+Ir2,'b:');hold on
% scatter(partial_data(:,1),partial_data(:,2),100,'filled');hold on 
% xlabel('Postinfection (day)')
% ylabel('Total Cells')
% legend('True','full (starting)','QSSA (starting)','Data')
% title('Before parameter estimate')
% str ={'True:, \beta=1.26e-9, k=6, \delta=2.1, p=11000'};
% text(t(10),vr(5),str);
% str1 ={['Starting (f):, \beta=', num2str(params1(1)), ' k=', num2str(params1(2)),...
%     ' \delta=', num2str(params1(3)), ' p=', num2str(params1(4))]};
% text(t(10),vr(7),str1);
% text(t(10),vr(8),['Rel.err (Q)= ', num2str(rel.target_q)])
% text(t(10),vr(9),['Rel.err (f)= ', num2str(rel.target_f)])

     %text(L_s_tot(i,1),Rel_m_tot(i,1),num2str(i),'Color','red');
figure(3)

 semilogy(t,v,'r-.',t,vr,'b:');hold on
 ylim([0.1 inf])
scatter(partial_data(:,3),partial_data(:,4),100,'filled')
text(t(10),vr(3),['Rel.err (Q)= ', num2str(rel.viral_q)],'FontSize',12)
text(t(10),vr(2),['Rel.err (Q_{il})= ', num2str(rel.viral_q_il)],'FontSize',12)
text(t(10),vr(4),['Rel.err (f)= ', num2str(rel.viral_f)],'FontSize',12)
str ={['True:, \beta=3.15e-7, \delta=2.1, p=11000, c=', num2str(c)]};
text(t(10),vr(5),str,'FontSize',12);
str1 ={['Starting:, \beta=', num2str(params1(1)), '\in[3e-7, 3.5e-7]', ...
    ', \delta=', num2str(params1(3)), '\in[1.8,3], p=', num2str(params1(4)), '\in[5000,15000], c=', num2str(c)]};
text(t(10),vr(7),str1,'FontSize',12);

xlabel('Postinfection (day)')
ylabel('V')
legend('Basic viral','QSSA','Data')
title('Before parameter estimate')
%% lsqnonlin estimation 
lower =[3e-7      0  1.8      5000   params1(5)]; %beta,k,delta,p,c
upper =[3.5e-7    0  3        15000  params1(5)];
lower1 =[3e-7     0  1.8      5000   params1(5) ]; 
upper1 =[3.5e-7   0  3        15000  params1(5)];
lower2 =[3e-7     0  1.8       5000   params1(5)];
upper2 =[3.5e-7   0  3         15000  params1(5)];
x0 = [params1(1) 0 3 12000 params1(5)]; % beta,k,delta,p,p'
y0 = params1;
y0 = x0;
z0 = params1;



mmyfun = @(x)myfun(x,partial_data);
mmyfun1 =@(x)myfun1(x,partial_data);
mmyfun2 =@(x)myfun2(x,partial_data);

x=lsqnonlin(mmyfun,x0,lower,upper) %beta,k,delta,p,p' and c=40 fixed
y=lsqnonlin(mmyfun1,y0,lower1,upper1) %c=10
z=lsqnonlin(mmyfun2,z0,lower2,upper2) %c=10
%
opts  = optimoptions(@lsqnonlin, ...
          'Display','iter', ...
          'TolX',1e-9, ...
          'MaxFunctionEvaluations',5e3, ...
          'OutputFcn',@outfun);  % 

%
global trace_full
trace_full = [];

% 
fun = @(par) myfun(par, partial_data);

% lsqnonlin 
[x_est,~,~,~,~,~,~] = lsqnonlin(fun, x0, lower, upper, opts);

%%



% --- 
% trace_qssa 
global trace_qssa
trace_qssa = [];

opts_q = optimoptions(opts,'OutputFcn',@outfun_q);

% 잔차 함수
fun_q = @(par) myfun1(par, partial_data);
  
[y_est,~,~,~,~,~,~] = lsqnonlin(fun_q, x0, lower1, upper1, opts_q);
iters = 1:size(trace_full,1);

figure;
subplot(1,3,1)
semilogx(iters, trace_full(:,1), iters, trace_qssa(:,1))
xlabel('Iteration');
legend('Basic Viral','QSSA','Location','best');
title('\beta')
ylim([2.8e-7 3.2e-7])
subplot(1,3,2)
semilogx(iters, trace_full(:,3), iters, trace_qssa(:,3))
xlabel('Iteration');
legend('Basic Viral','QSSA','Location','best');
title('\delta')
ylim([2 3])
subplot(1,3,3)
semilogx(iters, trace_full(:,4), iters, trace_qssa(:,4))
xlabel('Iteration');
legend('Basic Viral','QSSA','Location','best');
title('p')
ylim([10000 13000])


%% figure after estimation
% full model
T=zeros(length(t),1);I=zeros(length(t),1);I2 = zeros(length(t),1);v=zeros(length(t),1);
T(1)= Initials(1);
v(1) = Initials(3);
for i=1:length(t)-1
T(i+1)=T(i)-(x(1)*T(i)*v(i))*dt;
I(i+1)=I(i)+(x(1)*T(i)*v(i)-x(3)*I(i))*dt;
%I2(i+1)=I2(i)+(x(2)*I(i)-x(3)*I2(i))*dt;
v(i+1)=v(i)+(x(4)*I(i)-x(5)*v(i))*dt;
T1 = [T(i+1)>0];I1 = [I(i+1)>0];I3 = [I2(i+1)>0];v1 = [v(i+1)>0];
T(i+1) = T(i+1)*T1;I(i+1) = I(i+1)*I1;I2(i+1) = I2(i+1)*I3;v(i+1) = v(i+1)*v1;
end
% qssa model 
Tr=zeros(length(t),1);Ir=zeros(length(t),1);Ir2 = zeros(length(t),1);vr=zeros(length(t),1);
Tr(1)=Initials(1);
vr(1) = Initials(3);
% Ir(1) = y(1)*Initials(1)*Initials(3)*exp(-1)/(y(1)*Initials(3)*exp(-1)+y(2)).*...
%     (1-exp(-(y(1)*Initials(3)*exp(-1)+y(2))./c));
%Ir(1) = y(1)*Initials(1)*Initials(3)*exp(-1)/c;
Ir(1) = y(1)*Initials(1)*Initials(3)*exp((y(1)*Initials(3)*exp(-1)-y(1)*Initials(3)-y(3))/c)/c;
for i=1:length(t)-1    
Tr(i+1)=Tr(i)-(y(1)*y(4)/y(5).*Tr(i).*Ir(i))*dt ;
Ir(i+1)=Ir(i)+(y(1)*y(4)/y(5).*Tr(i).*Ir(i)-y(3)*Ir(i))*dt;
%Ir2(i+1)=Ir2(i)+(y(2)*Ir(i)-y(3).*Ir2(i))*dt;
Tr1 = [Tr(i+1)>0];Ir1= [Ir(i+1)>0];Ir3 = [Ir2(i+1)>0];
Tr(i+1) = Tr(i+1).*Tr1;Ir(i+1) = Ir(i+1).*Ir1;Ir2(i+1) = Ir2(i+1).*Ir3;
end
vr=y(4)/y(5)*Ir;

%% qssa_il model
% qssa model 
Tr_il=zeros(length(t),1);Ir_il=zeros(length(t),1);Ir2_il = zeros(length(t),1);vr_il=zeros(length(t),1);
Tr_il(1)=Initials(1);
vr_il(1)=Initials(3);

for i=1:length(t)-1    
Tr_il(i+1)=Tr_il(i)-(z(1).*Tr_il(i).*vr_il(i))*dt ;
%Ir_il(i+1)=Ir_il(i)+(z(1).*Tr_il(i).*vr_il(i)-z(2)*Ir_il(i))*dt;
vr_il(i+1)=vr_il(i)+(z(1)*z(4)/z(5)*Tr_il(i).*vr_il(i)-z(3).*vr_il(i))*dt;
Tr1 = [Tr_il(i+1)>0];Ir1 = [Ir_il(i+1)>0];vr3 = [vr_il(i+1)>0];
Tr_il(i+1) = Tr_il(i+1).*Tr1;Ir_il(i+1) = Ir_il(i+1).*Ir1;vr_il(i+1) = vr_il(i+1).*vr3;

end
Ir2_il=vr_il*z(5)/z(4);
%%
m1= sqrt(squeeze(sum((T(idx:end)+I(idx:end)+I2(idx:end)-CTOT(idx:end)).^2,1)));
m2= sqrt(squeeze(sum((CTOT(idx:end)).^2,1)));
m3= sqrt(squeeze(sum((I(idx:end)+I2(idx:end)-ITOT(idx:end)).^2,1)));
m4= sqrt(squeeze(sum((ITOT(idx:end)).^2,1)));
m5= sqrt(squeeze(sum((v(idx:end)-vTOT(idx:end)).^2,1)));
m6= sqrt(squeeze(sum((vTOT(idx:end)).^2,1)));

m11= sqrt(squeeze(sum((CTOT(idx:end)-Tr(idx:end)-Ir(idx:end)-Ir2(idx:end)).^2,1)));
m22= sqrt(squeeze(sum((CTOT(idx:end)).^2,1)));
m33= sqrt(squeeze(sum((ITOT(idx:end)-Ir(idx:end)-Ir2(idx:end)).^2,1)));
m44= sqrt(squeeze(sum((ITOT(idx:end)).^2,1)));
m55= sqrt(squeeze(sum((vTOT(idx:end)-vr(idx:end)).^2,1)));
m66= sqrt(squeeze(sum((vTOT(idx:end)).^2,1)));

m111= sqrt(squeeze(sum((CTOT(idx:end)-Tr_il(idx:end)-Ir_il(idx:end)-Ir2_il(idx:end)).^2,1)));
m222= sqrt(squeeze(sum((CTOT(idx:end)).^2,1)));
m333= sqrt(squeeze(sum((ITOT(idx:end)-Ir_il(idx:end)-Ir2_il(idx:end)).^2,1)));
m444= sqrt(squeeze(sum((ITOT(idx:end)).^2,1)));
m555= sqrt(squeeze(sum((vTOT(idx:end)-vr_il(idx:end)).^2,1)));
m666= sqrt(squeeze(sum((vTOT(idx:end)).^2,1)));

rel.target_f = round(m1./m2,2);
rel.infect_f = round(m3./m4,2);
rel.viral_f = round(m5./m6,2);

rel.target_q = round(m11./m22,2);
rel.infect_q = round(m33./m44,2);
rel.viral_q = round(m55./m66,2);

rel.target_q_il = round(m111./m222,2);
rel.infect_q_il = round(m333./m444,2);
rel.viral_q_il = round(m555./m666,2);
%%


figure(4)
semilogy (t,T+I,'r-.',t,Tr+Ir,'b:',t,Tr_il+Ir2_il,'g-.');hold on
scatter(partial_data(:,1),partial_data(:,2),100,'filled');

xlabel('Postinfection (Day)')
ylabel('Total Cells')
legend('full','QSSA','QSSA_{il}','Data')
title('After parameter estimation')
str ={['Estimate (f):, \beta=', num2str(x(1)), ' k=', num2str(x(2)),...
    ' \delta=', num2str(x(3)), ' p=', num2str(x(4))]};
% text(t(10),vr(5),str);
% str1 ={['Estimate (Q):, \beta=', num2str(y(1)), ' k=', num2str(y(2)),...
%      ' \delta=', num2str(y(3)), ' p=', num2str(y(5))]};
% text(t(10),vr(10),str1);
%   text(t(10),vr(8),['Rel.err (Q)= ', num2str(rel.target_q)])
% text(t(10),vr(9),['Rel.err (f)= ', num2str(rel.target_f)])
 figure(5)
semilogy(t,v,'r-.',t,vr,'b:'); hold on
 scatter(partial_data(:,3),partial_data(:,4),100,'filled');
 xlabel('Postinfection (day)')
 ylabel('V')

legend('Basic viral','QSSA','QSSA_{il}','Data')
title('After parameter estimation')
str ={['Estimate (f):, \beta=', num2str(x(1)), ...
    ', \delta=', num2str(x(3)), ', p=', num2str(x(4))]};
text(t(10),vr(5),str,'FontSize',12);
str1 ={['Estimate (Q):, \beta=', num2str(y(1)), ...
     ', \delta=', num2str(y(3)), ', p=', num2str(y(4))]};
text(t(10),vr(10),str1,'FontSize',12);
% str2 ={['Estimate (Q_{il}):, \beta=', num2str(z(1)), ...
%      ', \delta=', num2str(z(3)), ', p=', num2str(z(4))]};
%text(t(10),vr(13),str2,'FontSize',12);
%text(t(10),vr(2),['Rel.err (Q_{il})= ', num2str(rel.viral_q_il)],'FontSize',12)
text(t(10),vr(3),['Rel.err (Q)= ', num2str(rel.viral_q)],'FontSize',12)
text(t(10),vr(4),['Rel.err (f)= ', num2str(rel.viral_f)],'FontSize',12)

ylim([0.1 inf])

%%

% %% --- 공통 세팅 ---
% % (이미 lsqnonlin 으로 구해놓은) 초기 추정치
% beta0_viral = x(1);   % viral 모델에서 추정된 β
% beta0_qssa  = y(1);   % QSSA 모델에서 추정된 β
% delta0      = 2.1;    % 초기 δ
% p0          = 1e4;    % 초기 p
% c_fixed     = 60;     % 고정된 c
% 
% % 최적화 옵션
% opts = optimoptions(@lsqnonlin, ...
%     'Display','off', ...
%     'TolX',1e-8, ...
%     'MaxFunctionEvaluations',5000);
% 
% % 프로파일용 격자
% beta_grid  = linspace(1e-8,5e-7,30);   
% delta_grid = linspace(1.5,3,30);      
% p_grid     = linspace(10000,13000,30);  
% 
% % 로그‐변환용 함수 핸들 (xp = [ log(delta); log(p) ])
% viral_obj = @(b, xp, data) myfun( [ b, 0, exp(xp(1)), exp(xp(2)), c_fixed ], data );
% qssa_obj  = @(b, xp, data) myfun1([ b, 0, exp(xp(1)), exp(xp(2)), c_fixed ], data );
% 
% % 1) β‐프로파일 (δ,p warm‐start + 로그 re‐parameterization)
% SSE_beta_viral = nan(size(beta_grid));
% SSE_beta_qssa  = nan(size(beta_grid));
% 
% % warm‐start 변수 (xp = [ log(delta); log(p) ])
% last_xp_v = [ log(delta0); log(p0) ];
% last_xp_q = [ log(delta0); log(p0) ];
% 
% for i = 1:numel(beta_grid)
%     b = beta_grid(i);
% 
%     % Viral model: xp = [ log(delta); log(p) ] 최적화
%     obj_v = @(xp) viral_obj(b, xp, partial_data);
%     lb_xp = [ log(1.0);   log(5e3) ];   % δ ≥ 1, p ≥ 5e3
%     ub_xp = [ log(4.0);   log(2e4) ];   % δ ≤ 4, p ≤ 2e4
%     [xp_v, res_v] = lsqnonlin(obj_v, last_xp_v, lb_xp, ub_xp, opts);
%     SSE_beta_viral(i) = sum(res_v.^2);
%     last_xp_v = xp_v;   % warm‐start 업데이트
% 
%     % QSSA model: xp 최적화
%     obj_q = @(xp) qssa_obj(b, xp, partial_data);
%     [xp_q, res_q] = lsqnonlin(obj_q, last_xp_q, lb_xp, ub_xp, opts);
%     SSE_beta_qssa(i) = sum(res_q.^2);
%     last_xp_q = xp_q;   % warm‐start 업데이트
% end
% 
% figure; hold on;
% plot(beta_grid, log(SSE_beta_viral), 'r-o','LineWidth',1.5);
% plot(beta_grid, log(SSE_beta_qssa),  'b--s','LineWidth',1.5);
% xlabel('\beta'); ylabel('log(SSE)');
% legend('Viral','QSSA','Location','best');
% title('Profile‐Likelihood in \beta');
% grid on;
% 
% 
% % 2) δ‐프로파일 (β,p 재추정 → warm‐start & 로그 re‐parameterization)
% SSE_delta_viral = nan(size(delta_grid));
% SSE_delta_qssa  = nan(size(delta_grid));
% 
% % warm‐start 변수 (xp = [ log(beta); log(p) ])
% last_xp_v = [ log(beta0_viral); log(p0) ];
% last_xp_q = [ log(beta0_qssa);  log(p0) ];
% 
% viral_obj2 = @(d, xp, data) myfun( [ exp(xp(1)), 0, d, exp(xp(2)), c_fixed ], data );
% qssa_obj2  = @(d, xp, data) myfun1([ exp(xp(1)), 0, d, exp(xp(2)), c_fixed ], data );
% 
% for i = 1:numel(delta_grid)
%     d = delta_grid(i);
% 
%     % Viral model: xp = [ log(beta); log(p) ]
%     obj_v = @(xp) viral_obj2(d, xp, partial_data);
%     lb_xp = [ log(1e-7); log(5e3) ];   % β ≥ 1e-7, p ≥ 5e3
%     ub_xp = [ log(5e-7); log(2e4) ];   % β ≤ 5e-7, p ≤ 2e4
%     [xp_v, res_v] = lsqnonlin(obj_v, last_xp_v, lb_xp, ub_xp, opts);
%     SSE_delta_viral(i) = sum(res_v.^2);
%     last_xp_v = xp_v;
% 
%     % QSSA model
%     obj_q = @(xp) qssa_obj2(d, xp, partial_data);
%     [xp_q, res_q] = lsqnonlin(obj_q, last_xp_q, lb_xp, ub_xp, opts);
%     SSE_delta_qssa(i) = sum(res_q.^2);
%     last_xp_q = xp_q;
% end
% 
% figure; hold on;
% plot(delta_grid, log(SSE_delta_viral), 'r-o','LineWidth',1.5);
% plot(delta_grid, log(SSE_delta_qssa),  'b--s','LineWidth',1.5);
% xlabel('\delta'); ylabel('log(SSE)');
% legend('Viral','QSSA','Location','best');
% title('Profile‐Likelihood in \delta');
% grid on;
% 
% 
% % 3) p‐ (β,δ 재추정 → warm‐start & 로그 re‐parameterization)
% SSE_p_viral = nan(size(p_grid));
% SSE_p_qssa  = nan(size(p_grid));
% 
% % warm‐start  (xp = [ log(beta); log(delta) ])
% last_xp_v = [ log(beta0_viral); log(delta0) ];
% last_xp_q = [ log(beta0_qssa);  log(delta0) ];
% 
% viral_obj3 = @(pval, xp, data) myfun( [ exp(xp(1)), 0, exp(xp(2)), pval, c_fixed ], data );
% qssa_obj3  = @(pval, xp, data) myfun1([ exp(xp(1)), 0, exp(xp(2)), pval, c_fixed ], data );
% 
% for i = 1:numel(p_grid)
%     pval = p_grid(i);
% 
%     % Viral model
%     obj_v = @(xp) viral_obj3(pval, xp, partial_data);
%     lb_xp = [ log(1e-7); log(1.0) ];    % β ≥ 1e-7, δ ≥ 1
%     ub_xp = [ log(5e-7); log(4.0) ];    % β ≤ 5e-7, δ ≤ 4
%     [xp_v, res_v] = lsqnonlin(obj_v, last_xp_v, lb_xp, ub_xp, opts);
%     SSE_p_viral(i) = sum(res_v.^2);
%     last_xp_v = xp_v;
% 
%     % QSSA model
%     obj_q = @(xp) qssa_obj3(pval, xp, partial_data);
%     [xp_q, res_q] = lsqnonlin(obj_q, last_xp_q, lb_xp, ub_xp, opts);
%     SSE_p_qssa(i) = sum(res_q.^2);
%     last_xp_q = xp_q;
% end
% 
% figure; hold on;
% plot(p_grid, log(SSE_p_viral), 'r-o','LineWidth',1.5);
% plot(p_grid, log(SSE_p_qssa),  'b--s','LineWidth',1.5);
% xlabel('p'); ylabel('log(SSE)');
% legend('Viral','QSSA','Location','best');
% title('Profile‐Likelihood in p');
% grid on;


%%
%% ------------------------------------------------------------------------
% 3×2 2D profile-likelihood  (Viral vs QSSA )
%   (1) β–δ plane (p )
%   (2) β–p plane (δ )
%   (3) δ–p plane (β )
%   plane  ΔSSE ≤ th (≈95% CI) contour 

% --- 0)  --------------------------------------------------------
%  (lsqnonlin)
beta0_v   = x(1);      % Vira β
beta0_q   = y(1);      % QSSA β
delta0_v    = 2.6104;       % δ
delta0_q = 2.0999;
p0_v        = 11934;       %  p
p0_q = 10884;
c_fixed   = 60;        % c

% lsqnonlin 
opts = optimoptions(@lsqnonlin, ...
    'Display','off', ...
    'TolX',1e-8, ...
    'MaxFunctionEvaluations',5000);
%lower =[3e-7      0  1.8      5000   params1(5)]; %beta,k,delta,p,c
%upper =[3.5e-7    0  3        15000  params1(5)];
% 프로파일 격자
beta_grid  = linspace(2.5e-7,3.5e-7,150);
delta_grid = linspace(1.5,3,150);
p_grid     = linspace(8000,15000,150);2

% ODE residual 함수 핸들
viral_fun = @(pars,data)   myfun( pars,  data );    % pars = [β,0,δ,p,c]
qssa_fun  = @(pars,data) myfun1(pars,  data );

% ------------------------------ 1) β–δ plane -----------------------------
% p 재추정 → xp = log(p)
SSE_bd_v = nan(numel(delta_grid),numel(beta_grid));
SSE_bd_q = nan(size(SSE_bd_v));
xp_v = log(p0_v);
xp_q = log(p0_q);

for i = 1:numel(beta_grid)
  b = beta_grid(i);
  for j = 1:numel(delta_grid)
    d = delta_grid(j);

    % Viral
    objV = @(xp) viral_fun([b,0,d,exp(xp_v),c_fixed], partial_data);
    [xp_v, rV] = lsqnonlin(objV, xp_v, log(5e3), log(2e4), opts);
    SSE_bd_v(j,i) = sum(rV.^2);

    % QSSA
    objQ = @(xp) qssa_fun([b,0,d,exp(xp_q),c_fixed], partial_data);
    [xp_q, rQ] = lsqnonlin(objQ, xp_q, log(5e3), log(2e4), opts);
    SSE_bd_q(j,i) = sum(rQ.^2);
  end
end

% ------------------------------ 2) β–p plane -----------------------------
% δ 재추정 → xp = log(δ)
SSE_bp_v = nan(numel(p_grid),numel(beta_grid));
SSE_bp_q = nan(size(SSE_bp_v));
xp_v = log(delta0_v);
xp_q = log(delta0_q);

for i = 1:numel(beta_grid)
  b = beta_grid(i);
  for k = 1:numel(p_grid)
    pval = p_grid(k);

    % Viral
    objV = @(xp) viral_fun([b,0,exp(xp_v),pval,c_fixed], partial_data);
    [xp_v, rV] = lsqnonlin(objV, xp_v, log(1), log(4), opts);
    SSE_bp_v(k,i) = sum(rV.^2);

    % QSSA
    objQ = @(xp) qssa_fun([b,0,exp(xp_q),pval,c_fixed], partial_data);
    [xp_q, rQ] = lsqnonlin(objQ, xp_q, log(1), log(4), opts);
    SSE_bp_q(k,i) = sum(rQ.^2);
  end
end

% ------------------------------ 3) δ–p plane -----------------------------
% β  → xp = log(β)
SSE_dp_v = nan(numel(p_grid),numel(delta_grid));
SSE_dp_q = nan(size(SSE_dp_v));
xp_v = log(beta0_v);
xp_q = log(beta0_q);

for j = 1:numel(delta_grid)
  d = delta_grid(j);
  for k = 1:numel(p_grid)
    pval = p_grid(k);

    % Viral
    objV = @(xp) viral_fun([exp(xp_v),0,d,pval,c_fixed], partial_data);
    [xp_v, rV] = lsqnonlin(objV, xp_v, log(1e-7), log(5e-7), opts);
    SSE_dp_v(k,j) = sum(rV.^2);

    % QSSA
    objQ = @(xp) qssa_fun([exp(xp_q),0,d,pval,c_fixed], partial_data);
    [xp_q, rQ] = lsqnonlin(objQ, xp_q, log(1e-7), log(5e-7), opts);
    SSE_dp_q(k,j) = sum(rQ.^2);
  end
end

% ------------------------- 4) ΔSSE  & masking -----------------------
th = 5.99;  % χ²(2;0.95) = 5.99

% 1) ΔSSE  (95% CI )
% — β–δ
min_bd_v = min(SSE_bd_v(:));
min_bd_q = min(SSE_bd_q(:));
dSSE_bd_v = SSE_bd_v - min_bd_v;  dSSE_bd_v(dSSE_bd_v>th) = NaN;
dSSE_bd_q = SSE_bd_q - min_bd_q;  dSSE_bd_q(dSSE_bd_q>th) = NaN;

% — β–p
min_bp_v = min(SSE_bp_v(:));
min_bp_q = min(SSE_bp_q(:));
dSSE_bp_v = SSE_bp_v - min_bp_v;  dSSE_bp_v(dSSE_bp_v>th) = NaN;
dSSE_bp_q = SSE_bp_q - min_bp_q;  dSSE_bp_q(dSSE_bp_q>th) = NaN;

% — δ–p
min_dp_v = min(SSE_dp_v(:));
min_dp_q = min(SSE_dp_q(:));
dSSE_dp_v = SSE_dp_v - min_dp_v;  dSSE_dp_v(dSSE_dp_v>th) = NaN;
dSSE_dp_q = SSE_dp_q - min_dp_q;  dSSE_dp_q(dSSE_dp_q>th) = NaN;

% 2) meshgrid
[BB,  DD ] = meshgrid(beta_grid,  delta_grid);
[BB2, PP ] = meshgrid(beta_grid,  p_grid    );
[DD2, PP2] = meshgrid(delta_grid, p_grid    );

% 3) min
[~,idx_bd_v] = min(SSE_bd_v(:)); [r_bd_v,c_bd_v] = ind2sub(size(SSE_bd_v),idx_bd_v);
[~,idx_bd_q] = min(SSE_bd_q(:)); [r_bd_q,c_bd_q] = ind2sub(size(SSE_bd_q),idx_bd_q);

[~,idx_bp_v] = min(SSE_bp_v(:)); [r_bp_v,c_bp_v] = ind2sub(size(SSE_bp_v),idx_bp_v);
[~,idx_bp_q] = min(SSE_bp_q(:)); [r_bp_q,c_bp_q] = ind2sub(size(SSE_bp_q),idx_bp_q);

[~,idx_dp_v] = min(SSE_dp_v(:)); [r_dp_v,c_dp_v] = ind2sub(size(SSE_dp_v),idx_dp_v);
[~,idx_dp_q] = min(SSE_dp_q(:)); [r_dp_q,c_dp_q] = ind2sub(size(SSE_dp_q),idx_dp_q);

% 4)plot
figure('Position',[100 100 1200 900]);

% -- (1) β–δ plane --------------------------------------------------------
% Viral
ax1 = subplot(3,2,1);
pcolor(ax1, BB, DD, log10(SSE_bd_v)); shading(ax1,'interp'); hold(ax1,'on');
contour(ax1, BB, DD, dSSE_bd_v, [0 th], 'k-', 'LineWidth',2);
plot( ax1, beta_grid(c_bd_v), delta_grid(r_bd_v), 'kp', ...
      'MarkerSize',14, 'MarkerFaceColor','y');
title(ax1,'Viral: \beta–\delta SSE'); xlabel(ax1,'\beta'); ylabel(ax1,'\delta');
set(ax1,'Layer','top'); colorbar(ax1); grid(ax1,'on');

% QSSA
ax2 = subplot(3,2,2);
pcolor(ax2, BB, DD, log10(SSE_bd_q)); shading(ax2,'interp'); hold(ax2,'on');
contour(ax2, BB, DD, dSSE_bd_q, [0 th], 'k--','LineWidth',2);
plot( ax2, beta_grid(c_bd_q), delta_grid(r_bd_q), 'kp', ...
      'MarkerSize',14, 'MarkerFaceColor','y');
title(ax2,'QSSA:  \beta–\delta SSE'); xlabel(ax2,'\beta'); ylabel(ax2,'\delta');
set(ax2,'Layer','top'); colorbar(ax2); grid(ax2,'on');

% -- (2) β–p plane --------------------------------------------------------
% Viral
ax3 = subplot(3,2,3);
pcolor(ax3, BB2, PP, log10(SSE_bp_v)); shading(ax3,'interp'); hold(ax3,'on');
contour(ax3, BB2, PP, dSSE_bp_v, [0 th], 'k-', 'LineWidth',2);
plot( ax3, beta_grid(c_bp_v), p_grid(r_bp_v), 'kp', ...
      'MarkerSize',14, 'MarkerFaceColor','y');
title(ax3,'Viral: \beta–p SSE'); xlabel(ax3,'\beta'); ylabel(ax3,'p');
set(ax3,'Layer','top'); colorbar(ax3); grid(ax3,'on');

% QSSA
ax4 = subplot(3,2,4);
pcolor(ax4, BB2, PP, log10(SSE_bp_q)); shading(ax4,'interp'); hold(ax4,'on');
contour(ax4, BB2, PP, dSSE_bp_q, [0 th], 'k--','LineWidth',2);
plot( ax4, beta_grid(c_bp_q), p_grid(r_bp_q), 'kp', ...
      'MarkerSize',14, 'MarkerFaceColor','y');
title(ax4,'QSSA:  \beta–p SSE'); xlabel(ax4,'\beta'); ylabel(ax4,'p');
set(ax4,'Layer','top'); colorbar(ax4); grid(ax4,'on');

% -- (3) δ–p plane --------------------------------------------------------
% Viral
ax5 = subplot(3,2,5);
pcolor(ax5, DD2, PP2, log10(SSE_dp_v)); shading(ax5,'interp'); hold(ax5,'on');
contour(ax5, DD2, PP2, dSSE_dp_v, [0 th], 'k-', 'LineWidth',2);
plot( ax5, delta_grid(c_dp_v), p_grid(r_dp_v), 'kp', ...
      'MarkerSize',14, 'MarkerFaceColor','y');
title(ax5,'Viral: \delta–p SSE'); xlabel(ax5,'\delta'); ylabel(ax5,'p');
set(ax5,'Layer','top'); colorbar(ax5); grid(ax5,'on');

% QSSA
ax6 = subplot(3,2,6);
pcolor(ax6, DD2, PP2, log10(SSE_dp_q)); shading(ax6,'interp'); hold(ax6,'on');
contour(ax6, DD2, PP2, dSSE_dp_q, [0 th], 'k--','LineWidth',2);
plot( ax6, delta_grid(c_dp_q), p_grid(r_dp_q), 'kp', ...
      'MarkerSize',14, 'MarkerFaceColor','y');
title(ax6,'QSSA:  \delta–p SSE'); xlabel(ax6,'\delta'); ylabel(ax6,'p');
set(ax6,'Layer','top'); colorbar(ax6); grid(ax6,'on');






%%
function F = myfun(params1, partial_data) %full identifiability

tspan = linspace(0,10,20000); 
dt = tspan(2)-tspan(1);
t=tspan;

V1=300;
T_est=zeros(length(t),1);
I_est=zeros(length(t),1);
T_est(1)=30000000/V1;
v_est = zeros(length(t),1);
v_est(1) = 3000000/V1;
for i=1:length(t)-1

T_est(i+1)=T_est(i)-(params1(1)*T_est(i)*v_est(i))*dt;
I_est(i+1)=I_est(i)+(params1(1)*T_est(i)*v_est(i)-params1(3)*I_est(i))*dt;
v_est(i+1)=v_est(i)+(params1(4)*I_est(i)-params1(5)*v_est(i))*dt;
T1 = [T_est(i+1)>0];I1 = [I_est(i+1)>0];v1 = [v_est(i+1)>0];
T_est(i+1) = T_est(i+1)*T1;I_est(i+1) = I_est(i+1)*I1;v_est(i+1) = v_est(i+1)*v1;
end

for i=1:length(partial_data(:,1))
time(i) = find(abs(t-partial_data(i,3))<=0.5*dt);
end
totv = v_est;
totv = totv(time);

F=(abs(partial_data(:,4)-totv));

end


function F= myfun1(params1, partial_data) %qssa identifiability
tspan = linspace(0,10,20000); 
dt = tspan(2)-tspan(1);
t=tspan;



% qssa model 
Tr_est=zeros(length(t),1);Ir_est=zeros(length(t),1);vr_est=zeros(length(t),1);
V1=300;
Tr_est(1)=30000000/V1;
vrk = 3000000/V1;

Ir_est(1) = params1(1)*Tr_est(1)*vrk*exp((params1(1)*vrk*exp(-1)-params1(1)*vrk-params1(3))/params1(5))/params1(5);

for i=1:length(t)-1    
Tr_est(i+1)=Tr_est(i)-(params1(1)*params1(4)/params1(5)*Tr_est(i).*Ir_est(i))*dt ;
Ir_est(i+1)=Ir_est(i)+(params1(1)*params1(4)/params1(5)*Tr_est(i).*Ir_est(i)-params1(3)*Ir_est(i))*dt;
Tr1 = [Tr_est(i+1)>0];Ir1 = [Ir_est(i+1)>0];
Tr_est(i+1) = Tr_est(i+1).*Tr1;Ir_est(i+1) = Ir_est(i+1).*Ir1;
end
vr_est=params1(4)/params1(5)*Ir_est;

for i=1:length(partial_data(:,1))
time(i) = find(abs(t-partial_data(i,3))<=0.5*dt);
end

totv = vr_est;
totv = totv(time);

F=(abs(partial_data(:,4)-totv));
end

function F= myfun2(params1, partial_data) 
tspan = linspace(0,10,20000); 
dt = tspan(2)-tspan(1);
t=tspan;



% qssa_il model 
Tr_est_il=zeros(length(t),1);Ir_est_il=zeros(length(t),1);Ir2_est_il = zeros(length(t),1);vr_est_il=zeros(length(t),1);
V1=300;
Tr_est_il(1)=30000000/V1;
vr_est_il(1) = 3000000/V1;


for i=1:length(t)-1    
Tr_est_il(i+1)=Tr_est_il(i)-(params1(1)*Tr_est_il(i).*vr_est_il(i))*dt ;
vr_est_il(i+1)=vr_est_il(i)+(params1(1)*params1(4)/params1(5)*Tr_est_il(i).*vr_est_il(i)-params1(3)*vr_est_il(i))*dt;
Tr1 = [Tr_est_il(i+1)>0];Ir1 = [Ir_est_il(i+1)>0];vr3 = [vr_est_il(i+1)>0];
Tr_est_il(i+1) = Tr_est_il(i+1).*Tr1;Ir_est_il(i+1) = Ir_est_il(i+1).*Ir1;vr_est_il(i+1) = vr_est_il(i+1).*vr3;
end
Ir_est_il=vr_est_il*params1(5)/params1(4);

for i=1:length(partial_data(:,1))
time(i) = find(abs(t-partial_data(i,3))<=0.5*dt);
end

totv = vr_est_il;
totv = totv(time);

F=(abs(partial_data(:,4)-totv));
end

function stop = outfun(x,optimValues,state)
    % full model trace 
    global trace_full
    stop = false;
    if strcmp(state,'iter')
        trace_full = [trace_full; x];
    end
end

function stop = outfun_q(x,optimValues,state)
    % QSSA model trace 
    global trace_qssa
    stop = false;
    if strcmp(state,'iter')
        trace_qssa = [trace_qssa; x];
    end
end