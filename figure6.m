%% **data generation and Basic vs QSSA fit(Small & Large C_v)**
clc;
clear;
close all;
set(0, 'DefaultAxesFontSize', 16, 'DefaultAxesFontWeight', 'bold');
set(0, 'DefaultLineLineWidth', 3);
rng(11);

%
delta = 2.1;
beta_true = 3.15e-7;  
p_true = 11000; %  p 
T0 = 100000; % Target cells 
I0 = 0;
V0 = 10000; % Virus
V = 1; % 5kg * 70ml/kg = 350ml
Initials = [T0 * V, I0 * V, V0 * V];

% C_v 
c_small = delta / 60; % Small C_v (<= 0.03)
c_large = delta / 10; % Large C_v (> 0.03)
cases = {'C_v= 0.035', 'C_v= 0.21'};
Cv_vec = [c_small, c_large];
c_values = [60, 10];

% time
tspan = linspace(0, 10, 200);
sample_idx = round(linspace(1, length(tspan), 10)); % 10

% data storage
T_data = cell(2,1);
I_data = cell(2,1);
V_data = cell(2,1);
partial_data = cell(2,1);

% 
num_early = 6; %
num_late = 4;  % 

% 
sample_idx_early = round(linspace(1, round(length(tspan) * 0.3), num_early));

% 
sample_idx_late = round(linspace(round(length(tspan) * 0.3), length(tspan), num_late));

% 
sample_idx = unique([sample_idx_early, sample_idx_late]);

% 
for i = 1:2
    c = c_values(i);
    Cv = Cv_vec(i);

    params = [beta_true, c*Cv, p_true, c];

    % Basic Viral Model ODE45 
    [t, Y] = ode45(@(t,y) viral_model(t, y, params), tspan, Initials);

    % 
    num_samples = length(sample_idx); % 
    T_data{i} = Y(sample_idx,1) .* (1 + 0.05 * randn(num_samples,1)); % 5% noise
    I_data{i} = Y(sample_idx,2) .* (1 + 0.05 * randn(num_samples,1));
    V_data{i} = Y(sample_idx,3) .* (1 + 0.05 * randn(num_samples,1));
    partial_data{i} = [t(sample_idx), T_data{i}, I_data{i}, V_data{i}];
end
%% Basic Viral Model & QSSA Model fit
params_est_basic = zeros(2, 4);
params_est_qssa = zeros(2, 4);  % QSSA

for i = 1:2
    c = c_values(i);
    Cv = Cv_vec(i);

    data = partial_data{i};
    
 
    % x0 = [beta_true, delta, p_true, c];
    % lower = [3e-7, 1.8, 5000, c * 0.5];
    % upper = [3.5e-7, 3, 15000, c * 1.5];
 x0 = [beta_true, 1.2*Cv*c, 8000, c*0.8];
    lower = [beta_true, 1.8, 5000, c * 0.5];
    upper = [beta_true, 3, 15000, c * 1.5];



   fun_basic = @(x) error_basic(x, data, T0, V0, V);
params_est_basic(i, :) = lsqnonlin(fun_basic, x0, lower, upper);

% **QSSA Model fit (β, δ, p', c)**
    beta_init = beta_true;
    delta_init = delta;
    p_prime_init = p_true / c;
    c_init = c;
     x0_qssa = [beta_true, 0.95*c_init*Cv, 0.95*p_prime_init, 1.05*c_init];
    lower_qssa = [beta_true, 2, 105000/c, 0.5*c ];
    upper_qssa = [beta_true, 2.2, 11500/c, 1.5*c];

    % x0_qssa = [beta_init, delta_init, p_prime_init, c_init];
    % lower_qssa = [3.13e-7, 2, 105000/c, 0.5*c ];
    % upper_qssa = [3.16e-7, 2.2, 11500/c, 1.5*c];

  fun_qssa = @(x) error_qssa(x, data, T0, V0, V);
params_est_qssa(i, :) = lsqnonlin(fun_qssa, x0_qssa, lower_qssa, upper_qssa);
end

%% Dynamics  (Basic vs QSSA, Relative Error )
for i = 1:2
    figure;
    c = c_values(i);
    data = partial_data{i};
    params_basic = params_est_basic(i, :);
    params_qssa = params_est_qssa(i, :);

    % ODE  (Basic Model)
    [t_basic, Y_basic] = ode45(@(t,y) viral_model(t, y, params_basic), tspan, Initials);

    % ODE  (QSSA Model)
Init_qssa = [T0 * V, params_qssa(1) * T0 * V0 * exp((params_qssa(1) * V0 * exp(-1) - params_qssa(1) * V0 - params_qssa(2)) / params_qssa(3)) / params_qssa(3)];
[t_qssa, Y_qssa] = ode45(@(t,y) qssa_model(t, y, params_qssa), tspan, Init_qssa);
    V_qssa = params_qssa(3) * Y_qssa(:,2); % QSSA Virus dynamics

    %  Relative Error 
    
    T_interp_basic = interp1(t_basic, Y_basic(:,1), data(:,1));  % Target Cell
    I_interp_basic = interp1(t_basic, Y_basic(:,2), data(:,1));  % Infected Cell
    V_interp_basic = interp1(t_basic, Y_basic(:,3), data(:,1));  % Virus Load

    T_interp_qssa = interp1(t_qssa, Y_qssa(:,1), data(:,1));  % QSSA Target Cell
    I_interp_qssa = interp1(t_qssa, Y_qssa(:,2), data(:,1));  % QSSA Infected Cell
    V_interp_qssa = interp1(t_qssa, V_qssa, data(:,1));       % QSSA Virus Load

    % 
    rel_error_T_basic = norm(T_interp_basic - data(:,2)) / norm(data(:,2));
    rel_error_I_basic = norm(I_interp_basic - data(:,3)) / norm(data(:,3));
    rel_error_V_basic = norm(V_interp_basic - data(:,4)) / norm(data(:,4));

    rel_error_T_qssa = norm(T_interp_qssa - data(:,2)) / norm(data(:,2));
    rel_error_I_qssa = norm(I_interp_qssa - data(:,3)) / norm(data(:,3));
    rel_error_V_qssa = norm(V_interp_qssa - data(:,4)) / norm(data(:,4));

    %  Target Cell (T) 
    subplot(3,1,1);
    hold on;
    plot(t_basic, Y_basic(:,1), 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Basic Model (Err: %.2f)', rel_error_T_basic));
    plot(t_qssa, Y_qssa(:,1), 'r--', 'LineWidth', 2, 'DisplayName', sprintf('QSSA Model (Err: %.2f)', rel_error_T_qssa));
    scatter(data(:,1), data(:,2), 100, 'ko', 'DisplayName', 'Data');
    xlabel('Time');
    ylabel('T');
    %title(sprintf('%s: T Dynamics', cases{i}));
    legend show;
    grid on;

    %Infected Cell (I) 
    subplot(3,1,2);
    hold on;
    plot(t_basic, Y_basic(:,2), 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Basic Model (Err: %.2f)', rel_error_I_basic));
    plot(t_qssa, Y_qssa(:,2), 'r--', 'LineWidth', 2, 'DisplayName', sprintf('QSSA Model (Err: %.2f)', rel_error_I_qssa));
    scatter(data(:,1), data(:,3), 100, 'ko', 'DisplayName', 'Data');
    xlabel('Time');
    ylabel('I');
    %title(sprintf('%s: I Dynamics', cases{i}));
    legend show;
    grid on;

     
    subplot(3,1,3);
    hold on;
    plot(t_basic, Y_basic(:,3), 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Basic Model (Err: %.2f)', rel_error_V_basic));
    plot(t_qssa, V_qssa, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('QSSA Model (Err: %.2f)', rel_error_V_qssa));
    scatter(data(:,1), data(:,4), 100, 'ko', 'DisplayName', 'Data');
    xlabel('Time');
    ylabel('V');
    %title(sprintf('%s: V Dynamics', cases{i}));
    legend show;
    grid on;
    sgtitle(cases{i});
end
%% 
true_params_basic = [delta, p_true, c_values(2);  % c small
                     delta, p_true, c_values(1)]; % c large

true_params_qssa = [delta, p_true / c_values(2), c_values(2);
                    delta, p_true / c_values(1), c_values(1)];

% 
figure;
subplot(1,2,1);
bar_params = categorical({'Delta', 'p', 'c'});
bar_values_basic = [true_params_basic(2,:);  params_est_basic(1, 2:end); 
                    true_params_basic(1,:); params_est_basic(2, 2:end)];

bar(bar_params, bar_values_basic');
xlabel('Parameter');
ylabel('Estimated Values');
legend('True (Small C_v)', 'Estimated (Small C_v)', ...
       'True (Large C_v)', 'Estimated (Large C_v)');
title('Basic Viral Model');
grid on;

% 
subplot(1,2,2);
bar_params_qssa = categorical({'Delta', 'p''', 'c'});  % p -> p'
bar_values_qssa = [ true_params_qssa(2, :);params_est_qssa(1, 2:end); 
                   true_params_qssa(1,:); params_est_qssa(2, 2:end)];

bar(bar_params_qssa, bar_values_qssa');
xlabel('Parameter');
ylabel('Estimated Values');
legend('True (Small C_v)', 'Estimated (Small C_v)', ...
       'True (Large C_v)', 'Estimated (Large C_v)');
title('QSSA Model');
grid on;

%%
true_params_basic = [delta, p_true, c_values(2);  % c small
                     delta, p_true, c_values(1)]; % c large

true_params_qssa = [delta, p_true / c_values(2), c_values(2);
                    delta, p_true / c_values(1), c_values(1)];

% 
figure;
subplot(1,2,1);
bar_params = categorical({'Delta', 'p', 'c'});
bar_values_basic = [true_params_basic(2,:); params_est_basic(1, 2:end); 
                    true_params_basic(1,:); params_est_basic(2, 2:end)];

barh(bar_params, bar_values_basic'); % 
ylabel('Parameter'); % 
xlabel('Estimated Values'); %
legend('True (Small C_v)', 'Estimated (Small C_v)', ...
       'True (Large C_v)', 'Estimated (Large C_v)', 'Location', 'best');
title('Basic Viral Model');
grid on;

%
subplot(1,2,2);
bar_params_qssa = categorical({'Delta', 'p''', 'c'});  % p -> p'
bar_values_qssa = [true_params_qssa(2, :); params_est_qssa(1, 2:end); 
                   true_params_qssa(1,:); params_est_qssa(2, 2:end)];

barh(bar_params_qssa, bar_values_qssa'); % 
ylabel('Parameter'); % 
xlabel('Estimated Values'); % 
legend('True (Small C_v)', 'Estimated (Small C_v)', ...
       'True (Large C_v)', 'Estimated (Large C_v)', 'Location', 'best');
title('QSSA Model');
grid on;

%% 
figure;
true_params_qssa(1,2) = 11000;
true_params_qssa(2,2) = 11000;
params_est_qssa(1,3) = params_est_qssa(1,3) * params_est_qssa(1,4);
params_est_qssa(2,3) = params_est_qssa(2,3) * params_est_qssa(2,4);

% Define parameters for labeling
bar_params = categorical({'Basic', 'QSSA'}); 

% Define colors
colors_true_basic = [0.5 0.5 1]; % Light Blue (Basic True)
colors_true_qssa  = [1 0.6 0.6]; % Light Red (QSSA True)
colors_est_basic  = [0 0 1];     % Dark Blue (Basic Estimated)
colors_est_qssa   = [1 0 0];     % Dark Red (QSSA Estimated)

% 
for i = 1:3
   
    subplot(3,2,2*i-1);
    
    % 
    b1 = barh(bar_params(1), true_params_basic(2, i), 'FaceColor', colors_true_basic); hold on;
    b2 = barh(bar_params(2), true_params_qssa(2, i), 'FaceColor', colors_true_qssa);

    % Estimated
    b3 = barh(bar_params(1), params_est_basic(1, i+1), 'FaceColor', colors_est_basic);
    b4 = barh(bar_params(2), params_est_qssa(1, i+1), 'FaceColor', colors_est_qssa);

    xlabel('Estimated Values');
    title(sprintf('%s (Small C_v)', char(bar_params_qssa(i))));
    legend([b1, b2, b3, b4], {'Basic True', 'QSSA True', 'Basic Estimated', 'QSSA Estimated'}, 'Location', 'best');
    grid on;

    % 
    subplot(3,2,2*i);
    
    % True
    b5 = barh(bar_params(1), true_params_basic(1, i), 'FaceColor', colors_true_basic); hold on;
    b6 = barh(bar_params(2), true_params_qssa(1, i), 'FaceColor', colors_true_qssa);

    % Estimated 
    b7 = barh(bar_params(1), params_est_basic(2, i+1), 'FaceColor', colors_est_basic);
    b8 = barh(bar_params(2), params_est_qssa(2, i+1), 'FaceColor', colors_est_qssa);

    xlabel('Estimated Values');
    title(sprintf('%s (Large C_v)', char(bar_params_qssa(i))));
    legend([b5, b6, b7, b8], {'Basic True', 'QSSA True', 'Basic Estimated', 'QSSA Estimated'}, 'Location', 'best');
    grid on;
end

% Adjust layout
sgtitle('Comparison of Basic and QSSA Models for Each Parameter');



%%  ODE Model
function dydt = viral_model(~, y, params)
    beta = params(1);
    delta = params(2);
    p = params(3);
    c = params(4);
    dydt = [-beta * y(1) * y(3); beta * y(1) * y(3) - delta * y(2); p * y(2) - c * y(3)];
end

function dydt = qssa_model(~, y, params)
    beta = params(1);
    delta = params(2);
    p_prime = params(3); % p' = p / c
    c = params(4);

    dydt = [-beta * p_prime * y(1) * y(2); 
            beta * p_prime * y(1) * y(2) - delta * y(2)];
end

function F = error_qssa(params, data, T0, V0, V)
    tspan = data(:,1);
    beta_prime = params(1);
    delta = params(2);
    c = params(3);
    p = params(4);

    % QSSA 
    eta = (-beta_prime * V0 + beta_prime * V0 * exp(-1) - delta) / c;
    I0_qssa = (beta_prime * T0 * V0 * exp(eta)) / c;

    Initials = [T0 * V, I0_qssa];

    % QSSA Model
    [~, Y] = ode45(@(t,y) qssa_model(t, y, params), tspan, Initials);
    
  
    V_fit = (p / c) * Y(:,2);  % V ≈ (p/c) * I

    F = abs(data(:,4) - V_fit);
end

function F = error_basic(params, data, T0, V0, V)
    tspan = data(:,1);
    beta = params(1);
    delta = params(2);
    p = params(3);
    c = params(4);
    
  
    Initials = [T0 * V, 0, V0 * V];

 
    [~, Y] = ode45(@(t,y) viral_model(t, y, params), tspan, Initials);
    
    V_fit = Y(:,3);

   
    F = abs(data(:,4) - V_fit);
end