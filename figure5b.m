%% **Validity Condition C_v = δ/c에 따른 다이나믹스 및 상대 오차 분석**
clear
clc
close all
set(groot, 'DefaultAxesFontSize', 16, 'DefaultAxesFontWeight', 'bold', ...
           'DefaultTextFontSize', 16, 'DefaultTextFontWeight', 'bold');
set(0, 'DefaultLineLineWidth', 3);
rng(11)

% Define fixed parameters
beta = 3.15e-7;  % Infection rate
p = 11000;       % Virus production rate
N = 100000;      % Total cell count
V0 = 10000;      % Initial virus particles

% Selected delta and c values
%delta_vals = [1, 4, 10];  % Clearance rate of infected cells
%c_vals = [10, 40, 100];   % Clearance rate of virus

delta_vals = [0.1, 0.4, 1];  % Clearance rate of infected cells
c_vals = [10, 40, 100];   % Clearance rate of virus
% Compute C_v = delta/c values
C_v_vals = delta_vals ./ c_vals;

% Time span
tspan = linspace(0, 10, 1000);

% Define colors and line styles
basic_colors = lines(length(delta_vals)); 
qssa_colors = basic_colors * 0.7;         
line_styles = {'-', '--', '-.'};          
markers = {'o', 's', 'd'};                

%% Figure for Target Cells (T, I, V)
figure;

for i = 1:length(delta_vals)
    delta = delta_vals(i);
    c = c_vals(i);
    pp = p / c;
    
    % Initial conditions
    y0_basic = [N; 0; V0];
    I_r0 = beta * N * V0 * exp((beta * V0 * exp(-1) - beta * V0 - delta) / c) / c;
    y0_qssa = [N; I_r0];

    % Parameters
    params_full = [beta, delta, p, c];
    params_qssa = [beta * pp, delta];

    % Solve Basic Viral Model
    [t, y_basic] = ode45(@(t,y) basic_viral_model(t,y,params_full), tspan, y0_basic);

    % Solve QSSA Model
    [t, y_qssa] = ode45(@(t,y) qssa_model(t,y,params_qssa), tspan, y0_qssa);
    v_qssa = pp * y_qssa(:,2);

    % Compute Norm-Based Relative Errors
    norm_rel_err_T = norm(y_basic(5:end,1) - y_qssa(5:end,1)) / norm(y_basic(5:end,1));
    norm_rel_err_I = norm(y_basic(5:end,2) - y_qssa(5:end,2)) / norm(y_basic(5:end,2));
    norm_rel_err_V = norm(y_basic(5:end,3) - v_qssa(5:end)) / norm(y_basic(5:end,3));

    % Plot Target Cells (T)
    subplot(3,1,1);
    semilogy(t, y_basic(:,1), line_styles{i}, 'Color', basic_colors(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Basic: \\delta=%.1f, c=%.1f, C_v=%.2f', delta, c, C_v_vals(i)));
    hold on;
    semilogy(t, y_qssa(:,1), line_styles{i}, 'Color', qssa_colors(i,:), 'LineWidth', 2, ...
        'Marker', markers{i}, 'MarkerIndices', 1:100:length(t), 'MarkerSize', 6, ...
        'DisplayName', sprintf('QSSA: \\delta=%.1f, c=%.1f, C_v=%.2f (RE: %.3f)', delta, c, C_v_vals(i), norm_rel_err_T));
ylim([1e-5 inf])
    % Plot Infected Cells (I)
    subplot(3,1,2);
    semilogy(t, y_basic(:,2), line_styles{i}, 'Color', basic_colors(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Basic: \\delta=%.1f, c=%.1f, C_v=%.2f', delta, c, C_v_vals(i)));
    hold on;
    semilogy(t, y_qssa(:,2), line_styles{i}, 'Color', qssa_colors(i,:), 'LineWidth', 2, ...
        'Marker', markers{i}, 'MarkerIndices', 1:100:length(t), 'MarkerSize', 6, ...
        'DisplayName', sprintf('QSSA: \\delta=%.1f, c=%.1f, C_v=%.2f (RE: %.3f)', delta, c, C_v_vals(i), norm_rel_err_I));

    % Plot Virus (V)
    subplot(3,1,3);
    semilogy(t, y_basic(:,3), line_styles{i}, 'Color', basic_colors(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Basic: \\delta=%.1f, c=%.1f, C_v=%.2f', delta, c, C_v_vals(i)));
    hold on;
    semilogy(t, v_qssa, line_styles{i}, 'Color', qssa_colors(i,:), 'LineWidth', 2, ...
        'Marker', markers{i}, 'MarkerIndices', 1:100:length(t), 'MarkerSize', 6, ...
        'DisplayName', sprintf('QSSA: \\delta=%.1f, c=%.1f, C_v=%.2f (RE: %.3f)', delta, c, C_v_vals(i), norm_rel_err_V));
end

% Final adjustments for each subplot
subplot(3,1,1);
ylabel('T');
legend show;
grid on;
title('Effect of C_v on T (Target Cells)');

subplot(3,1,2);
ylabel('I');
legend show;
grid on;
title('Effect of C_v on I (Infected Cells)');

subplot(3,1,3);
ylabel('V');
xlabel('Time');
legend show;
grid on;
title('Effect of C_v on V (Virus)');

sgtitle('Basic Viral Model vs. QSSA Model with Different Validity Conditions');

%% 
function dydt = basic_viral_model(t, y, params)
    beta = params(1);
    delta = params(2);
    p = params(3);
    c = params(4);

    dydt = [-beta * y(1) * y(3); 
             beta * y(1) * y(3) - delta * y(2); 
             p * y(2) - c * y(3)];
end

function dydt = qssa_model(t, y, params)
    beta_pp = params(1);
    delta = params(2);

    dydt = [-beta_pp * y(1) * y(2); 
             beta_pp * y(1) * y(2) - delta * y(2)];
end