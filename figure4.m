%% Viral dynamics comparison: Local Sensitivity Analysis
clear
clc
close all
set(groot, 'DefaultAxesFontSize', 16, 'DefaultAxesFontWeight', 'bold', ...
           'DefaultTextFontSize', 16, 'DefaultTextFontWeight', 'bold');
set(0, 'DefaultLineLineWidth', 3);
rng(11)

% Define parameters
beta = 3.15e-7;  % Infection rate
delta = 2.1;     % Clearance rate of infected cells
p = 11000;       % Production rate of virus
c = 10;          % Clearance rate of virus

% Initial conditions for Basic Viral Model
T0 = 100000;    
I0 = 0;        
V0 = 10000;    
y0_basic = [T0; I0; V0];

% Initial conditions for QSSA Model
T_r0 = 100000;  
I_r0 = beta * T0 * V0 * exp((beta * V0 * exp(-1) - beta * V0 - delta) / c) / c;
y0_qssa = [T_r0; I_r0];

% Initial conditions for QSSA_{il} Model
T_il0 = 100000;
V_il0 = 10000;
y0_qssa_il = [T_il0; V_il0];

% Time span
tspan = [0 10];
t_eval = linspace(tspan(1), tspan(2), 1000);  % Evaluation time points

%% Solve the ODEs for original parameters
[t_basic, y_basic] = ode45(@(t, y) basic_viral_model(t, y, beta, delta, p, c), tspan, y0_basic);
[t_qssa, y_qssa] = ode45(@(t, y) qssa_model(t, y, beta, delta, p, c), tspan, y0_qssa);
[t_qssa_il, y_qssa_il] = ode45(@(t, y) qssa_il_model(t, y, beta, delta, p, c), tspan, y0_qssa_il);

% Interpolate to the same time points
y_basic = interp1(t_basic, y_basic, t_eval);
y_qssa = interp1(t_qssa, y_qssa, t_eval);
y_qssa_il = interp1(t_qssa_il, y_qssa_il, t_eval);

% Parameters array
params = [beta, delta, p, c];
delta_param = 0.3;  % 20% variation
param_names = {'\beta', '\delta', 'p', 'c'};

%% Sensitivity analysis for Basic Viral Model
sensitivity_basic = struct('T', [], 'I', [], 'V', []);
for i = 1:4
    result = parameter_sensitivity(@basic_viral_model, y0_basic, t_eval, params, i, delta_param, 'basic', p, c);
    sensitivity_basic.T(:, i) = result.T;
    sensitivity_basic.I(:, i) = result.I;
    sensitivity_basic.V(:, i) = result.V;
end

%% Sensitivity analysis for QSSA Model
sensitivity_qssa = struct('T', [], 'I', [], 'V', []);
for i = 1:4
    result = parameter_sensitivity(@qssa_model, y0_qssa, t_eval, params, i, delta_param, 'qssa', p, c);
    sensitivity_qssa.T(:, i) = result.T;
    sensitivity_qssa.I(:, i) = result.I;
    sensitivity_qssa.V(:, i) = (p/c) * result.I; % 변환식 적용
end

%% Sensitivity analysis for QSSA_{il} Model
sensitivity_qssa_il = struct('T', [], 'I', [], 'V', []);
for i = 1:4
    result = parameter_sensitivity(@qssa_il_model, y0_qssa_il, t_eval, params, i, delta_param, 'qssa_il', p, c);
    sensitivity_qssa_il.T(:, i) = result.T;
    sensitivity_qssa_il.V(:, i) = result.V;  
    sensitivity_qssa_il.I(:, i) = (c/p) * result.V; % I = (c/p) * V 
end

%%  Local Sensitivity Comparison**
figure;

% T (Target Cells) Sensitivity
subplot(3,2,1);
plot_sensitivity_comparison(t_eval, y_basic(:,1), sensitivity_basic.T, param_names, '-k', 'Basic Viral Model - T');

subplot(3,2,2);
plot_sensitivity_comparison(t_eval, y_qssa(:,1), sensitivity_qssa.T, param_names, '-k', 'QSSA Model - T');

% subplot(3,3,3);
% plot_sensitivity_comparison(t_eval, y_qssa_il(:,1), sensitivity_qssa_il.T, param_names, '-k', 'QSSA_{il} Model - T');

% I (Infected Cells) Sensitivity
subplot(3,2,3);
plot_sensitivity_comparison(t_eval, y_basic(:,2), sensitivity_basic.I, param_names, '-k', 'Basic Viral Model - I');

subplot(3,2,4);
plot_sensitivity_comparison(t_eval, y_qssa(:,2), sensitivity_qssa.I, param_names, '-k', 'QSSA Model - I');

% subplot(3,3,6);
% plot_sensitivity_comparison(t_eval, y_qssa_il(:,2), sensitivity_qssa_il.I, param_names, '-k', 'QSSA_{il} Model - I');

% V (Virus) Sensitivity
subplot(3,2,5);
plot_sensitivity_comparison(t_eval, y_basic(:,3), sensitivity_basic.V, param_names, '-k', 'Basic Viral Model - V');

subplot(3,2,6);
plot_sensitivity_comparison(t_eval, y_qssa(:,2) * (p/c), sensitivity_qssa.V, param_names, '-k', 'QSSA Model - V');

% subplot(3,3,9);
% plot_sensitivity_comparison(t_eval, y_qssa_il(:,2), sensitivity_qssa_il.V, param_names, '-k', 'QSSA_{il} Model - V');
%% 
figure;

% T Mean Relative Change
subplot(1,3,1);
bar([mean(abs(sensitivity_basic.T - y_basic(:,1)), 1) ./ mean(y_basic(:,1)) * 100;
     mean(abs(sensitivity_qssa.T - y_qssa(:,1)), 1) ./ mean(y_qssa(:,1)) * 100]');
set(gca, 'XTickLabel', param_names);
xlabel('Parameter');
ylabel('Mean Relative Change in T (%)');
legend({'Basic Viral', 'QSSA'});
%title('Mean Relative Change in T');

% I Mean Relative Change
subplot(1,3,2);
bar([mean(abs(sensitivity_basic.I - y_basic(:,2)), 1) ./ mean(y_basic(:,2)) * 100;
     mean(abs(sensitivity_qssa.I - y_qssa(:,2)), 1) ./ mean(y_qssa(:,2)) * 100 ]');
set(gca, 'XTickLabel', param_names);
xlabel('Parameter');
ylabel('Mean Relative Change in I (%)');
legend({'Basic Viral', 'QSSA'});
%title('Mean Relative Change in I');

% % V Mean Relative Change
subplot(1,3,3);
bar([mean(abs(sensitivity_basic.V - y_basic(:,3)), 1) ./ mean(y_basic(:,3)) * 100;
     mean(abs(sensitivity_qssa.V - p/c * y_qssa(:,2)), 1) ./ mean(p/c * y_qssa(:,2)) * 100
     ]');
set(gca, 'XTickLabel', param_names);
xlabel('Parameter');
ylabel('Mean Relative Change in V (%)');
legend({'Basic Viral', 'QSSA'});
%title('Mean Relative Change in V');

%% Sensitivity Trend: More Detailed Virus Analysis
figure;
colors = {'r', 'g', 'b', 'm'};

for i = 1:4
    %Virus (V) Sensitivity Trend in Basic Viral Model
    subplot(2,4,i);
    plot(t_eval, abs(sensitivity_basic.V(:,i) - y_basic(:,3)) ./ y_basic(:,3) * 100, colors{i}, 'LineWidth', 2);
    xlabel('Time'); ylabel('Rel. Change (%)');
    title(['Basic Viral -' param_names{i}]);
    grid on;

    % Virus (V) Sensitivity Trend in QSSA Model
    subplot(2,4,i+4);
    plot(t_eval, abs(sensitivity_qssa.V(:,i) - (p/c) * y_qssa(:,2)) ./ (p/c * y_qssa(:,2)) * 100, colors{i}, 'LineWidth', 2);
    xlabel('Time'); ylabel('Rel. Change (%)');
    title(['QSSA - ' param_names{i}]);
    grid on;
    % 
    % % Virus (V) Sensitivity Trend in QSSA_{il} Model
    % subplot(3,4,i+8);
    % plot(t_eval, abs(sensitivity_qssa_il.V(:,i) - y_qssa_il(:,2)) ./ y_qssa_il(:,2) * 100, colors{i}, 'LineWidth', 2);
    % xlabel('Time'); ylabel('Rel. Change (%)');
    % title(['QSSA_{il} - V (' param_names{i} ')']);
    % grid on;
end

sgtitle('Parameter Sensitivity Trend (Virus V)');
%%
function dydt = basic_viral_model(t, y, beta, delta, p, c)
    dydt = [-beta * y(1) * y(3); 
             beta * y(1) * y(3) - delta * y(2); 
             p * y(2) - c * y(3)];
end

function dydt = qssa_model(t, y, beta, delta, p, c)
    dydt = [-beta * (p/c) * y(1) * y(2); 
             beta * (p/c) * y(1) * y(2) - delta * y(2)];
end

function dydt = qssa_il_model(t, y, beta, delta, p, c)
    dydt = [-beta * y(1) * y(2); 
             beta * y(1) * y(2) * (p/c) - delta * y(2)];
end

function sensitivity = parameter_sensitivity(model, y0, t_eval, params, param_index, delta_param, model_type, p, c)
    params(param_index) = params(param_index) * (1 + delta_param);
    [~, y] = ode45(@(t, y) model(t, y, params(1), params(2), params(3), params(4)), t_eval, y0);
    y = interp1(t_eval, y, t_eval);

    switch model_type
        case 'basic'
            sensitivity.T = y(:, 1);
            sensitivity.I = y(:, 2);
            sensitivity.V = y(:, 3);
        case 'qssa'
            sensitivity.T = y(:, 1);
            sensitivity.I = y(:, 2);
            sensitivity.V = p/c * y(:, 2); 
        case 'qssa_il'
            sensitivity.T = y(:, 1);
            sensitivity.V = y(:, 2);
            sensitivity.I = c/p * y(:, 2); 
    end
end
function plot_sensitivity_comparison(t_eval, y_original, sensitivity, param_names, original_style, title_text)
    plot(t_eval, y_original, original_style, 'LineWidth', 3, 'DisplayName', 'Original');
    hold on;
    colors = {'-r', '-g', '-b', '-m'};
    for i = 1:4
        plot(t_eval, sensitivity(:,i), colors{i}, 'LineWidth', 2, ...
             'DisplayName', [param_names{i} ' + 20%']);
    end
    xlabel('Time');
    ylabel(title_text);
    legend;
    title(title_text);
    grid on;
end