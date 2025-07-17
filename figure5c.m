%% **C_v vs. Relative Error Plots (Separated)**
clear
clc

set(groot, 'DefaultAxesFontSize', 16, 'DefaultAxesFontWeight', 'bold', ...
           'DefaultTextFontSize', 16, 'DefaultTextFontWeight', 'bold');
set(0, 'DefaultLineLineWidth', 3);
rng(11)

% Fixed parameters
delta = 0.1;  
c_values = linspace(10,100,20); % Varying clearance rates (c)
Cv_values = delta ./ c_values;  % Compute C_v values

% Storage for relative errors
rel_err_T = zeros(size(Cv_values));
rel_err_I = zeros(size(Cv_values));
rel_err_V = zeros(size(Cv_values));

b = 3.15e-7;
p = 11000;
Initials = [100000, 0, 10000];  % T(0), I(0), V(0)

% Compute relative errors for different c values
for i = 1:length(c_values)
    c = c_values(i);
    pp = p / c;
    params_full = [b, delta, p, c]; % Basic Viral Model parameters
    params_qssa = [b * pp, delta];  % QSSA Model parameters

    % Solve Basic Viral Model
    fun_full = @(t,y,params) [-params(1)*y(1)*y(3); 
                               params(1)*y(1)*y(3)-params(2)*y(2); 
                               params(3)*y(2)-params(4)*y(3)];
    tspan = linspace(0, 10, 2000);
    [t, viral_tot] = ode45(@(t,y) fun_full(t,y,params_full), tspan, Initials);

    % Solve QSSA Model
Init_qssa = [Initials(1), b * Initials(1) * Initials(3) * exp((b * Initials(3) * exp(-1) - b * Initials(3) - delta) / c) / c];
fun_qssa = @(t,y,params) [-params(1)*y(1)*y(2); params(1)*y(1)*y(2)-params(2)*y(2)];
    [t, vi_qssa] = ode45(@(t,y) fun_qssa(t,y,params_qssa), tspan, Init_qssa);
    v_qssa = pp * vi_qssa(:,2);

    % Compute Norm-Based Relative Errors
    rel_err_T(i) = norm(viral_tot(5:end,1) - vi_qssa(5:end,1)) / norm(viral_tot(5:end,1));
    rel_err_I(i) = norm(viral_tot(5:end,2) - vi_qssa(5:end,2)) / norm(viral_tot(5:end,2));
    rel_err_V(i) = norm(viral_tot(5:end,3) - v_qssa(5:end)) / norm(viral_tot(5:end,3));
end

% Plot: C_v vs. Relative Error for Target Cells (T)
figure;
scatter(Cv_values, rel_err_T, 100, 'b', 'filled', 'DisplayName', 'Relative Error (T)');
hold on;
xline(0.1, '--', 'C_v = 0.1', 'Color', 'k', 'LineWidth', 2, 'LabelVerticalAlignment', 'middle');
yline(0.1, '--', 'Rel. Error = 0.1', 'Color', [0.5 0 0.5], 'LineWidth', 2, 'LabelHorizontalAlignment', 'right');
xlabel('C_v = \delta / c');
ylabel('Relative Error (T)');
grid on;
legend show;
title('C_v vs. Relative Error for Target Cells (T)');

% Plot: C_v vs. Relative Error for Infected Cells (I)
figure;
scatter(Cv_values, rel_err_I, 100, 'r', 'filled', 'DisplayName', 'Relative Error (I)');
hold on;
xline(0.1, '--', 'C_v = 0.1', 'Color', 'k', 'LineWidth', 2, 'LabelVerticalAlignment', 'middle');
yline(0.1, '--', 'Rel. Error = 0.1', 'Color', [0.5 0 0.5], 'LineWidth', 2, 'LabelHorizontalAlignment', 'right');
xlabel('C_v = \delta / c');
ylabel('Relative Error (I)');
grid on;
legend show;
title('C_v vs. Relative Error for Infected Cells (I)');

% Plot: C_v vs. Relative Error for Virus (V)
figure;
scatter(Cv_values, rel_err_V, 100, 'g', 'filled', 'DisplayName', 'Relative Error (V)');
hold on;
xline(0.1, '--', 'C_v = 0.1', 'Color', 'k', 'LineWidth', 2, 'LabelVerticalAlignment', 'middle');
yline(0.1, '--', 'Rel. Error = 0.1', 'Color', [0.5 0 0.5], 'LineWidth', 2, 'LabelHorizontalAlignment', 'right');
xlabel('C_v = \delta / c');
ylabel('Relative Error (V)');
grid on;
legend show;
title('C_v vs. Relative Error for Virus (V)');