%% Viral dynamics comparison: Basic Viral Model vs. QSSA & QSSA_il
clear
clc
close all
set(groot, 'DefaultAxesFontSize', 16, 'DefaultAxesFontWeight', 'bold', ...
           'DefaultTextFontSize', 16, 'DefaultTextFontWeight', 'bold');
set(0, 'DefaultLineLineWidth', 3);
rng(11)

% Parameter Set
delta = 2.1;
c_values = [5,10,20,40,60,80,100,120,150]; % Varying clearance rates
b = 3.15e-7;
p = 11000;
Initials = [100000, 0, 10000];  % T(0), I(0), V(0)

% Compute eta and fix I(0) for QSSA
eta1 = (b * Initials(3) * exp(-1) - c_values(end)) / c_values(end);
eta2 = (-b * Initials(3) + b * Initials(3) * exp(-1) - delta) / c_values(end);
eta = linspace(eta1, eta2, 10);
I_qssa_fixed = b * Initials(1) * Initials(3) * exp(eta(end)) / c_values(end); % Fix I(0) for QSSA

qssa_styles = {'--', '-.', ':'};  % Different line styles for QSSA models
qssa_colors = lines(length(c_values)); % Different colors for each c

for c_idx = 1:length(c_values)
    c = c_values(c_idx);
    pp = p / c;
    params_full = [b, delta, p, c]; % Parameters for Basic Viral Model
    params_qssa = [b * pp, delta];  % Parameters for QSSA Model
    params_qssa_il = params_full;   % Parameters for QSSA_il Model (Same as Basic Viral)

    % Solve the Basic Viral Model
    fun_full = @(t,y,params) [-params(1)*y(1)*y(3); 
                               params(1)*y(1)*y(3)-params(2)*y(2); 
                               params(3)*y(2)-params(4)*y(3)];
    tspan = linspace(0, 10, 2000);
    [t, viral_tot] = ode45(@(t,y) fun_full(t,y,params_full), tspan, Initials);

    % Solve QSSA Model (Fixed I(0))
    Init_qssa = [Initials(1), I_qssa_fixed];
    fun_qssa = @(t,y,params) [-params(1)*y(1)*y(2); params(1)*y(1)*y(2)-params(2)*y(2)];
    [t, vi_qssa] = ode45(@(t,y) fun_qssa(t,y,params_qssa), tspan, Init_qssa);
    v_qssa = pp * vi_qssa(:,2);

    % Solve QSSA_il Model
    Init_qssa_il = [Initials(1), Initials(3)];
    fun_qssa_il = @(t,y,params) [-params(1)*y(1)*y(2); params(1) * (params(3) / params(4)) * y(1) * y(2) - params(2) * y(2)];
    [t, vi_qssa_il] = ode45(@(t,y) fun_qssa_il(t,y,params_qssa_il), tspan, Init_qssa_il);
    I_il = c / p * vi_qssa_il(:,2);

    % Compute Norm-Based Relative Errors
    norm_rel_err_target_qssa = norm(viral_tot(:,1) - vi_qssa(:,1)) / norm(viral_tot(:,1));
    norm_rel_err_infect_qssa = norm(viral_tot(:,2) - vi_qssa(:,2)) / norm(viral_tot(:,2));
    norm_rel_err_virus_qssa = norm(viral_tot(:,3) - v_qssa) / norm(viral_tot(:,3));

    norm_rel_err_target_il = norm(viral_tot(:,1) - vi_qssa_il(:,1)) / norm(viral_tot(:,1));
    norm_rel_err_infect_il = norm(viral_tot(:,2) - (vi_qssa_il(:,2) / pp)) / norm(viral_tot(:,2));
    norm_rel_err_virus_il = norm(viral_tot(:,3) - I_il) / norm(viral_tot(:,3));

    % **Plot Figure for Current c**
    figure;
    sgtitle(sprintf('Basic Viral Model vs. QSSA vs. QSSA_{il} (c = %d)', c))

    % **Plot Target Cells (T)**
    subplot(3,1,1)
    semilogy(t, viral_tot(:,1), '-', 'Color', 'b', 'LineWidth', 3, 'DisplayName', 'Basic Viral Model');
    hold on;
    semilogy(t, vi_qssa(:,1), '--', 'Color', '[1.0, 0.5, 0.0]', 'LineWidth', 2, ...
         'DisplayName', sprintf('QSSA (Rel.Err: %.3f)', norm_rel_err_target_qssa));
    semilogy(t, vi_qssa_il(:,1), '-.', 'Color', 'r', 'LineWidth', 2, ...
         'DisplayName', sprintf('QSSA_{il} (Rel.Err: %.3f)', norm_rel_err_target_il));
    ylabel('T');
    legend show;
    grid on;

    % **Plot Infected Cells (I)**
    subplot(3,1,2)
    semilogy(t, viral_tot(:,2), '-', 'Color', 'b', 'LineWidth', 3, 'DisplayName', 'Basic Viral Model');
    hold on;
    semilogy(t, vi_qssa(:,2), '--', 'Color', '[1.0, 0.5, 0.0]', 'LineWidth', 2, ...
         'DisplayName', sprintf('Rel.Err: %.3f', norm_rel_err_infect_qssa));
    semilogy(t, I_il, '-.', 'Color', 'r', 'LineWidth', 2, ...
         'DisplayName', sprintf('Rel.Err: %.3f', norm_rel_err_infect_il));
    ylabel('I');
    legend show;
    grid on;

    % **Plot Virus Concentration (V)**
    subplot(3,1,3)
    semilogy(t, viral_tot(:,3), '-', 'Color', 'b', 'LineWidth', 3, 'DisplayName', 'Basic Viral Model');
    hold on;
    semilogy(t, v_qssa, '--', 'Color', [1.0, 0.5, 0.0], 'LineWidth', 2, ...
         'DisplayName', sprintf('Rel.Err: %.3f', norm_rel_err_virus_qssa));
    semilogy(t, vi_qssa_il(:,2) / pp, '-.', 'Color', 'r', 'LineWidth', 2, ...
         'DisplayName', sprintf('Rel.Err: %.3f', norm_rel_err_virus_il));
    ylabel('V');
    xlabel('Time');
    legend show;
    grid on;
end