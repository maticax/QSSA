%% Viral dynamics comparison: Basic Viral Model vs. QSSA Model
clear
clc
close all
set(groot, 'DefaultAxesFontSize', 16, 'DefaultAxesFontWeight', 'bold', ...
           'DefaultTextFontSize', 16, 'DefaultTextFontWeight', 'bold');
set(0, 'DefaultLineLineWidth', 3);
rng(11)

% Parameter Set
delta = 2.1;
c = 100;  % Viral clearance rate
params_full = [3.15e-7, delta, 11000, c]; % [b, delta, p, c]

b = params_full(1);
p = params_full(3);
pp = p/c;
params_qssa = [b*pp, delta]; 
beta = b;

% Time span
tspan = linspace(0, 10, 2000);
dt = tspan(2) - tspan(1);

% Initial values for the Basic Viral Model
Initials = [100000, 0, 10000];  % T(0), I(0), V(0)

% Compute a range of initial values for I(0) in QSSA model
eta1 = (beta * Initials(3) * exp(-1) - c) / c;
eta2 = (-beta * Initials(3) + beta * Initials(3) * exp(-1) - delta) / c;
eta = linspace(eta1, eta2, 10);

% Define QSSA initial values
Init_qssa_vec = beta * Initials(1) * Initials(3) * exp(eta) / c; 

% Solve the Basic Viral Model
fun_full = @(t,y,params) [-params(1)*y(1)*y(3); 
                           params(1)*y(1)*y(3)-params(2)*y(2); 
                           params(3)*y(2)-params(4)*y(3)];
[t, viral_tot] = ode45(@(t,y) fun_full(t,y,params_full), tspan, Initials);

% Figure setup
figure; % Single figure

% Color setup: Basic Viral Model (blue), QSSA (various colors)
basic_color = [0, 0, 1]; % Blue for Basic Viral Model
qssa_colors = lines(3);  % Using MATLAB's 'lines' colormap for QSSA variants
qssa_styles = {'--', '-.', ':'};  % Different line styles for QSSA models

% Define handle variables for legend
h_basic = []; % Store handle for Basic Viral Model
h_qssa = [];  % Store handles for QSSA models

% Define selected initial values for QSSA
selected_indices = [1, 5, 10];  % Selecting three different QSSA initial values

% Storage for norm-based errors
norm_rel_err_target = zeros(1, length(selected_indices));
norm_rel_err_infect = zeros(1, length(selected_indices));
norm_rel_err_virus = zeros(1, length(selected_indices));

for i_idx = 1:length(selected_indices)  
    i = selected_indices(i_idx);
    Init_qssa = [Initials(1), Init_qssa_vec(i)];

    % Solve QSSA Model
    fun_qssa = @(t,y,params) [-params(1)*y(1)*y(2); params(1)*y(1)*y(2)-params(2)*y(2)];
    [t, vi_qssa] = ode45(@(t,y) fun_qssa(t,y,params_qssa), tspan, Init_qssa);

    % Reconstruct Virus Concentration
    v = pp * vi_qssa(:,2);

    % Compute Norm-Based Relative Errors
    norm_rel_err_target(i_idx) = norm(viral_tot(:,1) - vi_qssa(:,1)) / norm(viral_tot(:,1));
    norm_rel_err_infect(i_idx) = norm(viral_tot(:,2) - vi_qssa(:,2)) / norm(viral_tot(:,2));
    norm_rel_err_virus(i_idx) = norm(viral_tot(:,3) - v) / norm(viral_tot(:,3));

    % **Plot Target Cells (T)**
    subplot(3,1,1) 
    if i_idx == 1
        h_basic = plot(t, viral_tot(:,1), '-', 'Color', basic_color, 'LineWidth', 3, 'DisplayName', 'Basic Viral Model'); 
        hold on;
    end
    h_qssa(i_idx) = semilogy(t, vi_qssa(:,1), qssa_styles{i_idx}, 'Color', qssa_colors(i_idx, :), 'LineWidth', 2, ...
                          'DisplayName', sprintf('QSSA I_{0,%d}, (Rel.Err: %.3f)', i, norm_rel_err_target(i_idx)));
    ylabel('T');
    legend show;
    grid on;

    % **Plot Infected Cells (I)**
    subplot(3,1,2)
    if i_idx == 1
        plot(t, viral_tot(:,2), '-', 'Color', basic_color, 'LineWidth', 3, 'DisplayName', 'Basic Viral Model'); 
        hold on;
    end
    semilogy(t, vi_qssa(:,2), qssa_styles{i_idx}, 'Color', qssa_colors(i_idx, :), 'LineWidth', 2, ...
         'DisplayName', sprintf('Rel.Err: %.3f',  norm_rel_err_infect(i_idx)));
    ylabel('I');
    legend show;
    grid on;

    % **Plot Virus Concentration (V)**
    subplot(3,1,3)
    if i_idx == 1
        plot(t, viral_tot(:,3), '-', 'Color', basic_color, 'LineWidth', 3, 'DisplayName', 'Basic Viral Model'); 
        hold on;
    end
    semilogy(t, v, qssa_styles{i_idx}, 'Color', qssa_colors(i_idx, :), 'LineWidth', 2, ...
         'DisplayName', sprintf('Rel.Err: %.3f',  norm_rel_err_virus(i_idx)));
    ylabel('V');
    xlabel('Time');
    legend show;
    grid on;
end

% Final adjustments for legends
subplot(3,1,1) % Target Cells
legend([h_basic, h_qssa(1), h_qssa(2), h_qssa(3)], ...
       {'Basic Viral Model', ...
        sprintf('QSSA (I_{0,1}) Rel.Err: %.3f', norm_rel_err_target(1)), ...
        sprintf('QSSA (I_{0,5}) Rel.Err: %.3f', norm_rel_err_target(2)), ...
        sprintf('QSSA (I_{0,10}) Rel.Err: %.3f', norm_rel_err_target(3))});
grid on;

subplot(3,1,2) % Infected Cells
legend([h_basic, h_qssa(1), h_qssa(2), h_qssa(3)], ...
       {'Basic Viral Model', ...
        sprintf('QSSA (I_{0,1}) Rel.Err: %.3f', norm_rel_err_infect(1)), ...
        sprintf('QSSA (I_{0,5}) Rel.Err: %.3f', norm_rel_err_infect(2)), ...
        sprintf('QSSA (I_{0,10}) Rel.Err: %.3f', norm_rel_err_infect(3))});
grid on;

subplot(3,1,3) % Virus
legend([h_basic, h_qssa(1), h_qssa(2), h_qssa(3)], ...
       {'Basic Viral Model', ...
        sprintf('QSSA (I_{0,1}) Rel.Err: %.3f', norm_rel_err_virus(1)), ...
        sprintf('QSSA (I_{0,5}) Rel.Err: %.3f', norm_rel_err_virus(2)), ...
        sprintf('QSSA (I_{0,10}) Rel.Err: %.3f', norm_rel_err_virus(3))});
grid on;
sgtitle('Viral model vs. QSSA model')



%% **Second Figure: Basic Viral Model vs. QSSA_il**
figure;

% Solve QSSA_il Model (Single Run)
Init_qssa_il = [Initials(1), Initials(3)];
fun_qssa_il = @(t,y,params) [-params(1)*y(1)*y(2); params(1) * (params(3) / params(4)) * y(1) * y(2) - params(2) * y(2)];
[t, vi_qssa_il] = ode45(@(t,y) fun_qssa_il(t,y,params_full), tspan, Init_qssa_il);
I_il = c / p * vi_qssa_il(:,2);

% Compute Norm-Based Relative Errors
norm_rel_err_target = norm(viral_tot(:,1) - vi_qssa_il(:,1)) / norm(viral_tot(:,1));
norm_rel_err_infect = norm(viral_tot(:,2) - (vi_qssa_il(:,2) / pp)) / norm(viral_tot(:,2));
norm_rel_err_virus = norm(viral_tot(:,3) - I_il) / norm(viral_tot(:,3));

% **Plot Target Cells (T)**
subplot(3,1,1)
semilogy(t, viral_tot(:,1), '-', 'Color', 'b', 'LineWidth', 3, 'DisplayName', 'Basic Viral Model');
hold on;
plot(t, vi_qssa_il(:,1), '--', 'Color', 'r', 'LineWidth', 2, ...
     'DisplayName', sprintf('QSSA_{il} Rel.Err: %.3f', norm_rel_err_target));
ylabel('T');
legend show;
grid on;

% **Plot Infected Cells (I)**
subplot(3,1,2)
semilogy(t, viral_tot(:,2), '-', 'Color', 'b', 'LineWidth', 3, 'DisplayName', 'Basic Viral Model');
hold on;
plot(t, I_il, '--', 'Color', 'r', 'LineWidth', 2, ...
     'DisplayName', sprintf('Rel.Err: %.3f', norm_rel_err_virus));
ylabel('I');
legend show;
grid on;

% **Plot Virus Concentration (V)**
subplot(3,1,3)
semilogy(t, viral_tot(:,3), '-', 'Color', 'b', 'LineWidth', 3, 'DisplayName', 'Basic Viral Model');
hold on;
plot(t, vi_qssa_il(:,2) / pp, '--', 'Color', 'r', 'LineWidth', 2, ...
     'DisplayName', sprintf('Rel.Err: %.3f', norm_rel_err_infect));
ylabel('V');
xlabel('Time');
legend show;
grid on;

sgtitle('Viral Model vs. QSSA_{il}')