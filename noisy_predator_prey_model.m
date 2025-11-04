% Code corresponding to the noisy predator-prey model in Section SI.4.2 - "A 
% Noisy Predator-Prey Model" of the paper "Assimilative Causal Inference".
%
% Authors: Marios Andreou, Nan Chen, Erik Bollt.
%
% Code Information: Application of Assimilative Causal Inference (ACI) to a 
% noisy version of the predator-prey (Lotka-Volterra) model with intermittent 
% regime switches:
%
%       (Predator Population) dx/dt = βxy - αx + σ_x·dot(W)_x,
%       (Prey Population) dy/dt = γy - δxy + σ_y·dot(W)_y.
%
% This is a stochastic conceptual model explaining the cyclical fluctuations 
% between the populations of species like lynx and hares. x is the population 
% density of a predator species, while y that of its prey. Their time 
% derivatives describe their instantaneous population growth rates, α is the 
% predator's natural death rate, β quantifies how prey availability enhances 
% predator growth, γ is the maximum intrinsic growth rate of prey, and δ 
% captures the negative impact of predators on the prey population. Here, since 
% β,δ>0, larger prey population y (y > α/β) enhances x via antidamping 
% mechanisms, while a larger predator population x (x < γ/δ) naturally
% intensifies the damping in y. The resulting coupled variations in x and y
% produce intermittent phase alternations and causal regime switches over time 
% t∈[0,T], from y(t) → x to x(t) → y and vice versa.  ACI is employed to study
% both causal relationships, therefore both x and y are assumed to be 
% observable variables. Since this model is a conditional Gaussian nonlinear 
% system (CGNS) for both x|y and y|x, then both of these posterior distributions 
% (filter and smoother) are Gaussian. This code uses the same model parameter 
% values as those cited in the paper.
%
% Written and tested in MATLAB R2024b.
%
% MATLAB Toolbox and M-file Requirements:
% 
% Code used to obtain the required m-file scripts and MATLAB toolboxes:
% [fList, pList] = matlab.codetools.requiredFilesAndProducts('noisy_predator_prey_model.m');
%
% M-file Scripts:
%   ➤ noisy_predator_prey_model.m
%   ➤ progress_bar.m
%   ➤ simps.m (https://www.mathworks.com/matlabcentral/fileexchange/25754-simpson-s-rule-for-numerical-integration)
%
% Data:
%   ➤ N/A
%
% Toolboxes:
%   ➤ N/A
%
% GitHub Repository: https://github.com/marandmath/ACI_code
% MIT License Information: https://github.com/marandmath/ACI_code/blob/main/LICENSE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MODEL SETUP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fixing the random number seed for reproducibility across each simulation.
rng(42)

% Total number of time steps within the given time interval.
N = 12000;
% Numerical integration time step.
dt = 0.005;
% Total simulation time.
T = N*dt;

% Predator population.
x = zeros(1, N+1);
% Prey population.
y = zeros(1, N+1);

% Predator's natural death rate.
alpha = 0.4;
% Effect of prey availability to predator growth.
beta = 0.1;
% Effect of predator presence to prey decline.
gamma = 1.1;
% Prey's natural growth rate.
delta = 0.4;

% Additive noise feedback in x.
sigma_x = 0.3;
% Additive noise feedback in y.
sigma_y = 0.3;

% Initial conditions.
x(1) = 4;
y(1) = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN THIS CODE SECTION TO STUDY THE CAUSAL RELATIONSHIP: x(t) → y
% y (Prey): Observed || x (Predator): Unobserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unobservable coefficient matrix: Feedback of x in x.
L_x = zeros(1, N+1);
% Forcing in the unobservable process.
f_x = 0;
% Noise feedback matrices in the unobservable process.
Sx_1 = sigma_x;
Sx_2 = 0;

% Observable coefficient matrix: Feedback of x in y.
L_y = zeros(1, N+1);
% Forcing in the observable process.
f_y = zeros(1, N+1);
% Noise feedback matrices in the observable process.
Sy_1 = 0;
Sy_2 = sigma_y;

% Initiating the time-dependent model components.
L_y(1) = - delta * y(1);
f_y(1) = gamma * y(1);
L_x(1) = beta * y(1) - alpha;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN THIS CODE SECTION TO STUDY THE CAUSAL RELATIONSHIP: y(t) → x
% x (Predator): Observed | y (Prey): Unobserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Observable coefficient matrix: Feedback of y in x.
L_x = zeros(1, N+1);
% Forcing in the observable process.
f_x = zeros(1, N+1);
% Noise feedback matrices in the observable process.
Sx_1 = sigma_x;
Sx_2 = 0;

% Unobservable coefficient matrix: Feedback of y in y.
L_y = zeros(1, N+1);
% Forcing in the unobservable process.
f_y = 0;
% Noise feedback matrices in the unobservable process.
Sy_1 = 0;
Sy_2 = sigma_y;

% Initiating the time-dependent model components.
L_x(1) = beta * x(1);
f_x(1) = - alpha * x(1);
L_y(1) = gamma - delta * x(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  GENERATING THE TRUE SIGNALS  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN THIS CODE SECTION TO STUDY THE CAUSAL RELATIONSHIP: x(t) → y
% y (Prey): Observed || x (Predator): Unobserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 2:N+1 

    % Wiener increments.
    dW_x = randn;
    dW_y = randn;

    % Updating the state variables based on the dynamical system equations.
    x(j) = x(j-1) + (L_x(j-1) * x(j-1) + f_x) * dt ...
                  + Sx_1 * sqrt(dt) * dW_x + Sx_2 * sqrt(dt) * dW_y;
    y(j) = y(j-1) + (L_y(j-1) * x(j-1) + f_y(j-1)) * dt ... 
                  + Sy_1 * sqrt(dt) * dW_x + Sy_2 * sqrt(dt) * dW_y;

    % Updating the time-dependent model components.
    L_y(j) = - delta * y(j);
    f_y(j) = gamma * y(j);
    L_x(j) = beta * y(j) - alpha;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN THIS CODE SECTION TO STUDY THE CAUSAL RELATIONSHIP: y(t) → x
% x (Predator): Observed | y (Prey): Unobserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 2:N+1 

    % Wiener increments.
    dW_x = randn;
    dW_y = randn;

    % Updating the state variables based on the dynamical system equations.
    x(j) = x(j-1) + (L_x(j-1) * y(j-1) + f_x(j-1)) * dt ...
                  + Sx_1 * sqrt(dt) * dW_x + Sx_2 * sqrt(dt) * dW_y;
    y(j) = y(j-1) + (L_y(j-1) * y(j-1) + f_y) * dt ... 
                  + Sy_1 * sqrt(dt) * dW_x + Sy_2 * sqrt(dt) * dW_y;

    % Updating the time-dependent model components.
    L_x(j) = beta * x(j);
    f_x(j) = - alpha * x(j);
    L_y(j) = gamma - delta * x(j);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FILTERING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN THIS CODE SECTION TO STUDY THE CAUSAL RELATIONSHIP: x(t) → y
% y (Prey): Observed || x (Predator): Unobserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sum of Grammians of the observable process' noise feedbacks.
S_yoS_y = Sy_1^2 + Sy_2^2;
% Inverse of the sum of the observational noise coefficients' Grammians.
S_yoS_y_inv = 1/S_yoS_y;
% Sum of Grammians of the unobservable process' noise feedbacks.
S_xoS_x =  Sx_1^2 + Sx_2^2;
% Noise cross-interactions between the observable and unobservable processes.
S_xoS_y = Sx_1*Sy_1 + Sx_2*Sy_2;
S_yoS_x = S_xoS_y.';

% Posterior filter mean of the latent variable x.
filter_mean = zeros(1, N+1);
% Initial value of the posterior filter mean.
filter_mean(1) = x(1);
% Posterior filter covariance matrix of the latent variable x.
filter_cov = zeros(1, N+1);
% Initial value of the posterior filter covariance. Choosing a positive definite 
% matrix as to preserve the positive-definiteness of the posterior covariance 
% matrices over time.
filter_cov(1) = 0.1;

mu0 = filter_mean(1);
R0 = filter_cov(1);

for j = 2:N+1
    
    dy = y(j) - y(j-1);
    aux = S_xoS_y + filter_cov(j-1) * L_y(j-1);

    % Update the posterior filter mean and posterior filter covariance using the
    % optimal nonlinear filter state estimation equations for CGNSs; See Section 
    % 2.1.2 of the Supplementary Information.
    mu = mu0 + (L_x(j-1) * mu0 + f_x) * dt ...
             + aux * S_yoS_y_inv * (dy - (L_y(j-1) * mu0 + f_y(j-1)) * dt);
    R = R0 + (L_x(j-1) * R0 + R0 * L_x(j-1) + S_xoS_x - aux * S_yoS_y_inv * aux) * dt;
    
    filter_mean(j) = mu;
    filter_cov(j) = R;
    mu0 = mu;
    R0 = R;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN THIS CODE SECTION TO STUDY THE CAUSAL RELATIONSHIP: y(t) → x
% x (Predator): Observed | y (Prey): Unobserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sum of Grammians of the observable process' noise feedbacks.
S_xoS_x = Sx_1^2 + Sx_2^2;
% Inverse of the sum of the observational noise coefficients' Grammians.
S_xoS_x_inv = 1/S_xoS_x;
% Sum of Grammians of the unobservable process' noise feedbacks.
S_yoS_y =  Sy_1^2 + Sy_2^2;
% Noise cross-interactions between the observable and unobservable processes.
S_yoS_x = Sy_1*Sx_1 + Sy_2*Sx_2;
S_xoS_y = S_yoS_x.';

% Posterior filter mean of the latent variable y.
filter_mean = zeros(1, N+1);
% Initial value of the posterior filter mean.
filter_mean(1) = y(1);
% Posterior filter covariance matrix of the latent variable y.
filter_cov = zeros(1, N+1);
% Initial value of the posterior filter covariance. Choosing a positive definite 
% matrix as to preserve the positive-definiteness of the posterior covariance 
% matrices over time.
filter_cov(1) = 0.1;

mu0 = filter_mean(1);
R0 = filter_cov(1);

for j = 2:N+1
    
    dx = x(j) - x(j-1);
    aux = S_yoS_x + filter_cov(j-1) * L_x(j-1);

    % Update the posterior filter mean and posterior filter covariance using the
    % optimal nonlinear filter state estimation equations for CGNSs; See Section 
    % 2.1.2 of the Supplementary Information.
    mu = mu0 + (L_y(j-1) * mu0 + f_y) * dt ...
             + aux * S_xoS_x_inv * (dx - (L_x(j-1) * mu0 + f_x(j-1)) * dt);
    R = R0 + (L_y(j-1) * R0 + R0 * L_y(j-1) + S_yoS_y - aux * S_xoS_x_inv * aux) * dt;
    
    filter_mean(j) = mu;
    filter_cov(j) = R;
    mu0 = mu;
    R0 = R;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SMOOTHING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN THIS CODE SECTION TO STUDY THE CAUSAL RELATIONSHIP: x(t) → y
% y (Prey): Observed || x (Predator): Unobserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Posterior smoother mean of the latent variable x.
smoother_mean = zeros(1, N+1); 
% Posterior smoother covariance matrix of the latent variable x.
smoother_cov = zeros(1, N+1);

% Smoother runs backwards: "Initial" values of the smoother statistics (i.e., at
% the last time instant T) are the corresponding posterior filter statistics.
smoother_mean(N+1) = filter_mean(N+1);
smoother_cov(N+1) = filter_cov(N+1);

% Auxiliary matrices used for the calculation of the online smoother for this 
% CGNS. The online smoother is required for the calculation of the subjective 
% and objective causal influence range (CIR) lengths of x(t) → y at each time 
% t∈[0,T]. Notation used is consistent with that of the original CGNS online 
% smoother work and the accompanying martingale-free introduction to CGNSs 
% paper: 
%   10.48550/arXiv.2411.05870 and 10.48550/arXiv.2410.24056
E_j_matrices = zeros(1, N+1);
F_j_matrices = zeros(1, N+1);
G_y_j = L_y(N+1) + S_yoS_x / filter_cov(N+1);
G_x_j = L_x(N+1) + S_xoS_x / filter_cov(N+1);
C_jj = 1 - G_x_j * dt;
H_j = filter_cov(N+1) \ (L_x(N+1) * filter_cov(N+1) + filter_cov(N+1) * L_x(N+1) + S_xoS_x);
K_j = S_yoS_y_inv * G_y_j;
E_j_matrices(N+1) = C_jj + S_xoS_y * K_j * dt;
F_j_matrices(N+1) = - filter_cov(N+1) * ( ...
                        K_j + (G_y_j * K_j * filter_cov(N+1) * K_j - filter_cov(N+1) \ H_j * filter_cov(N+1) * K_j + L_x(N+1) * K_j) * dt ...
                        - L_y(N+1) * (S_yoS_y_inv + K_j * filter_cov(N+1) * K_j * dt) ...
                      );

muT = smoother_mean(N+1);
RT = smoother_cov(N+1);

for j = N:-1:1

    % Calculation of the online smoother auxiliary matrices.
    G_y_j = L_y(j) + S_yoS_x / filter_cov(j);
    G_x_j = L_x(j) + S_xoS_x / filter_cov(j);
    C_jj = 1 - G_x_j * dt;
    H_j = filter_cov(j) \ (L_x(j) * filter_cov(j) + filter_cov(j) * L_x(j) + S_xoS_x);
    K_j = S_yoS_y_inv * G_y_j;
    E_j_matrices(j) = C_jj + S_xoS_y * K_j * dt;
    F_j_matrices(j) = - filter_cov(j) * ( ...
                            K_j + (G_y_j * K_j * filter_cov(j) * K_j - filter_cov(j) \ H_j * filter_cov(j) * K_j + L_x(j) * K_j) * dt ...
                            - L_y(j) * (S_yoS_y_inv + K_j * filter_cov(j) * K_j * dt) ...
                          );

    dy = y(j+1) - y(j);
    A_j = L_x(j) - S_xoS_y * S_yoS_y_inv * L_y(j);
    B_j = S_xoS_x - S_xoS_y * S_yoS_y_inv * S_yoS_x;


    % Update the posterior smoother mean and posterior smoother covariance using
    % the optimal nonlinear smoother state estimation backward equations for 
    % CGNSs; See Section 2.1.2 of the Supplementary Information.
    mu = muT - (L_x(j) * muT + f_x - B_j / filter_cov(j) * (filter_mean(j) - muT)) * dt ...
             + S_xoS_y * S_yoS_y_inv * (-dy + (L_y(j) * muT + f_y(j)) * dt);
    R = RT - ((A_j + B_j / filter_cov(j)) * RT + RT * (A_j + B_j / filter_cov(j)) - B_j) * dt;
    
    smoother_mean(j) = mu;
    smoother_cov(j) = R;
    muT = mu;
    RT = R;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN THIS CODE SECTION TO STUDY THE CAUSAL RELATIONSHIP: y(t) → x
% x (Predator): Observed | y (Prey): Unobserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Posterior smoother mean of the latent variable y.
smoother_mean = zeros(1, N+1); 
% Posterior smoother covariance matrix of the latent variable y.
smoother_cov = zeros(1, N+1);

% Smoother runs backwards: "Initial" values of the smoother statistics (i.e., at
% the last time instant T) are the corresponding posterior filter statistics.
smoother_mean(N+1) = filter_mean(N+1);
smoother_cov(N+1) = filter_cov(N+1);

% Auxiliary matrices used for the calculation of the online smoother for this 
% CGNS. The online smoother is required for the calculation of the subjective 
% and objective causal influence range (CIR) lengths of y(t) → x at each time 
% t∈[0,T]. Notation used is consistent with that of the original CGNS online 
% smoother work and the accompanying martingale-free introduction to CGNSs 
% paper: 
%   10.48550/arXiv.2411.05870 and 10.48550/arXiv.2410.24056
E_j_matrices = zeros(1, N+1);
F_j_matrices = zeros(1, N+1);
G_x_j = L_x(N+1) + S_xoS_y / filter_cov(N+1);
G_y_j = L_y(N+1) + S_yoS_y / filter_cov(N+1);
C_jj = 1 - G_y_j * dt;
H_j = filter_cov(N+1) \ (L_y(N+1) * filter_cov(N+1) + filter_cov(N+1) * L_y(N+1) + S_yoS_y);
K_j = S_xoS_x_inv * G_x_j;
E_j_matrices(N+1) = C_jj + S_yoS_x * K_j * dt;
F_j_matrices(N+1) = - filter_cov(N+1) * ( ...
                        K_j + (G_x_j * K_j * filter_cov(N+1) * K_j - filter_cov(N+1) \ H_j * filter_cov(N+1) * K_j + L_y(N+1) * K_j) * dt ...
                        - L_x(N+1) * (S_xoS_x_inv + K_j * filter_cov(N+1) * K_j * dt) ...
                      );

muT = smoother_mean(N+1);
RT = smoother_cov(N+1);

for j = N:-1:1

    % Calculation of the online smoother auxiliary matrices.
    G_x_j = L_x(j) + S_xoS_y / filter_cov(j);
    G_y_j = L_y(j) + S_yoS_y / filter_cov(j);
    C_jj = 1 - G_y_j * dt;
    H_j = filter_cov(j) \ (L_y(j) * filter_cov(j) + filter_cov(j) * L_y(j) + S_yoS_y);
    K_j = S_xoS_x_inv * G_x_j;
    E_j_matrices(j) = C_jj + S_yoS_x * K_j * dt;
    F_j_matrices(j) = - filter_cov(j) * ( ...
                        K_j + (G_x_j * K_j * filter_cov(j) * K_j - filter_cov(j) \ H_j * filter_cov(j) * K_j + L_y(j) * K_j) * dt ...
                        - L_x(j) * (S_xoS_x_inv + K_j * filter_cov(j) * K_j * dt) ...
                      );

    dx = x(j+1) - x(j);
    A_j = L_y(j) - S_yoS_x * S_xoS_x_inv * L_x(j);
    B_j = S_yoS_y - S_yoS_x * S_xoS_x_inv * S_xoS_y;

    % Update the posterior smoother mean and posterior smoother covariance using
    % the optimal nonlinear smoother state estimation backward equations for 
    % CGNSs; See Section 2.1.2 of the Supplementary Information.
    mu = muT - (L_y(j) * muT + f_y - B_j / filter_cov(j) * (filter_mean(j) - muT)) * dt ...
             + S_yoS_x * S_xoS_x_inv * (-dx + (L_x(j) * muT + f_x(j)) * dt);
    R = RT - ((A_j + B_j / filter_cov(j)) * RT + RT * (A_j + B_j / filter_cov(j)) - B_j) * dt;
    
    smoother_mean(j) = mu;
    smoother_cov(j) = R;
    muT = mu;
    RT = R;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  PLOTTING FILTER & SMOOTHER RESULTS  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN THIS CODE SECTION TO STUDY THE CAUSAL RELATIONSHIP: x(t) → y
% y (Prey): Observed || x (Predator): Unobserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting interval to be used throughout the run.
time_start_plot = 0;
time_end_plot = 50;

% Phase plot and time series of x and y.
figure('WindowState', 'maximized');

subplot(1, 3, 1)
plot(x(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), y(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'k', LineWidth=2)
title('Phase Plot of (x(t),y(t))')
xlabel('x (Predator)')
ylabel('y (Prey)')
grid on
fontsize(16, 'points')

subplot(1, 3, [2, 3])
plot(time_start_plot:dt:time_end_plot, x(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'm', LineWidth=2)
hold on
plot(time_start_plot:dt:time_end_plot, y(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'b', LineWidth=2)
yline(gamma/delta, 'k--', LineWidth=2)
title('Time Series of x (Predator) and y (Prey)')
xlabel('t')
legend('x', 'y', 'Antidamping Threshold for y: x < γ/δ', NumColumns=3)
grid on
fontsize(16, 'points')

% Time series of the posterior filter and smoother Gaussian statistics of x.
figure('WindowState', 'maximized');

subplot(2, 1, 1)
plot(time_start_plot:dt:time_end_plot, x(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'm', LineWidth=2)
hold on
plot(time_start_plot:dt:time_end_plot, filter_mean(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'g', LineWidth=2)
plot(time_start_plot:dt:time_end_plot, smoother_mean(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'r', LineWidth=2)
yline(gamma/delta, 'k--', LineWidth=2)
title('True and Posterior Mean Time Series of x')
xlabel('t')
ylabel('x')
legend('Truth', 'Filter', 'Smoother', 'Antidamping Threshold for y: x < γ/δ', NumColumns=4)
grid on
fontsize(16, 'points')

subplot(2, 1, 2)
plot(time_start_plot:dt:time_end_plot, filter_cov(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'g', LineWidth=2)
hold on
plot(time_start_plot:dt:time_end_plot, smoother_cov(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'r', LineWidth=2)
title('Posterior Variance of x')
xlabel('t')
ylabel('Var(x|y)')
legend('Filter', 'Smoother', NumColumns=2)
grid on
fontsize(16, 'points')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN THIS CODE SECTION TO STUDY THE CAUSAL RELATIONSHIP: y(t) → x
% x (Predator): Observed | y (Prey): Unobserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting interval to be used throughout the run.
time_start_plot = 0;
time_end_plot = 50;

% Phase plot and time series of x and y.
figure('WindowState', 'maximized');

subplot(1, 3, 1)
plot(x(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), y(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'k', LineWidth=2)
title('Phase Plot of (x(t),y(t))')
xlabel('x (Predator)')
ylabel('y (Prey)')
grid on
fontsize(16, 'points')

subplot(1, 3, [2, 3])
plot(time_start_plot:dt:time_end_plot, x(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'm', LineWidth=2)
hold on
plot(time_start_plot:dt:time_end_plot, y(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'b', LineWidth=2)
yline(alpha/beta, 'k--', LineWidth=2)
title('Time Series of x (Predator) and y (Prey)')
xlabel('t')
legend('x', 'y', 'Antidamping Threshold for x: y > α/β', NumColumns=3)
grid on
fontsize(16, 'points')

% Time series of the posterior filter and smoother Gaussian statistics of y.
figure('WindowState', 'maximized');

subplot(2, 1, 1)
plot(time_start_plot:dt:time_end_plot, y(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'b', LineWidth=2)
hold on
plot(time_start_plot:dt:time_end_plot, filter_mean(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'g', LineWidth=2)
plot(time_start_plot:dt:time_end_plot, smoother_mean(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'r', LineWidth=2)
yline(alpha/beta, 'k--', LineWidth=2)
title('True and Posterior Mean Time Series of y')
xlabel('t')
ylabel('y')
legend('Truth', 'Filter', 'Smoother', 'Antidamping Threshold for x: y > α/β', NumColumns=4)
grid on
fontsize(16, 'points')

subplot(2, 1, 2)
plot(time_start_plot:dt:time_end_plot, filter_cov(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'g', LineWidth=2)
hold on
plot(time_start_plot:dt:time_end_plot, smoother_cov(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'r', LineWidth=2)
title('Posterior Variance of y')
xlabel('t')
ylabel('Var(y|x)')
legend('Filter', 'Smoother', NumColumns=2)
grid on
fontsize(16, 'points')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ACI ANALYSIS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN THIS CODE SECTION TO STUDY THE CAUSAL RELATIONSHIP: x(t) → y
% y (Prey): Observed || x (Predator): Unobserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating the ACI metric for x(t) → y at each time t∈[0,T].
signal_smoother_filter = 0.5 * (smoother_mean - filter_mean).^2 ./ filter_cov;
cov_ratio_smoother_filter = smoother_cov ./ filter_cov;
dispersion_smoother_filter = 0.5 * (-log(cov_ratio_smoother_filter) + cov_ratio_smoother_filter - 1);
ACI_metric = signal_smoother_filter + dispersion_smoother_filter;

% Implementation of the fixed-lag online smoother for CGNSs: 
%   10.48550/arXiv.2411.05870
% The online smoother is required for the calculation of the subjective and 
% objective causal influence range (CIR) lengths of x(t) → y at each time 
% t∈[0,T]. Notation used is consistent with that of the original CGNS online 
% smoother work and the accompanying martingale-free introduction to CGNSs 
% paper: 
%   10.48550/arXiv.2411.05870 and 10.48550/arXiv.2410.24056

% The fixed-lag parameter is set equal to the total number of observations, 
% N = ⌈T/Δt⌉, such that at each time instant the full backward algorithm is 
% carried out. This is because each online smoother distribution, pₙ(xʲ|y), is 
% needed for the calculation of the subjective and objective CIRs; See Section 
% 2.3 of the Supplementary Information.
fixed_lag = N+1;

% Saving the online smoother mean, online smoother covariance matrices, and 
% update matrices in a cell array where each row is another cell array with as 
% many columns as the cardinal number of the current row. Using such nested cell
% arrays efficiently simulates staggered arrays in MATLAB. This approach 
% efficiently preserves space in memory without defining unnecessarily large 
% high-order tensors to store the online smoother estimations and update 
% matrices. In these nested cell arrays, the first/parent index corresponds to 
% n∈{j,j+1,...,N} for the current observation yⁿ, while the second/child index 
% corresponds to j∈{0,1,...,N} for the time instant tⱼ at which we carry out the
% online smoother state estimation for xʲ=x(tⱼ).
online_fixed_mean = cell(N+1, 1);
online_fixed_cov = cell(N+1, 1);
update_matrices_fixed = cell(N-1, 1);
for n = 1:(N-1)
    update_matrices_fixed{n} = zeros(1, n+1);
    online_fixed_mean{n} = zeros(1, n);
    online_fixed_cov{n} = zeros(1, n);
end
for n = N:N+1
    online_fixed_mean{n} = zeros(1, n);
    online_fixed_cov{n} = zeros(1, n);
end

% Details of the online smoother algorithm for CGNSs are briefly reviewed in 
% Section 2.2 of the Supplementary Information.

% Need to do the first two observations manually.

% A single observation (n=1).
online_fixed_mean{1}(1) = filter_mean(1);
online_fixed_cov{1}(1) = filter_cov(1);

% Two observations (n=2).
online_fixed_mean{2}(2) = filter_mean(2);
online_fixed_cov{2}(2) = filter_cov(2);
if fixed_lag == 0

    online_fixed_mean{2}(1) = online_fixed_mean{1}(1);
    online_fixed_cov{2}(1) = online_fixed_cov{1}(1);

else

    aux_vec = filter_mean(1) ...
              - E_j_matrices(1) * ((1 + L_x(1) * dt) * filter_mean(1) + f_x * dt) ...
              + F_j_matrices(1) * (y(2) - y(1) - (L_y(1) * filter_mean(1) + f_y(1)) * dt);
    online_fixed_mean{2}(1) = E_j_matrices(1) * filter_mean(2) + aux_vec;
    aux_mat = filter_cov(1) ...
              - E_j_matrices(1) * (1 + L_x(1) * dt) * filter_cov(1) ... 
              - F_j_matrices(1) * L_y(1) * filter_cov(1) * dt;
    online_fixed_cov{2}(1) = E_j_matrices(1) * filter_cov(2) * E_j_matrices(1) + aux_mat;

end

% Used for the text-based progress bar.
start_time = tic;

for n = 3:N+1

    % Text-based progress bar.
    progress_bar('Online Smoother Algorithm', n-2, length(3:N+1), start_time);

    online_fixed_mean{n}(n) = filter_mean(n);
    online_fixed_cov{n}(n) = filter_cov(n);

    if fixed_lag == 0

        online_fixed_mean{n}(n-1) = online_fixed_mean{n-1}(n-1);
        online_fixed_cov{n}(n-1) = online_fixed_cov{n-1}(n-1);

    else

        aux_vec = filter_mean(n-1) ...
                  - E_j_matrices(n-1) * ((1 + L_x(n-1) * dt) * filter_mean(n-1) + f_x * dt) ...
                  + F_j_matrices(n-1) * (y(n) - y(n-1) - (L_y(n-1) * filter_mean(n-1) + f_y(n-1)) * dt);
        online_fixed_mean{n}(n-1) = E_j_matrices(n-1) * filter_mean(n) + aux_vec;
        aux_mat = filter_cov(n-1) ...
                  - E_j_matrices(n-1) * (1 + L_x(n-1) * dt) * filter_cov(n-1) ...
                  - F_j_matrices(n-1) * L_y(n-1) * filter_cov(n-1) * dt;
        online_fixed_cov{n}(n-1) = E_j_matrices(n-1) * filter_cov(n) * E_j_matrices(n-1) + aux_mat;

    end

    for j = (n-1):-1:1

        if  (1 <= j) && (j <= n-1-fixed_lag)

            online_fixed_mean{n}(j) = online_fixed_mean{n-1}(j);
            online_fixed_cov{n}(j) = online_fixed_cov{n-1}(j);

        elseif (n-fixed_lag <= j) && (j <= n-1)

            if j == n-1
                update_matrices_fixed{n-2}(n-1) = 1;
            elseif j == n-2
                update_matrices_fixed{n-2}(n-2) = E_j_matrices(n-2);
            else
                update_matrices_fixed{n-2}(j) = update_matrices_fixed{n-3}(j) * E_j_matrices(n-2);
            end
            online_mean_inov = online_fixed_mean{n}(n-1) - filter_mean(n-1);
            online_fixed_mean{n}(j) = online_fixed_mean{n-1}(j) + update_matrices_fixed{n-2}(j) * online_mean_inov;
            online_cov_inov = online_fixed_cov{n}(n-1) - filter_cov(n-1);
            online_fixed_cov{n}(j) = online_fixed_cov{n-1}(j) + update_matrices_fixed{n-2}(j) * online_cov_inov * update_matrices_fixed{n-2}(j);

        end        
    end
end

% Calculating the subjective CIR length for x(t) → y at each time t∈[0,T] and 
% for various orders O(10⁻ᵏ) of ε values. The associated objective CIR length is 
% also calculated using its computationally efficient underestimating 
% approximation. The theory behind the subjective and objective CIR length is 
% given in Section 1.5 of the Supplementary Information, while their 
% computational details for CGNSs are given in Section 2.3. 

% Letting 10⁻⁶ ≤ ε ≤ 10⁰‧⁵ with a resolution of 513 points.
epsilon_resolution = 513;
lowest_order = -6;
highest_order = 0.5;
eps_ord_values = flip(linspace(lowest_order, highest_order, epsilon_resolution));

% The following snippet instead uses a more adaptive mesh of ε, instead of a 
% purely logarithmic one as the above implementation, for a more realistic 
% integration of the subjective CIRs over ε for obtaining the exact 
% corresponding objective CIR.
% adaptive_point = -2;
% half_resolution = 250;
% eps_ord_values = unique([
%     flip(log10(linspace(10^adaptive_point, 10^highest_order, half_resolution))), ...
%     flip(linspace(lowest_order, adaptive_point, half_resolution))
% ], 'stable');

% Calculating the subjective and objective CIRs over the plotting time interval 
% of choice. We add a lookahead tolerance for the lagged observational time: 
%   T'∈[t,time_end_plot+lookahead_tolerance],
% to avoid observational saturation as t approaches time_end_plot.
lookahead_tolerance = 5;
first_idx = round(time_start_plot/dt)+1;
if time_end_plot+lookahead_tolerance < T
    last_idx = round((time_end_plot+lookahead_tolerance)/dt)+1;
else
    last_idx = round(time_end_plot/dt)+1;
end
plot_len = length(first_idx:last_idx);

subjective_CIR = zeros(length(eps_ord_values), plot_len);
approx_objective_CIR = zeros(1, plot_len);
% CIR relative entropy metric δ(T';t) (See Section 1.5.1 of the Supplementary 
% Information) used to calculate the approximate objective CIR via a time 
% integral over the lagged observational time T' instead of integrating 
% the associated subjective CIR over ε as in the definition to get the exact 
% objective CIR length; See Section 1.5.3 of the Supplementary Information. 
% Using the notation from Section 2.3 of the Supplementary Information, 
% RE_metric is Pₙʲ, with the rows of RE_metric corresponding to the natural time 
% t (i.e., j∈{first_idx,first_idx+1,...,last_idx} index) while the columns 
% correspond to the lagged observational time T' (i.e., n∈{j,j+1,...,last_idx} 
% index).
RE_metric = zeros(plot_len, plot_len); 
max_RE_metric = zeros(1, plot_len);

% Used for the text-based progress bar.
start_time = tic;

for eps_idx = 1:length(eps_ord_values)

    epsilon = 10^eps_ord_values(eps_idx);
    
    for j = first_idx:last_idx

        % Text-based progress bar.
        progress_bar('Calculation of the CIRs', j-first_idx+1+(eps_idx-1)*plot_len, length(eps_ord_values)*plot_len, start_time);
        
        % Calculation of the objective CIR length approximation.
        if eps_idx == 1

            % Calculating and storing Pₙʲ over n∈{j,j+1,...,last_idx} for a 
            % fixed j∈{first_idx,first_idx+1,...,last_idx}.
            RE_n = zeros(1, length(j:last_idx));
            for obs = j:last_idx
                cov_ratio = online_fixed_cov{end}(j) / online_fixed_cov{obs}(j);
                RE_n(obs-j+1) = 0.5 * (online_fixed_mean{end}(j) - online_fixed_mean{obs}(j))^2 / online_fixed_cov{obs}(j) ...
                                + 0.5 * (cov_ratio - 1 - log(cov_ratio));
            end
            max_RE_metric(j-first_idx+1) = max(RE_n);
            RE_metric(j-first_idx+1, 1:length(RE_n)) = RE_n;
            
            % If it is essentially zero then do not calculate the CIR and set 
            % equal to 0 instead (as to avoid operationally inflated CIR values
            % due to the normalization in the objective CIR approximation and
            % to numerical precision errors).
            RE_metric_threshold = 1e-5;
            if max(RE_n) > RE_metric_threshold

                % Estimation of the integral using a composite trapezoidal
                % rule.
                % approx_objective_CIR(j-first_idx+1) = trapz(RE_n)*dt/max(RE_n);

                % Estimation of the integral using a composite Simpson's 1/3 
                % rule for better accuracy.
                try
                    approx_objective_CIR(j-first_idx+1) = simps(RE_n)*dt/max(RE_n);
                catch
                    approx_objective_CIR(j-first_idx+1) = 0;
                end

            else

                approx_objective_CIR(j-first_idx+1) = 0;
            
            end

        end
        
        % Calculation of the subjective CIR length for this ε value.
        RE_n = RE_metric(j-first_idx+1, 1:length(j:last_idx));
        subj_CIR_idx = find(RE_n > epsilon, 1, 'last');
        if isempty(subj_CIR_idx)
            subj_CIR_idx = 0;
        end
        subjective_CIR(eps_idx, j-first_idx+1) = subj_CIR_idx*dt;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  NOTE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% When choosing a sufficiently large epsilon_resolution and a sufficiently large 
% k where 10⁻ᵏ ≤ ε, then the objective CIR at each time tⱼ can be calculated 
% through the definition by averaging the corresponding subjective CIR over ε 
% via a numerical quadrature method using the following command (the flip 
% operations are needed to put the ε (see eps_ord_values variable) interval in 
% ascending order):

% Estimation of the integral using a composite trapezoidal rule.
% defn_objective_CIR = trapz(10.^flip(eps_ord_values), flipud(subjective_CIR(:, 1:end-lookahead_tolerance/dt)), 1)./max_RE_metric(1:end-lookahead_tolerance/dt);

% Estimation of the integral using a composite Simpson's 1/3 rule for better 
% accuracy.
defn_objective_CIR = simps(10.^flip(eps_ord_values), flipud(subjective_CIR(:, 1:end-lookahead_tolerance/dt)), 1)./max_RE_metric(1:end-lookahead_tolerance/dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN THIS CODE SECTION TO STUDY THE CAUSAL RELATIONSHIP: y(t) → x
% x (Predator): Observed | y (Prey): Unobserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating the ACI metric for y(t) → x at each time t∈[0,T].
signal_smoother_filter = 0.5 * (smoother_mean - filter_mean).^2 ./ filter_cov;
cov_ratio_smoother_filter = smoother_cov ./ filter_cov;
dispersion_smoother_filter = 0.5 * (-log(cov_ratio_smoother_filter) + cov_ratio_smoother_filter - 1);
ACI_metric = signal_smoother_filter + dispersion_smoother_filter;

% Implementation of the fixed-lag online smoother for CGNSs: 
%   10.48550/arXiv.2411.05870
% The online smoother is required for the calculation of the subjective and 
% objective causal influence range (CIR) lengths of y(t) → x at each time 
% t∈[0,T]. Notation used is consistent with that of the original CGNS online 
% smoother work and the accompanying martingale-free introduction to CGNSs 
% paper: 
%   10.48550/arXiv.2411.05870 and 10.48550/arXiv.2410.24056

% The fixed-lag parameter is set equal to the total number of observations, 
% N = ⌈T/Δt⌉, such that at each time instant the full backward algorithm is 
% carried out. This is because each online smoother distribution, pₙ(yʲ|x), is 
% needed for the calculation of the subjective and objective CIRs; See Section 
% 2.3 of the Supplementary Information.
fixed_lag = N+1;

% Saving the online smoother mean, online smoother covariance matrices, and 
% update matrices in a cell array where each row is another cell array with as 
% many columns as the cardinal number of the current row. Using such nested cell
% arrays efficiently simulates staggered arrays in MATLAB. This approach 
% efficiently preserves space in memory without defining unnecessarily large 
% high-order tensors to store the online smoother estimations and update 
% matrices. In these nested cell arrays, the first/parent index corresponds to 
% n∈{j,j+1,...,N} for the current observation xⁿ, while the second/child index 
% corresponds to j∈{0,1,...,N} for the time instant tⱼ at which we carry out the
% online smoother state estimation for yʲ=y(tⱼ).
online_fixed_mean = cell(N+1, 1);
online_fixed_cov = cell(N+1, 1);
update_matrices_fixed = cell(N-1, 1);
for n = 1:(N-1)
    update_matrices_fixed{n} = zeros(1, n+1);
    online_fixed_mean{n} = zeros(1, n);
    online_fixed_cov{n} = zeros(1, n);
end
for n = N:N+1
    online_fixed_mean{n} = zeros(1, n);
    online_fixed_cov{n} = zeros(1, n);
end

% Details of the online smoother algorithm for CGNSs are briefly reviewed in 
% Section 2.2 of the Supplementary Information.

% Need to do the first two observations manually.

% A single observation (n=1).
online_fixed_mean{1}(1) = filter_mean(1);
online_fixed_cov{1}(1) = filter_cov(1);

% Two observations (n=2).
online_fixed_mean{2}(2) = filter_mean(2);
online_fixed_cov{2}(2) = filter_cov(2);
if fixed_lag == 0

    online_fixed_mean{2}(1) = online_fixed_mean{1}(1);
    online_fixed_cov{2}(1) = online_fixed_cov{1}(1);

else

    aux_vec = filter_mean(1) ...
              - E_j_matrices(1) * ((1 + L_y(1) * dt) * filter_mean(1) + f_y * dt) ...
              + F_j_matrices(1) * (x(2) - x(1) - (L_x(1) * filter_mean(1) + f_x(1)) * dt);
    online_fixed_mean{2}(1) = E_j_matrices(1) * filter_mean(2) + aux_vec;
    aux_mat = filter_cov(1) ...
              - E_j_matrices(1) * (1 + L_y(1) * dt) * filter_cov(1) ... 
              - F_j_matrices(1) * L_x(1) * filter_cov(1) * dt;
    online_fixed_cov{2}(1) = E_j_matrices(1) * filter_cov(2) * E_j_matrices(1) + aux_mat;

end

% Used for the text-based progress bar.
start_time = tic;

for n = 3:N+1

    % Text-based progress bar.
    progress_bar('Online Smoother Algorithm', n-2, length(3:N+1), start_time);

    online_fixed_mean{n}(n) = filter_mean(n);
    online_fixed_cov{n}(n) = filter_cov(n);

    if fixed_lag == 0

        online_fixed_mean{n}(n-1) = online_fixed_mean{n-1}(n-1);
        online_fixed_cov{n}(n-1) = online_fixed_cov{n-1}(n-1);

    else

        aux_vec = filter_mean(n-1) ...
                  - E_j_matrices(n-1) * ((1 + L_y(n-1) * dt) * filter_mean(n-1) + f_y * dt) ...
                  + F_j_matrices(n-1) * (x(n) - x(n-1) - (L_x(n-1) * filter_mean(n-1) + f_x(n-1)) * dt);
        online_fixed_mean{n}(n-1) = E_j_matrices(n-1) * filter_mean(n) + aux_vec;
        aux_mat = filter_cov(n-1) ...
                  - E_j_matrices(n-1) * (1 + L_y(n-1) * dt) * filter_cov(n-1) ...
                  - F_j_matrices(n-1) * L_x(n-1) * filter_cov(n-1) * dt;
        online_fixed_cov{n}(n-1) = E_j_matrices(n-1) * filter_cov(n) * E_j_matrices(n-1) + aux_mat;

    end

    for j = (n-1):-1:1

        if  (1 <= j) && (j <= n-1-fixed_lag)

            online_fixed_mean{n}(j) = online_fixed_mean{n-1}(j);
            online_fixed_cov{n}(j) = online_fixed_cov{n-1}(j);

        elseif (n-fixed_lag <= j) && (j <= n-1)

            if j == n-1
                update_matrices_fixed{n-2}(n-1) = 1;
            elseif j == n-2
                update_matrices_fixed{n-2}(n-2) = E_j_matrices(n-2);
            else
                update_matrices_fixed{n-2}(j) = update_matrices_fixed{n-3}(j) * E_j_matrices(n-2);
            end
            online_mean_inov = online_fixed_mean{n}(n-1) - filter_mean(n-1);
            online_fixed_mean{n}(j) = online_fixed_mean{n-1}(j) + update_matrices_fixed{n-2}(j) * online_mean_inov;
            online_cov_inov = online_fixed_cov{n}(n-1) - filter_cov(n-1);
            online_fixed_cov{n}(j) = online_fixed_cov{n-1}(j) + update_matrices_fixed{n-2}(j) * online_cov_inov * update_matrices_fixed{n-2}(j);

        end        
    end
end

% Calculating the subjective CIR length for y(t) → x at each time t∈[0,T] and 
% for various orders O(10⁻ᵏ) of ε values. The associated objective CIR length is 
% also calculated using its computationally efficient underestimating 
% approximation. The theory behind the subjective and objective CIR length is 
% given in Section 1.5 of the Supplementary Information, while their 
% computational details for CGNSs are given in Section 2.3. 

% Letting 10⁻⁶ ≤ ε ≤ 10⁰‧⁵ with a resolution of 513 points.
epsilon_resolution = 513;
lowest_order = -6;
highest_order = 0.5;
eps_ord_values = flip(linspace(lowest_order, highest_order, epsilon_resolution));

% The following snippet instead uses a more adaptive mesh of ε, instead of a 
% purely logarithmic one as the above implementation, for a more realistic 
% integration of the subjective CIRs over ε for obtaining the exact 
% corresponding objective CIR.
% adaptive_point = -2;
% half_resolution = 250;
% eps_ord_values = unique([
%     flip(log10(linspace(10^adaptive_point, 10^highest_order, half_resolution))), ...
%     flip(linspace(lowest_order, adaptive_point, half_resolution))
% ], 'stable');

% Calculating the subjective and objective CIRs over the plotting time interval 
% of choice. We add a lookahead tolerance for the lagged observational time: 
%   T'∈[t,time_end_plot+lookahead_tolerance],
% to avoid observational saturation as t approaches time_end_plot.
lookahead_tolerance = 5;
first_idx = round(time_start_plot/dt)+1;
if time_end_plot+lookahead_tolerance < T
    last_idx = round((time_end_plot+lookahead_tolerance)/dt)+1;
else
    last_idx = round(time_end_plot/dt)+1;
end
plot_len = length(first_idx:last_idx);

subjective_CIR = zeros(length(eps_ord_values), plot_len);
approx_objective_CIR = zeros(1, plot_len);
% CIR relative entropy metric δ(T';t) (See Section 1.5.1 of the Supplementary 
% Information) used to calculate the approximate objective CIR via a time 
% integral over the lagged observational time T' instead of integrating 
% the associated subjective CIR over ε as in the definition to get the exact 
% objective CIR length; See Section 1.5.3 of the Supplementary Information. 
% Using the notation from Section 2.3 of the Supplementary Information, 
% RE_metric is Pₙʲ, with the rows of RE_metric corresponding to the natural time 
% t (i.e., j∈{first_idx,first_idx+1,...,last_idx} index) while the columns 
% correspond to the lagged observational time T' (i.e., n∈{j,j+1,...,last_idx} 
% index).
RE_metric = zeros(plot_len, plot_len); 
max_RE_metric = zeros(1, plot_len);

% Used for the text-based progress bar.
start_time = tic;

for eps_idx = 1:length(eps_ord_values)

    epsilon = 10^eps_ord_values(eps_idx);
    
    for j = first_idx:last_idx

        % Text-based progress bar.
        progress_bar('Calculation of the CIRs', j-first_idx+1+(eps_idx-1)*plot_len, length(eps_ord_values)*plot_len, start_time);
        
        % Calculation of the objective CIR length approximation.
        if eps_idx == 1

            % Calculating and storing Pₙʲ over n∈{j,j+1,...,last_idx} for a 
            % fixed j∈{first_idx,first_idx+1,...,last_idx}.
            RE_n = zeros(1, length(j:last_idx));
            for obs = j:last_idx
                cov_ratio = online_fixed_cov{end}(j) / online_fixed_cov{obs}(j);
                RE_n(obs-j+1) = 0.5 * (online_fixed_mean{end}(j) - online_fixed_mean{obs}(j))^2 / online_fixed_cov{obs}(j) ...
                                + 0.5 * (cov_ratio - 1 - log(cov_ratio));
            end
            max_RE_metric(j-first_idx+1) = max(RE_n);
            RE_metric(j-first_idx+1, 1:length(RE_n)) = RE_n;
            
            % If it is essentially zero then do not calculate the CIR and set 
            % equal to 0 instead (as to avoid operationally inflated CIR values
            % due to the normalization in the objective CIR approximation and
            % to numerical precision errors).
            RE_metric_threshold = 1e-5;
            if max(RE_n) > RE_metric_threshold

                % Estimation of the integral using a composite trapezoidal
                % rule.
                % approx_objective_CIR(j-first_idx+1) = trapz(RE_n)*dt/max(RE_n);

                % Estimation of the integral using a composite Simpson's 1/3 
                % rule for better accuracy.
                try
                    approx_objective_CIR(j-first_idx+1) = simps(RE_n)*dt/max(RE_n);
                catch
                    approx_objective_CIR(j-first_idx+1) = 0;
                end

            else

                approx_objective_CIR(j-first_idx+1) = 0;
            
            end

        end
        
        % Calculation of the subjective CIR length for this ε value.
        RE_n = RE_metric(j-first_idx+1, 1:length(j:last_idx));
        subj_CIR_idx = find(RE_n > epsilon, 1, 'last');
        if isempty(subj_CIR_idx)
            subj_CIR_idx = 0;
        end
        subjective_CIR(eps_idx, j-first_idx+1) = subj_CIR_idx*dt;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  NOTE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% When choosing a sufficiently large epsilon_resolution and a sufficiently large 
% k where 10⁻ᵏ ≤ ε, then the objective CIR at each time tⱼ can be calculated 
% through the definition by averaging the corresponding subjective CIR over ε 
% via a numerical quadrature method using the following command (the flip 
% operations are needed to put the ε (see eps_ord_values variable) interval in 
% ascending order):

% Estimation of the integral using a composite trapezoidal rule.
% defn_objective_CIR = trapz(10.^flip(eps_ord_values), flipud(subjective_CIR(:, 1:end-lookahead_tolerance/dt)), 1)./max_RE_metric(1:end-lookahead_tolerance/dt);

% Estimation of the integral using a composite Simpson's 1/3 rule for better 
% accuracy.
defn_objective_CIR = simps(10.^flip(eps_ord_values), flipud(subjective_CIR(:, 1:end-lookahead_tolerance/dt)), 1)./max_RE_metric(1:end-lookahead_tolerance/dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  PLOTTING ACI ANALYSIS RESULTS  %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN THIS CODE SECTION TO STUDY THE CAUSAL RELATIONSHIP: x(t) → y
% y (Prey): Observed || x (Predator): Unobserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time series of the ACI metric and objective CIR length for x(t) → y, as well 
% as heatmap of the subjective CIR length over time and ε.
figure('WindowState', 'maximized');

subplot(3, 1, 1)
plot(time_start_plot:dt:time_end_plot, ACI_metric(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'k', LineWidth=2)
title('ACI Metric for x(t) \rightarrow y')
xlabel('t')
ylabel('Relative Entropy')
grid on
fontsize(16, 'points')

subplot(3, 1, 2)
xx = time_start_plot:dt:time_end_plot;
yy = eps_ord_values;
[X, Y] = meshgrid(xx, yy);
% pcolor(X, Y, log10(subjective_CIR(:, 1:end-lookahead_tolerance/dt)./max_RE_metric(1:end-lookahead_tolerance/dt)));
pcolor(X, Y, log10(subjective_CIR(:, 1:end-lookahead_tolerance/dt)));
shading interp
colormap("jet")
clim([-1, 0.5])
cb = colorbar('eastoutside', Ticks=-1:0.5:0.5, TickLabels=arrayfun(@(x) sprintf('10^{%.1f}', x), -1:0.5:0.5, 'UniformOutput', false));
cb.Position(1) = cb.Position(1) + 0.1;
cb.Position(3) = 0.015;
ylim([lowest_order, highest_order])
set(gca, 'YDir','reverse')
title('Subjective CIR Length for x(t) \rightarrow y (Logarithmic Scale)')
xlabel('t')
ylabel('ε')
yticks(lowest_order:highest_order)
yticklabels(arrayfun(@(x) sprintf('10^{%d}', x), lowest_order:highest_order, 'UniformOutput', false))
fontsize(16, 'points')

subplot(3, 1, 3)
plot(time_start_plot:dt:time_end_plot, x(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'm', LineWidth=2)
hold on
plot(time_start_plot:dt:time_end_plot, y(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'b', LineWidth=2)
yline(gamma/delta, 'k--', LineWidth=2)
title('True and Posterior Mean Time Series of x')
% Plotting the objective CIR length whiskers extending forward in time from each 
% x(t) at uniform intervals for optimal plotting.
% Whisker-plotting interval measured in indices.
whisker_interval = 40; 
for j = round(time_start_plot/dt)+1:whisker_interval:round(time_end_plot/dt)+1
    plot([(j-1)*dt, (j-1)*dt+approx_objective_CIR(j-round(time_start_plot/dt))], [x(j), x(j)], 'color', [0.5410, 0.2470, 0, 0.4], 'LineWidth', 0.5)
    plot((j-1)*dt+approx_objective_CIR(j-round(time_start_plot/dt)), x(j), 'color', [0.5410, 0.2470, 0], 'marker', '|')
end
xlim([time_start_plot, time_end_plot])
title('Time Series and Objective CIR Length for x(t) \rightarrow y')
xlabel('t')
legend('x', 'y', 'Antidamping Threshold for y: x < γ/δ', 'Objective CIR Length', NumColumns=4)
set(gca,'XGrid', 'on', 'YGrid', 'off')
fontsize(16, 'points')

% Heatmap of the normalized CIR relative entropy metric for x(t) → y: 
%   δ(T';t)/max{δ(T';t): T'∈[t,T]},
% over natural time t∈[0,T] and lagged observational time after t, T'-t∈[0,T],
% and comparison of the time series of the objective CIR length values for 
% x(t) → y calculated using the definition and computationally efficient 
% underestimate approximation.
figure('WindowState', 'maximized');

subplot(3, 2, 1:4)
% Choosing a smaller lagged observational time window for optimal plotting:
%   lag_obs_time_end_plot∈[0,time_end_plot-time_start_plot]
% Note that T'∈[t,T] at each t∈[0,T].
lag_obs_time_end_plot = 10;
data = log10(RE_metric(1:end-lookahead_tolerance/dt, 1:end-lookahead_tolerance/dt).'./max_RE_metric(1:end-lookahead_tolerance/dt));
data(data == -Inf) = NaN;
h = imagesc([time_start_plot time_end_plot], [0, lag_obs_time_end_plot], data(1:round(lag_obs_time_end_plot/dt)+1,:));
set(h, 'AlphaData', ~isnan(data(1:round(lag_obs_time_end_plot/dt)+1,:)))
colormap("jet")
clim([-10, 0])
cb = colorbar('eastoutside', Ticks=-10:0, TickLabels=arrayfun(@(x) sprintf('10^{%.1f}', x), -10:0, 'UniformOutput', false));
cb.Position(1) = cb.Position(1) + 0.1;
cb.Position(3) = 0.015;
set(gca, 'Color', [0.7, 0.7, 0.7])
set(gca, 'YDir','reverse')
ylim([0, lag_obs_time_end_plot])
title("Normalized CIR Relative Entropy Metric for x(t) \rightarrow y: δ(T';t)/max_{T'\in[t,T]}\{δ(T';t)\} (Logarithmic Scale)")
xlabel('t (Natural Time)')
ylabel("T'-t (Lagged Observational Time After t)")
fontsize(16, 'points')

subplot(3, 2, [5, 6])
plot(time_start_plot:dt:time_end_plot, approx_objective_CIR(1:end-lookahead_tolerance/dt), 'color', [0.5410, 0.2470, 0], LineWidth=2)
hold on
plot(time_start_plot:dt:time_end_plot, defn_objective_CIR, 'r', LineWidth=2)
title('Comparison of the Objective CIR Length Algorithms for x(t) \rightarrow y: Definition vs Approximation')
xlabel('t')
ylabel('Objective CIR Length')
legend('Efficient Approximation', 'Definition', NumColumns=2)
grid on
fontsize(16, 'points')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN THIS CODE SECTION TO STUDY THE CAUSAL RELATIONSHIP: y(t) → x
% x (Predator): Observed | y (Prey): Unobserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time series of the ACI metric and objective CIR length for y(t) → x, as well 
% as heatmap of the subjective CIR length over time and ε.
figure('WindowState', 'maximized');

subplot(3, 1, 1)
plot(time_start_plot:dt:time_end_plot, ACI_metric(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'k', LineWidth=2)
title('ACI Metric for y(t) \rightarrow x')
xlabel('t')
ylabel('Relative Entropy')
grid on
fontsize(16, 'points')

subplot(3, 1, 2)
xx = time_start_plot:dt:time_end_plot;
yy = eps_ord_values;
[X, Y] = meshgrid(xx, yy);
% pcolor(X, Y, log10(subjective_CIR(:, 1:end-lookahead_tolerance/dt)./max_RE_metric(1:end-lookahead_tolerance/dt)));
pcolor(X, Y, log10(subjective_CIR(:, 1:end-lookahead_tolerance/dt)));
shading interp
colormap("jet")
clim([-1, 1.25])
cb = colorbar('eastoutside', Ticks=-1:0.5:1.25, TickLabels=arrayfun(@(x) sprintf('10^{%.1f}', x), -1:0.5:1.25, 'UniformOutput', false));
cb.Position(1) = cb.Position(1) + 0.1;
cb.Position(3) = 0.015;
ylim([lowest_order, highest_order])
set(gca, 'YDir','reverse')
title('Subjective CIR Length for x(t) \rightarrow y (Logarithmic Scale)')
xlabel('t')
ylabel('ε')
yticks(lowest_order:highest_order)
yticklabels(arrayfun(@(x) sprintf('10^{%d}', x), lowest_order:highest_order, 'UniformOutput', false))
fontsize(16, 'points')

subplot(3, 1, 3)
plot(time_start_plot:dt:time_end_plot, x(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'm', LineWidth=2)
hold on
plot(time_start_plot:dt:time_end_plot, y(round(time_start_plot/dt)+1:round(time_end_plot/dt)+1), 'b', LineWidth=2)
yline(alpha/beta, 'k--', LineWidth=2)
% Plotting the objective CIR length whiskers extending forward in time from each 
% y(t) at uniform intervals for optimal plotting.
% Whisker-plotting interval measured in indices.
whisker_interval = 40; 
for j = round(time_start_plot/dt)+1:whisker_interval:round(time_end_plot/dt)+1
    plot([(j-1)*dt, (j-1)*dt+approx_objective_CIR(j-round(time_start_plot/dt))], [y(j), y(j)], 'color', [0, 0.2470, 0.5410, 0.4], 'LineWidth', 0.5)
    plot((j-1)*dt+approx_objective_CIR(j-round(time_start_plot/dt)), y(j), 'color', [0 0.2470 0.5410], 'marker', '|')
end
xlim([time_start_plot, time_end_plot])
title('Time Series and Objective CIR Length for y(t) \rightarrow x')
xlabel('t')
legend('x', 'y', 'Antidamping Threshold for x: y > α/β', 'Objective CIR Length', NumColumns=4)
set(gca,'XGrid', 'on', 'YGrid', 'off')
fontsize(16, 'points')

% Heatmap of the normalized CIR relative entropy metric for y(t) → x: 
%   δ(T';t)/max{δ(T';t): T'∈[t,T]},
% over natural time t∈[0,T] and lagged observational time after t, T'-t∈[0,T],
% and comparison of the time series of the objective CIR length values for 
% y(t) → x calculated using the definition and computationally efficient 
% underestimate approximation.
figure('WindowState', 'maximized');

subplot(3, 2, 1:4)
% Choosing a smaller lagged observational time window for optimal plotting:
%   lag_obs_time_end_plot∈[0,time_end_plot-time_start_plot]
% Note that T'∈[t,T] at each t∈[0,T].
lag_obs_time_end_plot = 20;
data = log10(RE_metric(1:end-lookahead_tolerance/dt, 1:end-lookahead_tolerance/dt).'./max_RE_metric(1:end-lookahead_tolerance/dt));
data(data == -Inf) = NaN;
h = imagesc([time_start_plot time_end_plot], [0, lag_obs_time_end_plot], data(1:round(lag_obs_time_end_plot/dt)+1,:));
set(h, 'AlphaData', ~isnan(data(1:round(lag_obs_time_end_plot/dt)+1,:)))
colormap("jet")
clim([-10, 0])
cb = colorbar('eastoutside', Ticks=-10:0, TickLabels=arrayfun(@(x) sprintf('10^{%.1f}', x), -10:0, 'UniformOutput', false));
cb.Position(1) = cb.Position(1) + 0.1;
cb.Position(3) = 0.015;
set(gca, 'Color', [0.7, 0.7, 0.7])
set(gca, 'YDir','reverse')
ylim([0, lag_obs_time_end_plot])
title("Normalized CIR Relative Entropy Metric for y(t) \rightarrow x: δ(T';t)/max_{T'\in[t,T]}\{δ(T';t)\} (Logarithmic Scale)")
xlabel('t (Natural Time)')
ylabel("T'-t (Lagged Observational Time After t)")
fontsize(16, 'points')

subplot(3, 2, [5, 6])
plot(time_start_plot:dt:time_end_plot, approx_objective_CIR(1:end-lookahead_tolerance/dt), 'color', [0 0.2470 0.5410], LineWidth=2)
hold on
plot(time_start_plot:dt:time_end_plot, defn_objective_CIR, 'r', LineWidth=2)
title('Comparison of the Objective CIR Length Algorithms for y(t) \rightarrow x: Definition vs Approximation')
xlabel('t')
ylabel('Objective CIR Length')
legend('Efficient Approximation', 'Definition', NumColumns=2)
grid on
fontsize(16, 'points')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%