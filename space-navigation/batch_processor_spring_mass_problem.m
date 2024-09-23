clc;
clf;
close all;
clear all;

% Observation Vector
%   range               range-rate
%   [m]                 [m/s]
Y = [
    6.1773780845922     0;
    5.56327661282686    1.31285863495514;
    5.69420161397342    -1.54488114381612;
    6.15294262127432    0.534923988815733;
    5.46251322092491    0.884698415328368;
    5.83638064328625    -1.56123248918054;
    6.08236452736002    1.00979943157547;
    5.40737619817037    0.31705117039215;
    5.97065615746125    -1.37453070975606;
    5.97369258835895    1.36768169443236;
    5.40669060248179    -0.302111588503166
];
tspan = 0:10; %Time span [s]

% Constants
k1 = 2.5; %Spring constant 1 [N/m]
k2 = 3.7; %Spring constant 2 [N/m]
m = 1.5; %Mass [kg]
h = 5.4; %Height between mass and observation point [m]
x0 = 3.0; %Initial position [m]
v0 = 0.0; %Initial velocity [m/s]
omega = sqrt((k1 + k2)/m); %Frequency [1/s]

% A Priori Information
Xstar = [4.0; 0.2]; %Initial State Vector [m; m/s]
P0 = [1000 0; 0 100]; %Covariance Matrix [m -; - m/s]
deltaXbar = [0; 0]; %State Vector Deviation [m; m/s]

fprintf("Batch Processing Algorithm Results\n")
fprintf("%s",repelem("-",1,124))
fprintf("\n| Iteration | xhat0 [m] | vhat0 [m/s] | Information Matrix Condition Number | Range Residual RMS | Range-Rate Residual RMS |\n")

% Batch Processing Algorithm
j = 1;
while j<=4
    % Initialize the iteration
    Lambda = pinv(P0);
    Nu = pinv(P0)*deltaXbar;
    epsilon = [];
    for t=tspan
        % Read the next observation
        rho = Y((t+1),1);
        rhodot = Y((t+1),2);

        % Integrate the reference trajectory and state transition matrix
        if t == 0
            X = [Xstar; 1; 0; 0; 1];
        else
            [~, X] = ode89(@odefunc,[t-1 t],[x; v; STM(1,1); STM(2,1); STM(1,2); STM(2,2)]);
            X = X(end,:);
        end
        x = X(1);
        v = X(2);
        STM = reshape(X(3:end),2,2);

        % Accumulate current observations
        Htilde = [
            x/sqrt(x^2 + h^2) 0;
            (v/sqrt(x^2 + h^2) - x^2*v/sqrt(x^2 + h^2)^3) x/rho
        ];
        y = [rho;rhodot] - [sqrt(x^2 + h^2); x*v/sqrt(x^2 + h^2)];
        H = Htilde*STM;
        Lambda = Lambda + transpose(H)*H;
        Nu = Nu + transpose(H)*y;
        epsilon = [epsilon; transpose(y)];
    end

    % Solve normal equations
    deltaXhat = pinv(Lambda)*(Nu);

    % Update values for next iteration
    Xstar = Xstar + deltaXhat;
    deltaXbar = deltaXbar - deltaXhat;

    % Calculate RMS values
    rho_RMS = sqrt(sum(transpose(epsilon(:,1))*epsilon(:,1))/11);
    rhodot_RMS = sqrt(sum(transpose(epsilon(:,2))*epsilon(:,2))/11);

    % Plot the residuals
    figure()
    s = stackedplot(tspan,epsilon);
    s.FontSize = 12;
    s.Title = sprintf("Residuals, Iteration %i", j);
    s.DisplayLabels = ["Range [m]", "Range Rate [m/s]"];
    s.XLabel = "Time [s]";
    s.Marker = "x";
    s.LineStyle = "--";

    % Output table row
    fprintf("| %9i | %9.5f | %11.5e | %35.5f | %18.2e | %23.2e |\n", j, Xstar(1), Xstar(2), cond(Lambda), rho_RMS, rhodot_RMS)
    j = j + 1;
end
fprintf("%s",repelem("-",1,124))
fprintf("\nBest Estimate of x0 = %.5f m\n", Xstar(1))
fprintf("Best Estimate of v0 = %.5e m/s\n", Xstar(2))
P0 = pinv(Lambda);
sigma_x = sqrt(P0(1,1));
sigma_v = sqrt(P0(2,2));
mu_xv = P0(2,1);
rho_xv = mu_xv/(sigma_x*sigma_v);
fprintf("Standard Deviation of x0 = %.3f\n", sigma_x)
fprintf("Standard Deviation of v0 = %.3f\n", sigma_v)
fprintf("Correlation Coefficient = %.4f\n", rho_xv)

function dXdt = odefunc(~,X)
% Constants
k1 = 2.5; %Spring constant 1 [N/m]
k2 = 3.7; %Spring constant 2 [N/m]
m = 1.5; %Mass [kg]
omega = sqrt((k1 + k2)/m); %Frequency [1/s]
A = [0 1; -omega^2 0];

% Seperate input vector
x = X(1);
v = X(2);
phi = reshape(X(3:end),2,2);
dphi = A*phi;

% Set derivative values
dXdt = zeros(6,1);
dXdt(1) = v;
dXdt(2) = -omega^2*x;
dXdt(3:end) = reshape(dphi,4,1);
end

