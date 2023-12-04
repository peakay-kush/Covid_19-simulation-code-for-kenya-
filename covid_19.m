% Parameters
beta = 0.3;   % infection rate
gamma = 0.1;  % recovery rate
sigma = 0.1;  % incubation rate

% Initial conditions
I0 = 110;        % initial number of infected individuals
E0 = 0;        % initial number of exposed individuals
R0 = 2;        % initial number of recovered individuals
S0 = 53010000;     % initial number of susceptible individuals (approximate population of Kenya)

% Time vector (days)
t = linspace(0, 365*5, 365*5);

% SEIR model differential equations
ode = @(t, y) [ -beta * y(1) * y(3) / sum(y);           % dS/dt
                 beta * y(1) * y(3) / sum(y) - sigma * y(2);  % dE/dt
                 sigma * y(2) - gamma * y(3);                % dI/dt
                 gamma * y(3)];                             % dR/dt

% Initial vector of state variables
y0 = [S0; E0; I0; R0];

% Solve the ODE system
[t, Y] = ode45(ode, t, y0);

% Plot the results
figure;
plot(t, Y(:,1), 'LineWidth', 2, 'DisplayName', 'Susceptible');
hold on;
plot(t, Y(:,2), 'LineWidth', 2, 'DisplayName', 'Exposed');
plot(t, Y(:,3), 'LineWidth', 2, 'DisplayName', 'Infected');
plot(t, Y(:,4), 'LineWidth', 2, 'DisplayName', 'Recovered');
xlabel('Days');
ylabel('Population');
title('SEIR Model for COVID-19 Simulation in Kenya for five years from 2nd April 2020');
legend('Location', 'best');
grid on;
hold off;
