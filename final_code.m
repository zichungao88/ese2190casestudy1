%% Task 2: Initial Visualizations of Solutions & Phase Portraits
% These plots correspond to Figures 2-4 in the paper.

% We will rewrite the VDP equation as a first order differential equation
% We can write the VDP equation in vector form as:
% d/dt [x, y] = [y, mu*(1 - x^2)y - x]
% This can be modelled by using ode45 on a 2-vector x where
% x(1) = x
% x(2) = y

% both functions defined below
vdp_solution([0.25 0.5 1 4], [0 2]); % plot x and dx/dt
vdp_field([0.25 0.5 1 4], [0 1; 0 2; 0 3;]) % plot system as vector field

vdp_solution(1000, [0 2]); % comment out line 30 for mu = 1000

% Stores the functions used in this program in order to make space more
% efficient
function vdp_solution(mu_array, x0)
    time_range = [0, 20]; % time range from 0 to 20 (longer --> better deduction of trends
    figure; % solution over time plots (4 total)
    for i = 1:length(mu_array)
        mu = mu_array(i); % select each value of mu
        [t, x] = ode45(@(t, x) [x(2); mu * (1 - x(1)^2) * x(2) - x(1)], time_range, x0); % solve
        
        % 4-in-1 plot for concision, neatness, & ease of visualization
        subplot(ceil(length(mu_array)/2), ceil(length(mu_array)/2), i);
        hold on;
        plot(t, x(:,1), 'b', 'LineWidth', 1.5);
        plot(t, x(:,2), 'r', 'LineWidth', 1.5); % comment out this line when plotting mu = 1000
        xlabel('Time t'); ylabel('Solution');
        title(['Van der Pol: \mu = ', num2str(mu)]);
        legend('x', 'dx/dt');
        grid on;
        hold off;
    end
end

function vdp_field(mu_array, initial_conditions)
    % define vector field grid range
    x_min = -3; x_max = 3;
    y_min = -3; y_max = 3;
    [x, y] = meshgrid(linspace(x_min, x_max, 20), linspace(y_min, y_max, 20));
    sub_length = ceil(length(mu_array)/2);
    
    figure; % phase portrait plots (4 total)
    for i = 1:length(mu_array)
        mu = mu_array(i); % select each value of mu
        dx = y;
        dy = mu * (1 - x^2) * y - x;
        L = sqrt(dx.^2 + dy.^2);

        % normalize vectors for better visualization
        dx = dx ./ L;
        dy = dy ./ L;

        subplot(sub_length, sub_length, i);
        quiver(x, y, dx, dy, 'b', 'LineWidth', 1.2, 'DisplayName', "Vector Field"); hold on;
        xlabel('x'); ylabel('dx/dt');
        title(['Phase Portrait: \mu = ', num2str(mu)]);
        grid on;
        axis equal;
        
        time_range = [0, 20];
        
        % create limit cycles
        for j = 1:size(initial_conditions, 1)
            % solve system for given mu
            [~, z] = ode45(@(t, x) [x(2); mu * (1 - x(1)^2) * x(2) - x(1)], ...
                time_range, initial_conditions(j, :));
            % plot limit cycle
            plot(z(: ,1), z(: ,2), 'LineWidth', 1.5, 'DisplayName', ['IC: [', ...
                num2str(initial_conditions(j, 1)), ', ', num2str(initial_conditions(j, 2)), ']']);
        end
        legend;
    end
end


%% Task 3: Vector Field & Trajectory for mu = 0 (Simple Harmonic motion)
% These plots correspond to Figure 6 in the paper.

vdp = @(t, y, mu) [y(2); mu * (1 - y(1)^2) * y(2) - y(1)];

time_range = [0, 20];
y0 = [1; 0]; % initial conditions: x(0) = 1, dx/dt(0) = 0
mu = 0; % only value of mu we wish to plot is 0
[x, y] = meshgrid(linspace(-3, 3, 20), linspace(-3, 3, 20)); % create the vector field

[~, z] = ode45(@(t, y) vdp(t, y, mu), time_range, y0); % solve

% initializing the vector field in x & y directions
dxdt = y;
dydt = mu * (1 - x^2) * y - x;

% normalizing the vectors to better fit on the grid
dxdt = dxdt ./ sqrt(dxdt.^2 + dydt.^2);
dydt = dydt ./ sqrt(dxdt.^2 + dydt.^2);

% plotting
figure;
hold on;
quiver(x, y, dxdt, dydt, 'r', 'DisplayName', 'Vector Field');
plot(z(:,1), z(:,2), 'b-', 'DisplayName', 'Trajectory');

title('Van der Pol Oscillator for \mu = 0');
xlabel('x');
ylabel('dx/dt');
legend('Location', 'Best');
grid on;
axis equal;
hold off;


%% Task 4: Limit Cycle Oscillations from mu = 0.01 to 100
% These plots correspond to Figures 7-16 in the paper.

vdp = @(t, y, mu) [y(2); mu * (1 - y(1)^2) * y(2) - y(1)];

time_range = [0 20]; % time range (exactly 1 full traversal)
y0 = [1; 0]; % initial condition: x(0) = 1, dx/dt(0) = 0
mu_range = [0, 0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 100]; % one plot for one mu
[x, y] = meshgrid(linspace(-3, 3, 20), linspace(-3, 3, 20)); % vector field

for i = 1:length(mu_range)
    mu = mu_range(i); % for each mu
    [t, z] = ode45(@(t, y) vdp(t, y, mu), time_range, y0); % solve
    
    % initialize vector field in x & y directions
    dxdt = y;
    dydt = mu * (1 - x^2) * y - x;

    % normalize vectors to better fit on grid
    dxdt = dxdt ./ sqrt(dxdt.^2 + dydt.^2);
    dydt = dydt ./ sqrt(dxdt.^2 + dydt.^2);

    % animation plot
    figure;
    hold on;
    quiver(x, y, dxdt, dydt, 'r', 'DisplayName', 'Vector Field'); % vector field lines
    plot(z(:, 1), z(:, 2), 'b-', 'DisplayName', 'Trajectory'); % full trajectory
    dot = plot(z(1, 1), z(1, 2), 'ro', 'MarkerFaceColor', 'g', ...
         'DisplayName', 'Current Position'); % moving point

    title(['Van der Pol Oscillator for \mu = ', num2str(mu)]);
    xlabel('x');
    ylabel('dx/dt');
    legend('Location', 'Best');
    grid on;
    axis equal;
    
    % animation loop
    for j = 1:length(t)
        set(dot, 'XData', z(j, 1), 'YData', z(j, 2)); % update moving point
        pause(0.1); % control animation speed
    end
    hold off;
end


%% Task 6: Verification of Calculated Limit Cycle Radius

mu = 0.01; % small mu

vanderpol = @(t, y) [y(2); mu * (1 - y(1)^2) * y(2) - y(1)];

time_range = [0 1000]; % r converges to theoretical value as time range increases
y0 = [1; 0];

[t, y] = ode45(vanderpol, time_range, y0);

r = sqrt(y(:, 1).^2 + y(:, 2).^2); % radius

% estimate limit cycle radius as stable peak value
% take max of r for stability purposes
radius_estimate = max(r);

disp(radius_estimate);