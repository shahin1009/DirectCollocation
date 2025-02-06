clear; clc;

% Parameters
m_c = 1.0;          % Cart mass
m_p = 0.1;          % Pole mass
l = 0.5;            % Pole length
g = 9.81;           % gravity
track_limit = 2.0;  % cart position bound

x0=[1;pi;0;0];
xf=[0;0;0;0];

scale.pos = 1;      % position scale
scale.ang = 1;      % angle scale
scale.vel = 1;      % velocity scale
scale.force = 1;    % force scale


Q = diag([10, 10, 1, 1]);
R = 1;

params.x0=x0;
params.xf=xf;
params.scale=scale;
params.Q=Q;  
params.R=R;

N = 20;             % number of main nodes
tf = 2.0;           % final time
t = linspace(0, tf, N);
dt = tf/(N-1);

nStates = 4;        % number of states
nInputs = 1;        % number of inputs



% X = [states(N nodes); inputs(N nodes)]
nVars = N*nStates + N*nInputs;
x_init = zeros(nVars, 1);

% Set initial guess
for i = 1:N
    state_idx = (1:nStates) + (i-1)*nStates;
    t_norm = (i-1)/(N-1);  
    input_idx = nStates*N + i;
    x_init(input_idx) = 0.1;  % normalized force
end


%%

% Bounds (normalized)
lb = -inf(nVars, 1);
ub = inf(nVars, 1);

for i = 1:N
    state_idx = (1:nStates) + (i-1)*nStates;
    % Cart Position Bound
    lb(state_idx(1)) = -track_limit/scale.pos;
    ub(state_idx(1)) = track_limit/scale.pos;
    
    % Velocity Bound
    lb(state_idx(3:4)) = -10/scale.vel;
    ub(state_idx(3:4)) = 10/scale.vel;
end

for i = 1:N
    input_idx = nStates*N + i;
    lb(input_idx) = -10/scale.force;  
    ub(input_idx) = 10/scale.force;
end
%%
% Solve optimization problem
options = optimoptions('fmincon', ...
                       'Display', 'iter', ...
                       'MaxFunctionEvaluations', 1e5, ...
                       'MaxIterations', 1000, ...
                       'OptimalityTolerance', 1e-8, ...
                       'ConstraintTolerance', 1e-8);

[X_opt, fval] = fmincon(@(X) objective(X, N, nStates,params), x_init, [], [], [], [], lb, ub, ...
                        @(X) constraints(X, N, dt, m_c, m_p, l, g, nStates, params), options);

% Extract optimal trajectories
x_traj = zeros(N, nStates);
u_traj = zeros(N, 1);
for i = 1:N
    state_idx = (1:nStates) + (i-1)*nStates;
    input_idx = nStates*N + i;
    x_traj(i,:) = X_opt(state_idx)';
    u_traj(i) = X_opt(input_idx);
end
%%
% Visualization
figure('Position', [50 100 600 400]);

% Plot trajectories (denormalized)
subplot(311);
plot(t, x_traj(:,1)*scale.pos, 'b-', 'LineWidth', 2);
hold on;
plot(t, [-track_limit*ones(size(t)); track_limit*ones(size(t))], 'r--');
xlabel('Time (s)');
ylabel('Cart Position (m)');
title('Cart Position');
grid on;

subplot(312);
plot(t, mod(x_traj(:,2)*scale.ang, 2*pi), 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Pole Angle (rad)');
title('Pole Angle');
grid on;

subplot(313);
plot(t, u_traj*scale.force, 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Force (N)');
title('Control Input');
grid on;

% Animate cart-pole
animation_subplot = figure('Position', [100 100 1200 400]);
for i = 1:N
    cla(animation_subplot);
    
    % Cart
    cart_width = 0.4;
    cart_height = 0.2;
    cart_x = x_traj(i,1)*scale.pos;
    cart_y = 0;
    rectangle('Position', [cart_x-cart_width/2, cart_y-cart_height/2, cart_width, cart_height], ...
              'FaceColor', 'b');
    hold on;
    
    % Pole
    pole_x = [cart_x, cart_x + l*sin(x_traj(i,2)*scale.ang)];
    pole_y = [cart_y, cart_y + l*cos(x_traj(i,2)*scale.ang)];
    line(pole_x, pole_y, 'Color', 'r', 'LineWidth', 2);
    
    % Track limits
    line([-track_limit -track_limit], [-0.5 0.5], 'Color', 'k', 'LineStyle', '--');
    line([track_limit track_limit], [-0.5 0.5], 'Color', 'k', 'LineStyle', '--');
    
    axis equal;
    xlim([-track_limit-0.5 track_limit+0.5]);
    ylim([-0.5 1.5]);
    title(sprintf('Time: %.2f s', t(i)));
    grid on;
    drawnow;
    pause(0.05);
end




function [J] = objective(X, N, nStates,params)
    
    Q = params.Q;  
    R = params.R;                   
    
    J = 0;
    
    % State cost - minimize deviation from final upright position
    for i = 1:N
        state_idx = (1:nStates) + (i-1)*nStates;
        x_i = X(state_idx);

        
        % Calculate desired state at current time
        x_desired = zeros(4,1);      % Target upright position at origin

        % State error from desired position
        state_error = x_i - x_desired;

        if i~=N
            u_i = X(nStates*N + i);
            u_ip1 = X(nStates*N + i + 1);
            J = J + 0.5 * state_error' * Q * state_error+0.5*R*(u_i^2 + u_ip1^2);
        else
            J = J + 0.5 * state_error' * Q * state_error;
        end
    end
    
end



function [dx] = cartpole_dynamics(x_norm, u_norm, m_c, m_p, l, g, scale)
    x = [x_norm(1)*scale.pos;
         x_norm(2)*scale.ang;
         x_norm(3)*scale.vel;
         x_norm(4)*scale.vel];
    u = u_norm * scale.force;
    
    s = sin(x(2));
    c = cos(x(2));
    
    M = [m_c + m_p, m_p*l*c;
         m_p*l*c, m_p*l^2];
    C = [0; -m_p*l*x(4)^2*s];
    G = [0; m_p*g*l*s];
    B = [1; 0];
    
    acc = M \ (B*u - C - G);
    
    dx = [x(3)/scale.pos;
          x(4)/scale.ang;
          acc(1)/scale.vel;
          acc(2)/scale.vel];
end

% Constraint function with non-scaled time
function [c, ceq] = constraints(X, N, dt, m_c, m_p, l, g, nStates, params)
    ceq = [];
    
    scale=params.scale;
    x0=params.x0;
    xf=params.xf;

    % Collocation constraints
    for i = 1:N-1
        state_idx_i = (1:nStates) + (i-1)*nStates;
        state_idx_ip1 = (1:nStates) + i*nStates;
        x_i = X(state_idx_i);
        x_ip1 = X(state_idx_ip1);
        
        input_idx_i = nStates*N + i;
        input_idx_ip1 = nStates*N + i + 1;
        u_i = X(input_idx_i);
        u_ip1 = X(input_idx_ip1);
        
        dx_i = cartpole_dynamics(x_i, u_i, m_c, m_p, l, g, scale);
        dx_ip1 = cartpole_dynamics(x_ip1, u_ip1, m_c, m_p, l, g, scale);
        
        t_i = (i-1)*dt;     % Start time of interval
        t_ip1 = i*dt;       % End time of interval
        t_mid = t_i + dt/2; % Midpoint time
        
        % For each state, compute polynomial interpolation
        x_mid = zeros(nStates, 1);
        dx_mid = zeros(nStates, 1);
        
        for j = 1:nStates
            
            
            % For a cubic polynomial: a*t^3 + b*t^2 + c*t + d
            A = [t_i^3     t_i^2    t_i     1;      % x(t_i) = x_i
                 t_ip1^3   t_ip1^2  t_ip1   1;      % x(t_ip1) = x_ip1
                 3*t_i^2   2*t_i    1       0;      % dx(t_i) = dx_i
                 3*t_ip1^2 2*t_ip1  1       0];     % dx(t_ip1) = dx_ip1
            
            b = [x_i(j);
                 x_ip1(j);
                 dx_i(j);
                 dx_ip1(j)];
            
            % Solve for polynomial coefficients
            coeff = A\b;
            
            % Evaluate polynomial at collocation point
            x_mid(j) = coeff(1)*t_mid^3 + coeff(2)*t_mid^2 + coeff(3)*t_mid + coeff(4);
            
            % Evaluate derivative at collocation point
            dx_mid(j) = 3*coeff(1)*t_mid^2 + 2*coeff(2)*t_mid + coeff(3);
        end
        
        % Control input at collocation point
        u_mid = u_i + 0.5*(u_ip1 - u_i);
        
        % System dynamics at collocation point
        dx_sys = cartpole_dynamics(x_mid, u_mid, m_c, m_p, l, g, scale);
        
        
        ceq = [ceq; dx_mid - dx_sys];
    end
    
    % Boundary constraints
    x_initial = X(1:nStates);
    x_final = X((N-1)*nStates+1:N*nStates);
    
    ceq = [ceq;
           x_initial(1) - x0(1);   
           x_initial(2) - x0(2);   
           x_initial(3) - x0(3);   
           x_initial(4) - x0(4);   
           x_final(1) - xf(1);     
           x_final(2) - xf(2);     
           x_final(3) - xf(3);     
           x_final(4) - xf(4)];    
    
    c = []; 
end

