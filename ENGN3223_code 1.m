%% Aircraft Pitch Control using Transformation Matrix Method
% Based on provided longitudinal dynamics and state-space formulation

clc; clear; close all;

%% State-space matrices from the document (Eq. 14)
A = [-2.02     1       0;
     -6.986  -2.9476   0;
      0       1        0];

B = [0.16; 
     11.7304; 
     0];

C = eye(3);  % should be [0 0 1], but for the simulink to work, it needs to be doing the y isolation instead
D = zeros(1, 3); % again, should be [0], but size is adjusted for the simulink

%% Step 1: Check controllability
Co = ctrb(A,B);
rankCo = rank(Co);
n = size(A,1);

fprintf('Rank of controllability matrix: %d (system order = %d)\n', rankCo, n);

if rankCo == n
    disp('System is controllable.');
else
    disp('System is NOT controllable.');
end

%% Step 2: Characteristic polynomial of A
char_poly = poly(A);
disp('Characteristic polynomial of A:');
disp(char_poly);

% This gives coefficients of |sI - A| = s^n + a1*s^(n-1) + ... + an

%% Step 3: Desired pole placement
% Design requirements
target_pitch_deg = 0.2;   % (for reference, scaling the step input)
rise_time = 0.771;       % seconds
overshoot = 10;          % percent

% Compute damping ratio and natural frequency
zeta = sqrt((-1/pi*log(overshoot/100))^2/(1+(-1/pi*log(overshoot/100))^2));
wn   = 1.8 / rise_time;

% Dominant pole pair
s1 = -zeta*wn + 1i*wn*sqrt(1-zeta^2);
s2 = conj(s1);

% Place third pole further left (e.g., 3x farther than real part of s1)
s3 = 3*real(s1);

desired_poles = [s1 s2 s3];
disp('Desired poles:');
disp(desired_poles);

% Compute desired characteristic polynomial
desired_char_poly = poly(desired_poles);
disp('Desired characteristic polynomial:');
disp(desired_char_poly);

%% Step 4: Transformation matrix T (Eq. 16)
% Canonical transformation matrix for controllability form
% T = [B  AB  A^2B ... A^(n-1)B]
T = Co;

%% Step 5: Compute feedback gain K
% Ackermann's formula: K = [0 ... 0 1] * inv(T) * desired_char_poly(A)
% MATLAB has built-in 'acker' for this
K = acker(A,B,desired_poles);

disp('State feedback gain matrix K:');
disp(K);

%% Step 6: Closed-loop system
A_cl = A - B*K;

%% Step 7: Verify closed-loop poles
cl_poles = eig(A_cl);
disp('Closed-loop poles:');
disp(cl_poles);

% Step 8: Simulate the closed-loop system
out = sim('ENGN3223_system');

% Extract timeseries from SimulationOutput
y_ts = out.y_out;   % timeseries object

% Get time and data arrays
t = y_ts.Time;
y = y_ts.Data;

% Plot output (red) and setpoint (blue, zero line)
figure; clf; hold on;

pitch = plot(t, y, 'r', 'LineWidth', 1.5);       
set = plot(t, target_pitch_deg*ones(size(t)), 'b', 'LineWidth', 1.5);

legend([pitch(1); set(1)], {'Response System', 'Setpoint'});
xlabel('Time (s)');
ylabel('Pitch Angle (deg)');
grid on;
