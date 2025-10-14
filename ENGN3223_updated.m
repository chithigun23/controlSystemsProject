%% Aircraft Pitch Control using Transformation Matrix Method
% Based on provided longitudinal dynamics and state-space formulation

clc; clear; close all;

%% Variables

% Aircraft coefficients (Convair 880)
m = 126000*0.453592; % lb to kg
I_y = 2450000*1.35581795; % Slug*ft^2 to kg*m^2
C_D_u = 0.01;
C_L_u = -0.01;
C_D_0 = 0.02;
C_L_0 = 0.2;
C_D_al = 0.27;
C_L_al = 4.52;
C_m_al = -0.903;
C_m_aldot = -4.13;
C_m_q = -12.1;

% Dynamic pressure and other variables
rho = 1.225;      % kg/m^3 (air density)
u_0  = 230;        % m/s (trim speed)
S   = 2000*0.092903; % m^2 (wing area)
Q = 0.5*rho*u_0^2; % Pa (dynamic pressure)
cbar = 18.94*0.3048; % feet to metres (mean chord length)
g = 9.81; % m/s^2 (gravity)

% U-derivatives
X_u = -(C_D_u + 2*C_D_0)*Q*S/(u_0*m);
Z_u = -(C_L_u + 2*C_L_0)*Q*S/(u_0*m);
M_u = 0;

% W-derivatives
X_w = -(C_D_al - C_L_0)*Q*S/(u_0*m);
Z_w = -(C_L_al + C_D_0)*Q*S/(u_0*m);
M_w = C_m_al*Q*S*cbar/(u_0*I_y);

% Q-derivative
M_q = C_m_q*cbar/(2*u_0)*Q*S*cbar/I_y;

% Wdot-derivatives
X_wdot = 0;
Z_wdot = 0;
M_wdot = C_m_aldot*cbar/(2*u_0)*Q*S*cbar/(u_0*I_y);

% delta-derivatives
C_L_del  = 0.213;    
C_m_del  = -0.637;   
C_X_del  = -0.005;   % guessed/small (axial)
% Note: C_Z_delta (body axis) often = -C_L_delta depending on sign convention:
C_Z_del  = -C_L_del;

Z_del = Q * S * C_Z_del / m;       
X_del = Q * S * C_X_del / m;       
M_del = Q * S * cbar * C_m_del / I_y;   

% Thrust-derivatives
dT_dthrottle = 4e4;   % N per (throttle unit 0..1) -- set to your engine slope
l_T = 2.0;            % m moment arm about CG (positive = nose-up moment when thrust increases)
X_del_T = dT_dthrottle / m;     % N per throttle
Z_del_T = dT_dthrottle * sin(0) / m;  % approx 0 at small alpha
M_del_T = dT_dthrottle * l_T / I_y;     % N*m per throttle unit

%% State-space matrices from the document (Eq. 14)
A = [X_u              X_w            0              -g;
     Z_u              Z_w            u_0             0; 
     M_u+M_wdot*Z_u  M_w+M_wdot*Z_w  M_q+M_wdot*u_0  0
     0                0              1               0];

B = -[X_del               X_del_T; 
      Z_del               Z_del_T; 
      M_del+M_wdot*Z_del  M_del_T+M_wdot*Z_del_T;
      0                   0];

C = eye(4);  % should be [0 0 0 1], but for the simulink to work, it needs to be doing the y isolation instead
D = zeros(4, 2); % again, should be [0], but size is adjusted for the simulink

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
target_pitch_rad = 0.3;   % (for reference, scaling the step input)
rise_time = 2.113;       % seconds
overshoot = 5;     % percent

% Compute damping ratio and natural frequency
zeta = sqrt((-1/pi*log(overshoot/100))^2/(1+(-1/pi*log(overshoot/100))^2));
wn   = 1.8 / rise_time;

% Dominant pole pair
s1 = -zeta*wn + 1i*wn*sqrt(1-zeta^2);
s2 = conj(s1);

% Place remaining poles further left (faster dynamics)
s3 = 3 * real(s1);
s4 = 30 * real(s1);   % new fourth pole (further left than s3)

desired_poles = [s1 s2 s3 s4];
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
K = place(A,B,desired_poles);

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
set = plot(t, target_pitch_rad*ones(size(t)), 'b', 'LineWidth', 1.5);

legend([pitch(1); set(1)], {'Response System', 'Setpoint'});
xlabel('Time (s)');
ylabel('Pitch Angle (rad)');
grid on;