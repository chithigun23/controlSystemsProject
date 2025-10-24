% Reset command window, workspace and figures
clc; clear; close all;

% Aircraft coefficients (Convair 880)
m = 126000*0.453592; % mass (lb to kg)
I_y = 2450000*1.35581795; % pitch moment of inertia (slug*ft^2 to kg*m^2)
C_D_u = 0.01; % drag due to velocity
C_L_u = -0.01; % lift due to velocity
C_D_0 = 0.02; % reference drag coefficient
C_L_0 = 0.2; % reference lift coefficient
C_D_al = 0.27; % drag due to angle of attack
C_L_al = 4.52; % lift due to angle of attack
C_m_al = -0.903; % pitch moment due to angle of attack
C_m_aldot = -4.13; % pitch moment due to angle of attack rate
C_m_q = -12.1; % pitch moment due to pitch rate
C_L_del  = 0.213; % lift due to elevator deflection
C_m_del  = -0.637; % pitch moment due to elevator deflection
C_X_del  = -0.005; % axial force due to elevator deflection
C_Z_del  = -0.213; % normal force due to elevator deflection
dT_dthrottle = 4e4;   % Thrust derivative per throttle unit (N)
l_T = 2.0;            % Thrust moment arm (m)

% Dynamic pressure and other variables
rho = 1.225; % air density (kg/m^3)
u_0  = 230; % trim speed (m/s)
S   = 2000*0.092903; % wing area (m^2) 
Q = 0.5*rho*u_0^2; % dynamic pressure (Pa) 
cbar = 18.94*0.3048; % mean chord length (ft to m)
g = 9.81; % gravity (m/s^2)

% Calculate flight derivatives using equations in (12)
% U derivatives
X_u = -(C_D_u + 2*C_D_0)*Q*S/(u_0*m);
Z_u = -(C_L_u + 2*C_L_0)*Q*S/(u_0*m);
M_u = 0;

% W derivatives
X_w = -(C_D_al - C_L_0)*Q*S/(u_0*m);
Z_w = -(C_L_al + C_D_0)*Q*S/(u_0*m);
M_w = C_m_al*Q*S*cbar/(u_0*I_y);

% q derivative
M_q = C_m_q*cbar/(2*u_0)*Q*S*cbar/I_y;

% Wdot derivatives
X_wdot = 0;
Z_wdot = 0;
M_wdot = C_m_aldot*cbar/(2*u_0)*Q*S*cbar/(u_0*I_y);

% Delta derivatives
Z_del = Q * S * C_Z_del / m;       
X_del = Q * S * C_X_del / m;       
M_del = Q * S * cbar * C_m_del / I_y;   

% Thrust derivatives
X_del_T = dT_dthrottle / m; 
Z_del_T = dT_dthrottle * sin(0) / m;  % assuming small angle of attack
M_del_T = dT_dthrottle * l_T / I_y; 

% State-space matrices from the equations in (11)
A = [X_u              X_w            0              -g;
     Z_u              Z_w            u_0             0; 
     M_u+M_wdot*Z_u  M_w+M_wdot*Z_w  M_q+M_wdot*u_0  0
     0                0              1               0];

B = -[X_del               X_del_T; 
      Z_del               Z_del_T; 
      M_del+M_wdot*Z_del  M_del_T+M_wdot*Z_del_T;
      0                   0];

C = eye(4);  % Needs to be 4x4 for the simulink to work
D = zeros(4, 2); % Again, the size is adjusted for the simulink

% Check controllability and display results
W = ctrb(A,B);
rankW = rank(W); 
d_min = min(size(W, 1), size(W, 2)); % smallest dimension of W

if rankW == d_min
    disp('System is controllable.');
else
    disp('System is NOT controllable.');
end

disp('Controllability matrix W:');
disp(W);

% Determine pole placement based on design requirements
target_pitch_rad = 10*pi/180; % 10 degrees target (0.1745 rad)
rise_time = 11; % seconds
overshoot = 5; % percent
steady_state_error = 0; % percent

% Calculate the damping ratio and natural frequency using the equations in (13)
zeta = sqrt((-1/pi*log(overshoot/100))^2/(1+(-1/pi*log(overshoot/100))^2));
wn   = 1.8 / rise_time;

% Place the dominant pole pair using the equations in (14)
s1 = -zeta*wn + 1i*wn*sqrt(1-zeta^2);
s2 = conj(s1);

% Place the remaining poles further left on the real axis
s3 = 3 * real(s1);
s4 = 30 * real(s1); 

% Display the calculated poles
desired_poles = [s1 s2 s3 s4];
disp('Desired poles:');
disp(desired_poles);


% Compute gain matrix K using the desired poles using the place function
K = place(A,B,desired_poles);

% Display the resulting matrix K
disp('State feedback gain matrix K:');
disp(K);

% Closed-loop A matrix
A_cl = A - B*K;

% Calculate precompensator gain N_bar
B_r = -B(:, 1); % column of B corresponding to elevator input
C_r = C(4, :); % row of C corresponding to pitch angle output
D_r = D(4, 1); % corresponding row, column value from matrix D

% Compute the steady state gain with respect to the target setpoint
G0 = C_r * (A_cl \ B_r) + D_r;

% Compute N_bar for the given steady state error target
N_bar = (1-steady_state_error/100)/G0;

% Display resulting gain N_bar
disp('Computed N_bar:');
disp(N_bar);

% Step 9: Simulate the closed-loop system
out = sim('ENGN3223_final_system');

% Extract output from the simulink as a timeseries
y_ts = out.y_out; 

% Get the time and data arrays
t = y_ts.Time;
y = y_ts.Data;

% Plot output (red) and setpoint (blue)
figure; clf; hold on;

pitch = plot(t, y, 'r', 'LineWidth', 1.5);       
set = plot(t, target_pitch_rad*ones(size(t)), 'b', 'LineWidth', 1.5);

legend([pitch(1); set(1)], {'System output', 'Setpoint'});
xlabel('Time (s)');
ylabel('Pitch Angle (rad)');
grid on;

% Set y-axis range
ylim([0 0.25]);

% Display actual simulated characteristics
info = stepinfo(y, t, target_pitch_rad);
disp(info.RiseTime)
disp(info.Overshoot)
disp(info.SettlingTime)


