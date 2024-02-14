clc, clear all, close all
% Longitudinal states are: forward velocity (u), vertical velocity (w),
% pitch angle (theta), and pitch rate (q).
%
% Lateral states are: sideways velocity (v), roll rate(p), roll angle
% (phi), yaw rate (r).

%% Airplane Parameters
span = 3.067 ; %m
area = .6282 ; %m^2
chord = .208 ; %m
mass = 5.74 ; %kg
g = 9.81 ; %m/s^2
Ixx = 1.2009 ; %kg m^2
Iyy = .9318 ; %kg m^2
Izz = 2.0734 ; %kg m^2
Ixz = .0946 ; %kg m^2

%% Trim State
h = 1800 ; %m

%Longitudinal
Va = 18 ; %m/s
theta = .0515 ; %rad
u = 17.9762 ; %m/s
w = .9263 ; %m/s
q = 0 ;
delta_e = -.537 ; %rad
delta_t = 17.92 ; %percentage of max
du = 0 ;
dw = 0 ;
dq = 0 ;

%Lateral
delta_a = 0 ; %rad
delta_r = 0 ; %percentage of max
dv = 0 ;
dp = 0 ;
dr = 0 ;
v = 0 ; %m/s
p = 0 ; 
r = 0 ;
phi = 0 ;
B = 0 ;

% trim conditions for lon states: 
trim_lon_x = [17.9762 / 18; 0.9263; 0; 0.0515];
trim_lon_y = [17.9762 / 18; 0.0515];

%% Linear Model Parameters

%Force Parameters
Xu = -.1271 ; %s^-1
Xw = .6409 ; %s^-1
Xq = -.9106 ; %m/s
Yv = -.3714 ; %s^-1
Yp = .8254 ; %m/s
Yr = -17.6451 ; %m/s
Zu = -.07655 ; %s^-1
Zw = -6.3237 ; %s^-1
Zq = 16.9091 ; %m/s
X_delta_e = .0018 ; %m/s^2
X_delta_t = 3.3846 ; %m/s^2
Y_delta_a = -.0137 ; %m/s^2
Y_delta_r = .0556 ; %m/s^2
Z_delta_e = -.1234 ; %m/s^2

%Moment Parameters
Lv = -1.1467 ; %(m*s)^-1
Lp = -15.7093 ; %s^-1
Lr = 2.6774 ; %s^-1
Mu = .1090 ; %(m*s)^-1
Mw = -2.1148 ; %m*s)^-1
Mq = -3.2853 ; %s^-1
Nv = .6400 ; %(m*s)^-1
Np = -1.2356 ; %s^-1
Nr = -.5669 ; %s^-1
L_delta_a = -5.3580 ; %s^-2
L_delta_r = .0316 ; %s^-2
M_delta_e = -1.3996 ; %s^-2
N_delta_a = -.2566 ; %s^-2
N_delta_r = -.1309 ; %s^-2

%% State Space Matrices

A_lon = [Xu,    Xw/Va,  Xq/Va, -g*cos(theta)/Va; 
         Zu*Va,  Zw,      Zq,     -g*sin(theta); 
         Mu*Va,  Mw,      Mq,           0; 
         0,       0,       1,           0]  ;
A_lat = [Yv, Yp,   Yr,    g*cos(theta);
         Lv, Lp,   Lr,        0;
         Nv, Np,   Nr,        0;
         0,  1,   tan(theta), 0]  ;
B_lon = [X_delta_e/Va, X_delta_t/Va;
         Z_delta_e,         0      ;
         M_delta_e,         0      ;
            0,              0] ;
B_lat = [Y_delta_a, Y_delta_r;
         L_delta_a, L_delta_r;
         N_delta_a, N_delta_r;
             0,         0] ;
C_lon = [1, 0, 0, 0;
         0, 0, 0, 1] ;
C_lat = [1, 0, 0, 0;
         0, 0, 0, 1] ;
u_lon = [delta_e; delta_t] ;
u_lat = [delta_a; delta_r] ;
x_lon = [u/Va, w, q, theta]' ;
x_lat = [v, p, r, phi]' ;
ylon = [u/Va; theta] ;
ylat = [v; phi] ;
D = zeros(2,2) ;

states_lon = {'Forward Velocity' 'Vertical Velocity' 'Pitch Angular Velocity' 'Pitch Angle'};
inputs_lon = {'Elevator' 'Throttle'};
outputs_lon = {'Forward Velocity' 'Pitch'};
states_lat = {'Sideways Velocity' 'Roll Angular Velocity' 'Yaw Angular Velocity' 'Roll Angle'};
inputs_lat = {'Aileron' 'Rudder'};
outputs_lat = {'Sideways Velocity' 'Roll Angle'};
%% State Space Control

longitudinal = ss(A_lon, B_lon, C_lon, D, 'statename', states_lon, 'inputname', inputs_lon, 'outputname', outputs_lon) ;

poles_lon = eig(A_lon) ;

figure
step(longitudinal)
title('Uncompensated Longitudinal Step Response')
ylabel('Longitudinal State')
grid on

figure
bode(longitudinal)
title('Uncompensated Longitudinal Bode Response')
grid on

%% Proportional Control

% At this point we have a stable system for the longitudinal analysis.
% However, the poles are close enough to the origin that we'll move them anyway. 

%Controllability check
Co_lon = ctrb(A_lon,B_lon) ;
rank(Co_lon) ;
%Controllability matrix is full rank, and therefore is controllable. Now
%pick pole values st the lateral system dynamics are stable
poles_lon = eig(A_lon); k = 10;
corrected_poles_lon = [-1*k, -1*k, -.5*k, -.5*k] ;
[Kp_lon,prec,message] = place(A_lon,B_lon,corrected_poles_lon) %Kp_lon is proportional controller, prec is measurement of how closely eig(A-BKp) match specified pole location
A_lon_new = A_lon-B_lon*Kp_lon ;
sys_place = ss(A_lon_new,B_lon,C_lon,D,'statename', states_lon, 'inputname', inputs_lon, 'outputname', outputs_lon) ;

figure
step(sys_place)
grid on
title('Longitudinal Step Response')
ylabel('Longitudinal State')

figure
bode(sys_place)
title('Shifted Poles: Longitudinal Bode Response')
grid on

%% Feed Forward/Tracking Control
%We need to implement a control term s.t. the output can track the loop
%error from the input. This control should be s.t. u = kx+rF. The F term
%will allow y (output) to track r (input)


Z = [zeros(size(A_lon,1), size(B_lon,2)); eye(size(C_lon,1),size(C_lon,1))] ;
N = inv([A_lon, B_lon; C_lon, D])*Z ;
Nx = [N(1:size(A_lon,1),1:size(B_lon,2))] ;
Nu = [N((size(A_lon,1)+1):(size(A_lon,1)+size(C_lon,1)),1:size(B_lon,2))] ;
N_bar = Nu+Kp_lon*Nx ;


B_lon_new = B_lon*N_bar ;


sys_track = ss(A_lon_new,B_lon_new,C_lon,D,'statename', states_lon, 'inputname', inputs_lon, 'outputname', outputs_lon) ;
poles_lat_track = eig(A_lon-B_lon*Kp_lon) 

figure
step(sys_track)
grid on
title('Longitudinal Step Response with Reference Tracking')
ylabel('Longitudinal State')

figure
bode(sys_track)
title('Longitudinal Bode Response with Tracking')
grid on

%% Tuned Proportional Control Fast Response
fast_poles_lon = [-100 -100 -50 -50] ;
[Kp_lon_fast,prec,message] = place(A_lon,B_lon,fast_poles_lon) 

N_bar_fast = Nu+Kp_lon_fast*Nx ;

A_lon_fast = A_lon-B_lon*Kp_lon_fast ;
B_lon_fast = B_lon*N_bar_fast ;

sys_lon_fast = ss(A_lon_fast,B_lon_fast,C_lon,D,'statename', states_lon, 'inputname', inputs_lon, 'outputname', outputs_lon) ;
poles_lon_fast = eig(A_lon-B_lon*Kp_lon_fast)

figure
step(sys_lon_fast)
grid on
title('Lateral Step Response With Fast Poles')
ylabel('Lateral State')

figure
bode(sys_lon_fast)
title('Longitudinal Bode Response with Tracking')
grid on

%% Tuned Proportional Control Imaginary Poles
im_poles_lon = [-.5+j -.5-j -.950 -.950] ;
[Kp_lon_im,prec,message] = place(A_lon,B_lon,im_poles_lon) 

N_bar_im = Nu+Kp_lon_im*Nx ;
A_lon_im = A_lon-B_lat*Kp_lon_im ;
B_lon_im = B_lon*N_bar_im ;

sys_lon_im = ss(A_lon_im,B_lon_im,C_lon,D,'statename', states_lon, 'inputname', inputs_lon, 'outputname', outputs_lon) ;
poles_lon_im = eig(A_lon-B_lon*Kp_lon_im)

figure
step(sys_lon_im)
grid on
title('Lateral Step Response With Imaginary Poles')
ylabel('Lateral State')

figure
bode(sys_lon_im)
title('Longitudinal Bode Response with Tracking')
grid on

%% LQR Control - Changing Dynamics
% Prior control uses trial and error to tune gain. Now, let's use linear
% quadratic control to find optimal gain matrix to minimize \|\dot x \|^2.
% 
Q = A_lon_new'*A_lon_new
R = B_lon'*B_lon
N = 0.001*A_lon_new'*B_lon
[Kp_lqr, ARE, poles_lqr] = lqr(sys_place,Q,R,N)
% 
N_bar_lqr = Nu+Kp_lqr*Nx ;
A_lon_lqr = A_lon-B_lon*Kp_lqr ;
B_lon_lqr = B_lon*N_bar_lqr ;
% 
sys_lon_lqr = ss(A_lon_lqr,B_lon_lqr,C_lat,D,'statename', states_lon, 'inputname', inputs_lon, 'outputname', outputs_lon) ;
poles_lon_lqr = eig(A_lon_lqr)
% 
figure
step(sys_lon_lqr)
grid on
title('Longitudinal Step Response')
ylabel('Longitudinal State')

figure
bode(sys_lon_lqr)
title('Longitudinal Bode Response with Tracking using LQR')
grid on
%% LQR Tuner - choosing costs
%
Q = 200*eye(4) ; %The higher Q and lower R, better OS & rt
R = .02*eye(2) ;
N = 0 ;

[Kp_lqr_new, ARE_new, poles_lqr_new] = lqr(sys_place,Q,R,N)

N_bar_lqr_new = Nu+Kp_lqr_new*Nx ;
A_lon_lqr_new = A_lon-B_lon*Kp_lqr_new ;
B_lon_lqr_new = B_lon*N_bar_lqr_new ;

sys_lon_lqr_new = ss(A_lon_lqr_new,B_lon_lqr_new,C_lon,D,'statename', states_lon, 'inputname', inputs_lon, 'outputname', outputs_lon) ;
poles_lon_lqr_new = eig(A_lon_lqr_new)


figure
step(sys_lon_lqr_new)
grid on
title('Longitudinal Step Response')
ylabel('Longitudinal State')


figure
bode(sys_lon_lqr_new)
title('Longitudinal Bode Response with Tracking using LQR')
grid on

%% Plotting for checking trim conditions: 
xinit = x_dot_lon; t = linspace(0, 150, 10^3); tspan = [0, t(end)]; 
[t,x] = ode45(@sys_lon_control, tspan, xinit);
x_u = x(:, 1);
x_w = x(:, 2);
x_q = x(:, 3);
x_t = x(:, 4);  
plot_lon(t,x_u, x_t); 
%
figure;
step(longitudinal_new)
grid on
%% Functions
%
function dx = sys_lon_control(t, x)
% u and y are static, algebraic maps. So, we just need to integrate dx and
% then use its solution, x, to give the output map y. 
    global A_lon K_lon
    %
    A = A_lon; 
    B = [0.0001    0.1880;
        -0.1234         0;
        -1.3996         0;
         0         0]; 
    C = [1     0     0     0;
         0     0     0     1]; 
    D = [0     0;
         0     0];
    % K = K_lon;
    K = K_lon; % [K_lon(:,1), K_lon(:,2), zeros(2,1), zeros(2,1)];
    
    %
    G_lon = -C*(A-B*K)^-1*B; F= G_lon^-1; r = [0.9987; 0.0515];
    %
    u = -K*x + F*r;
    dx = A*x + B*u; 
end
%
function [x_u, x_w, x_q, x_t] = solver_RK4(xdot, t, yinit)
stepsize = t(2) - t(1);
m = length(yinit);
x = zeros(m, length(t) - 1);
x(:, 1) = yinit;
k1 = zeros(1, m); k2 = k1; k3 = k1; k4 = k1;
for ii=1:length(t)-1
    % calculate k values for T4
    k1 = xdot(t(ii),x(:, ii))';
    k2 = xdot(t(ii) + stepsize/2, x(:, ii) + (stepsize/2)*k1)';
    k3 = xdot(t(ii) + stepsize/2, x(:, ii) + (stepsize/2)*k2)';
    k4 = xdot(t(ii) + stepsize,x(:, ii) + stepsize*k3');
    T4 = (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    x(:, ii+1) = x(:, ii) + stepsize*T4;
end
%
x_u = x(1, :);
x_w = x(2, :);
x_q = x(3, :);
x_t = x(4, :);  
end
%
function plot_lon(t,x_u, x_t)
figure; 
subplot(2,1,1)
plot(t, x_u);
hold on
title('Forward Velocity Normalized to Air Speed')
plot(t, (17.9762/18)*ones(numel(t),1),'LineStyle','--')
legend({'Forward Velocity','Forward Velocity Reference at Trim Conditions'})
hold off
grid on
subplot(2,1,2)
plot(t, x_t);
hold on
title('Pitch')
plot(t, 0.0515.*ones(numel(t),1),'LineStyle','--')
legend({'Pitch','Pitch Reference at Trim Conditions'})
grid on
hold off
end