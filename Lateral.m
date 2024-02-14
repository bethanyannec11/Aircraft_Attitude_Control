clc, clear all, close all

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
A_lat = [Yv, Yp,   Yr,    g*cos(theta);
         Lv, Lp,   Lr,        0;
         Nv, Np,   Nr,        0;
         0,  1,   tan(theta), 0]  ;
B_lat = [Y_delta_a, Y_delta_r;
         L_delta_a, L_delta_r;
         N_delta_a, N_delta_r;
             0,         0] ;
C_lat = [1, 0, 0, 0;
         0, 0, 0, 1] ;
% C_lat = eye(4) ;
u_lat = [delta_a; delta_r] ;
x_lat = [v, p, r, phi]' ;
ylat = [v; phi] ;
D = zeros(2,2) ;

%% State Space Control

states_lat = {'Lateral Velocity' 'Roll Angular Velocity' 'Yaw Angular Velocity' 'Roll Angle'} ;
inputs_lat = {'Aileron' 'Rudder'} ;
outputs_lat = {'Lateral Velocity' 'Roll Angle'} ;
lateral = ss(A_lat, B_lat, C_lat, D, 'statename', states_lat, 'inputname', inputs_lat, 'outputname',outputs_lat) ;

[num1,den1] = ss2tf(A_lat, B_lat, C_lat, D,1) ;
[num2,den2] = ss2tf(A_lat, B_lat, C_lat, D,2) ;
yvel_aileron = tf(num1(1,:),den1) ;
yvel_rudder = tf(num2(1,:),den2) ;
roll_aileron = tf(num1(2,:),den1) ;
roll_rudder = tf(num2(2,:),den2) ;

poles_lat = eig(A_lat) ;

% One pole is in RHP, leading to instability in step response

figure 
step(lateral)
lateral_info = stepinfo(lateral) ;
title('Uncompensated Lateral Step Response')
ylabel('Lateral State')

figure 
bode(lateral)
title('Uncompensated Lateral Bode Response')

%% Proportional Control

%At this point we have a stable system for the longitudinal analysis.
%However, the lateral analysis contains one unstable pole. Use feedback
%control law u = kx in order to place poles, and change system stability.
%In order to do this, we need to check if the system itself is controllable

%Controllability check
Co_lat = ctrb(A_lat,B_lat) ;
rank(Co_lat) ;
%Controllability matrix is full rank, and therefore is controllable. Now
%pick pole values st the lateral system dynamics are stable
corrected_poles_lat = [-1, -1, -.5, -.5] ;
[Kp_lat,prec,message] = place(A_lat,B_lat,corrected_poles_lat) ; %Kp_lat is proportional controller, prec is measurement of how closely eig(A-BKp) match specified pole location
A_lat_new = A_lat-B_lat*Kp_lat ;

states_lat = {'Lateral Velocity' 'Roll Angular Velocity' 'Yaw Angular Velocity' 'Roll Angle'} ;
inputs_lat = {'Aileron' 'Rudder'} ;
outputs_lat = {'Lateral Velocity' 'Roll Angle'} ;
sys_place = ss(A_lat_new,B_lat,C_lat,D, 'statename', states_lat, 'inputname', inputs_lat, 'outputname',outputs_lat) ;

figure 
step(sys_place)
sys_place_info = stepinfo(sys_place) ;
title('Lateral Step Response with Pole Placement')
ylabel('Lateral State')

figure 
bode(sys_place)
title('Lateral Bode Response with Pole Placement')

%% Feed Forward/Tracking Control
% We need to implement a control term s.t. the output can track the loop
% error from the input. This control should be s.t. u = kx+rF. The F term
% will allow y (output) to track r (input)


Z = [zeros(size(A_lat,1),size(C_lat,1)); eye(2)] ;
N = inv([A_lat, B_lat; C_lat, D])*Z ;
Nx = [N(1:size(A_lat,1),1:size(C_lat,1))] ;
Nu = [N((size(A_lat,1)+1):(size(A_lat,1)+size(C_lat,1)),1:size(C_lat,1))] ;
N_bar = Nu+Kp_lat*Nx ;


B_lat_new = B_lat*N_bar ;

states_lat = {'Lateral Velocity' 'Roll Angular Velocity' 'Yaw Angular Velocity' 'Roll Angle'} ;
inputs_lat = {'Aileron' 'Rudder'} ;
outputs_lat = {'Lateral Velocity' 'Roll Angle'} ;
sys_track = ss(A_lat_new,B_lat_new,C_lat,D, 'statename', states_lat, 'inputname', inputs_lat, 'outputname',outputs_lat) ;
[num,den] = ss2tf(A_lat_new,B_lat_new,C_lat,D,1) ;


poles_lat_track = eig(A_lat-B_lat*Kp_lat)  ;

figure 
step(sys_track)
sys_track_info = stepinfo(sys_track) ;
title('Lateral Step Response with Reference Tracking')
ylabel('Lateral State')

figure 
bode(sys_place)
title('Lateral Bode Response with Reference Tracking')

%% Tuned Proportional Control Fast Response
fast_poles_lat = [-100 -100 -90 -50] ; %Poles are extremely fast. The larger the split between the poles, the greater the overshoot
[Kp_lat_fast,prec,message] = place(A_lat,B_lat,fast_poles_lat) ;

N_bar_fast = Nu+Kp_lat_fast*Nx ;

A_lat_fast = A_lat-B_lat*Kp_lat_fast ;
B_lat_fast = B_lat*N_bar_fast ;

states_lat = {'Lateral Velocity' 'Roll Angular Velocity' 'Yaw Angular Velocity' 'Roll Angle'} ;
inputs_lat = {'Aileron' 'Rudder'} ;
outputs_lat = {'Lateral Velocity' 'Roll Angle'} ;
sys_lat_fast = ss(A_lat_fast,B_lat,C_lat,D, 'statename', states_lat, 'inputname', inputs_lat, 'outputname',outputs_lat) ;
poles_lat_fast = eig(A_lat-B_lat*Kp_lat_fast) ;

figure
step(sys_lat_fast)
sys_fast_info = stepinfo(sys_lat_fast) ;
title('Lateral Step Response With Fast Poles')
ylabel('Lateral State')

figure 
bode(sys_place)
title('Lateral Bode Response with Fast Poles')

%% Tuned Proportional Control Imaginary Poles
im_poles_lat = [-1 -7+.7*j -7-.7*j -.5] ; % The higher the imaginary value, the worse the overshoot becomes on top right and bottom left
[Kp_lat_im,prec,message] = place(A_lat,B_lat,im_poles_lat) 

N_bar_im = Nu+Kp_lat_im*Nx ;
A_lat_im = A_lat-B_lat*Kp_lat_im ;
B_lat_im = B_lat*N_bar_im ;

states_lat = {'Lateral Velocity' 'Roll Angular Velocity' 'Yaw Angular Velocity' 'Roll Angle'} ;
inputs_lat = {'Aileron' 'Rudder'} ;
outputs_lat = {'Lateral Velocity' 'Roll Angle'} ;
sys_lat_im = ss(A_lat_im,B_lat_im,C_lat,D, 'statename', states_lat, 'inputname', inputs_lat, 'outputname',outputs_lat) ;
poles_lat_im = eig(A_lat-B_lat*Kp_lat_im) ;

figure
step(sys_lat_im)
sys_im_info = stepinfo(sys_lat_im) ;
title('Lateral Step Response With Imaginary Poles')
ylabel('Lateral State')

figure 
bode(sys_place)
title('Lateral Bode Response with Imaginary Poles')

%% LQR Control
%Prior control uses trial and error to tune gain. Now, let's use linear
%quadratic control to find optimal gain matrix.


Q = A_lat'*A_lat ;
R = B_lat'*B_lat ;
N = A_lat'*B_lat ;
[Kp_lqr, ARE, poles_lqr] = lqr(sys_place,Q,R,N) ;

N_bar_lqr = Nu+Kp_lqr*Nx ;
A_lat_lqr = A_lat-B_lat*Kp_lqr ;
B_lat_lqr = B_lat*N_bar_lqr ;

states_lat = {'Lateral Velocity' 'Roll Angular Velocity' 'Yaw Angular Velocity' 'Roll Angle'} ;
inputs_lat = {'Aileron' 'Rudder'} ;
outputs_lat = {'Lateral Velocity' 'Roll Angle'} ;
sys_lat_lqr = ss(A_lat_lqr,B_lat_lqr,C_lat,D, 'statename', states_lat, 'inputname', inputs_lat, 'outputname',outputs_lat) ;
poles_lat_lqr = eig(A_lat_lqr) ;

figure
step(sys_lat_lqr)
sys_lqr_info = stepinfo(sys_lat_lqr) ;
title('Lateral Step Response with LQR')
ylabel('Lateral State')

figure 
bode(sys_place)
title('Lateral Bode Response with LQR')

%% LQR Tuner
Q = 200*eye(4) ; %The higher Q and lower R, better OS & rt
R = .02*eye(2) ;
N = 0 ;

[Kp_lqr_new, ARE_new, poles_lqr_new] = lqr(sys_place,Q,R,N) ;

N_bar_lqr_new = Nu+Kp_lqr_new*Nx ;
A_lat_lqr_new = A_lat-B_lat*Kp_lqr_new ;
B_lat_lqr_new = B_lat*N_bar_lqr_new ;

states_lat = {'Lateral Velocity' 'Roll Angular Velocity' 'Yaw Angular Velocity' 'Roll Angle'} ;
inputs_lat = {'Aileron' 'Rudder'} ;
outputs_lat = {'Lateral Velocity' 'Roll Angle'} ;
sys_lat_lqr_new = ss(A_lat_lqr_new,B_lat_lqr_new,C_lat,D, 'statename', states_lat, 'inputname', inputs_lat, 'outputname',outputs_lat) ;
poles_lat_lqr_new = eig(A_lat_lqr_new) ;


figure
step(sys_lat_lqr_new)
sys_lqrnew_info = stepinfo(sys_lat_lqr_new) ;
title('Lateral Step Response with Tuned LQR')
ylabel('Lateral State')

figure 
bode(sys_place)
title('Lateral Bode Response with Tuned LQR')