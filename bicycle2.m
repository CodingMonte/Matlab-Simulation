function dx = bicycle2(x,t,u)

% This function defines the bicycle model. 
%
% Input
%   t: current time [s].
%   x: column array storing state values.
%     x(1): vx, velocity of vehicle c.g. point in body-fixed x axis [m/s].
%     x(2): vy, velocity of vehicle c.g. point in body-fixed y axis [m/s].
%     x(3): gamma, yaw rate at vehicle c.g. point [rad/s].
%     x(4): X, displacement of vehicle c.g. point in inertial x axis [m].
%     x(5): Y, displacement of vehicle c.g. point in inertial y axis [m].
%     x(6): psi, yaw angle at vehicle c.g. point [rad].
%
% Input by global variable
%   m: total vehicle mass [kg]
%   Izz: moment of inertia of vehicle [kg*m*m]
%   lf: length b/w c.g. and front axle along body-fixed x axis [m]
%   lr: length b/w c.g. and rear axle along body-fixed x axis [m]
%   Cf: cornering stiffness of a single front tire [N/rad]
%   Cr: cornering stiffness of a single rear tire [N/rad]
%
% Output
%   dx: column array that stores time derivatives of states.

global m Iz Lf Lr Cf Cr mu

%  Compute slip angle of front and rear tires.
%  Here approximation of tan(a) = a is not used.
alp_f = atan2(x(2)+Lf*x(3),x(1))-u;
alp_r = atan2(x(2)-Lr*x(3),x(1));

fyf = - mu * Cf * alp_f;
fyr = - mu * Cr * alp_r;
Fyf = 2*fyf;
Fyr = 2*fyr;

%===============================================================================
%  Compute time derivative of states.
%  1. dvx/dt    = 0
%  2. dvy/dt    = (F_yf + F_yr)/m - vx*gamma
%  3. dgamma/dt = (lf*F_yf - lr*F_yr)/Izz
%  4. dX/dt     = vx*cos(psi) - vy*sin(psi)
%  5. dY/dt     = vx*sin(psi) + vy*cos(psi)
%  6. dpsi/dt   = gamma

dx(1) = 0; 
dx(2) = (Fyf + Fyr)/m - x(1)*x(3);
dx(3) = (Lf*Fyf - Lr*Fyr)/Iz;
dx(4) = x(1)*cos(x(6)) - x(2)*sin(x(6));
dx(5) = x(1)*sin(x(6)) + x(2)*cos(x(6));
dx(6) = x(3);
