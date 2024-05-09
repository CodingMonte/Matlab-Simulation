clc;
close all;
clear all;

% Generation of Desired Path

% % time
% t0=0;
% ts=0.001;
% tf=16;
% t=t0:ts:tf;
% iter=fix((tf-t0)/ts); %16000
% 
% x_ref = 20*t; % x 좌표 설정
% y_ref = 5*cos(t) + 5*sin(t) + t; % y 좌표 설정
% 
% % x_ref = 30*sin(t) - 2*cos(3*t) + t; % x 좌표 설정
% % y_ref = 5*cos(t)+ 20*sin(2*t) - t; % y 좌표 설정

% time
t0=0;
ts=0.001;
tf=25;
t=t0:ts:tf;
iter=fix((tf-t0)/ts); %1501

syms a b c d e f

T = 10;
t1 = 4;

% condition
eq0 = a*t1^5 + b*t1^4 + c*t1^3 + d*t1^2 + e*t1 + f == 0;
eq1 = 5*a*t1^4 + 4*b*t1^3 + 3*c*t1^2 + 2*d*t1 + e == 0; % y'(T) = 0
eq2 = 20*a*t1^3 + 12*b*t1^2 + 6*c*t1 + 2*d == 0; % y''(T) = 0
eq3 = a*T^5 + b*T^4 + c*T^3 + d*T^2 + e*T + f == -2; % y(T) = -3.75
eq4 = 5*a*T^4 + 4*b*T^3 + 3*c*T^2 + 2*d*T + e == 0; % y'(T) = 0
eq5 = 20*a*T^3 + 12*b*T^2 + 6*c*T + 2*d == 0; % y''(T) = 0

% solution
solution = solve([eq0, eq1, eq2, eq3, eq4, eq5], [a, b, c, d, e, f]);

% parameter of path
A = double(solution.a);
B = double(solution.b);
C = double(solution.c);
D = double(solution.d);
E = double(solution.e);
F = double(solution.f);

% reference trajectory

y_ref = A*t.^5 + B*t.^4 + C*t.^3 + D*t.^2 + E*t + F;
y_ref(t <= t1) = 0;
y_ref(t >= T) = -2;

% Derivative of Desired Path
y_ref_dot = gradient(y_ref, t); % y의 시간에 대한 미분 (첫 번째 미분)

% Derivative of Desired Path
y_ref_ddot = gradient(y_ref_dot, t); % y의 시간에 대한 미분 (첫 번째 미분)

% Desired Heading angle
psi_ref = atan(y_ref_dot); % y의 기울기에 대한 아크탄젠트를 통해 헤딩 각도 계산

% Derivative of Heading angle
psi_ref_dot = gradient(psi_ref, t); % 헤딩 각도의 시간에 대한 미분 (첫 번째 미분)

% Second Derivative of Heading angle
psi_ref_ddot = gradient(psi_ref_dot, t); % 헤딩 각도의 시간에 대한 두 번째 미분

global m Iz Lf Lr Cf Cr mu

% Vehicle parameters(Global)
m  = 1515;          % [kg]      Total mass
Iz = 3392;          % [kgm^2]   Total inertia
Lf = 0.967;         % [m]       Distance between CoG and front tire
Lr = 1.673;         % [m]       Distance between Cog and rear tire
mu = 0.15;          %           Road friction coefficient 

% Linear tire parameters(Global)
Cf = 118800;        % [N/rad]   Front tire cornering stiffness
Cr = 165300;        % [N/rad]   Rear tire cornering stiffness

% Input saturation
u_bar = 5*pi/180;
u     = 0;
sat   = 0;

% Initial condition
x         = [120*1000/3600 0 0 0 (y_ref(1) - 0.1) psi_ref(1)];
z_hat     = [0 0];
z_hat_dot = [0 0];

% Reference of z
zd = zeros(1,iter+1);

% Control gain
k1 = 1;
k2 = 10;

% Observer gain
l1 = 15;     
l2 = 15;
p  = 0;     % Output feedback intial condition

% DSC gain
tau     = 0.05;
alpha_f = 0;        % Filltered virtual controller

% Weight function of RBFNN
l       = 5;
sigma   = 0.1;
Gamma   = 2;
N       = 2;
omega   = zeros(1,l);
W1_hat  = zeros(1,l);
W2_hat  = zeros(1,l);
W3_hat  = zeros(1,l);
zeta    = zeros(1,l);

% Error dynamics parameter
a22 = -2 * (mu*Cf + mu*Cr) / m / x(1);
a23 = -a22 * x(1);
a24 = -2 * (mu*Cf * Lf - mu*Cr * Lr)/ m / x(1) / x(1);
a42 = -2 * (Lf * mu * Cf - Lr * mu * Cr) / Iz;
a43 = -a42;
a44 = -2 * (Lf^2 * Cf + Lr^2 * Cr)/ Iz / x(1);
b1  =  2 * Cf / m / x(1);
b2  =  2 * Cf * Lf / Iz;

% Modified dynamics parameter
L       =  0.01;
a1      =  (a23 + L*a43)/L;
a2      =  (a22 + L*a42);
a3      = -(a23+L*a43)/L;
a4      =  (a24 + L*a44 - L*a22 - L^2*a42);
g       =  b1 + L*b2;
h       =  a24 + L*a44 - x(1);

g0      = mu*g;
d       = (1-mu)*g*u_bar;

% Initial condition of signals
varphi  = 0;
beta    = 0;
rho     = 4;
lambda  = 1;
rho_inf = 0.05;
gamma1  = 10;
gamma2  = 20;
gamma3  = 10;

%     x(1): vx, velocity of vehicle c.g. point in body-fixed x axis [m/s].
%     x(2): vy, velocity of vehicle c.g. point in body-fixed y axis [m/s].
%     x(3): gamma, yaw rate at vehicle c.g. point [rad/s].
%     x(4): X, displacement of vehicle c.g. point in inertial x axis [m].
%     x(5): Y, displacement of vehicle c.g. point in inertial y axis [m].
%     x(6): psi, yaw angle at  vehicle c.g. point [rad].

for n = 1:iter

    err(n,1) = x(n,5) - y_ref(n);           % position error
    err(n,2) = x(n,2) - y_ref_dot(n);
    err(n,3) = x(n,6) - psi_ref(n);         % heading angle error
    err(n,4) = x(n,3) - psi_ref_dot(n);

    z(n,1) = err(n,1) + L*err(n,3);         % modified error
    z(n,2) = err(n,2) + L*err(n,4);

    func = a1*z(n,1) + a2*z(n,2) + a3*err(n,1) + a4*err(n,4) + h*psi_ref_dot(n);

    xi(n)           = (z(n,1) - zd(n) - varphi(n))/rho(n);
    e1(n)           = log((1 + xi(n))/(1 - xi(n))) - zeta(n);
    z1_tilda(n)     = z(n,1) - z_hat(n,1);
    alpha(n)        = -l1*z1_tilda(n) - (1-xi(n)^2)/2*k1*log((1 + xi(n))/(1 - xi(n)));
    z_hat(n,2)      = p(n) + l2*z1_tilda(n);
    e2(n)           = z_hat(n,2) - alpha_f(n) - beta(n);


    % Layer Matrix
    G1(n,:)   = [z(n,:) psi_ref_dot(n)];
    G_hat(n,:)= [z_hat(n,:) psi_ref_dot(n)];
    G2(n,:)   = [varphi(n) rho(n) xi(n) e2(n) alpha(n) alpha_f(n)];
    G3(n,:)   = [z_hat(n,2) beta(n)];
    
    % Gaussian function
    for i = 1:l
        S1(n,i)     = exp(-(G1(n,:)*G1(n,:)')/N^2);
        S1_hat(n,i) = exp(-(G_hat(n,:)*G_hat(n,:)')/N^2);
        S2(n,i)     = exp(-(G2(n,:)*G2(n,:)')/N^2);
        S3(n,i)     = exp(-(G3(n,:)*G3(n,:)')/N^2);
    end
    
    if n >= 2
        S1_dot(n,:) = (S1_hat(n,:) - S1_hat(n-1,:))/ts;
    else
        S1_dot(n,:) = zeros(1,l);
    end

    % Weight vector
    W1_hat       = Gamma*(S1_hat(n,:) * z1_tilda(n) - omega(n,:));
    W2_hat_dot   = Gamma*(S2(n,:)*e1(n) - sigma*W2_hat(n));
    W3_hat_dot   = Gamma*(S3(n,:)*e2(n) - sigma*W3_hat(n));

    % Signal
    if n >= 2
        omega_dot = -l1*z1_tilda(n)*S1_hat(n,:) + S1_dot(n,:)*z1_tilda(n) + sigma*W1_hat;
    else
        omega_dot = zeros(1,l);
    end
    
    alpha_f_dot(n)  = (alpha(n) - alpha_f(n))/tau;    

    % Control input
    u(n) = (-k2*e2(n) -l2*e2(n)- W3_hat(n,:) * S3(n,:)' - W1_hat * S1_hat(n,:)' - z1_tilda(n) + z_hat(n,2)*d^2/2 + alpha_f_dot(n))/g0;

    if u(n) >= u_bar
        sat(n) = u_bar;

    elseif u(n) <= -u_bar
        sat(n) = -u_bar;

    else
        sat(n) = u(n);
    end
    % if u(n) >= 40*pi/180
    %     u(n) = 40*pi/180;
    % elseif u(n) <= -40*pi/180
    %     u(n) = -40*pi/180;
    % end

    % derivative of signal
    varphi_dot(n)   = -gamma2*varphi(n) + beta(n);
    zeta_dot(n)     = -k1/rho(n)*zeta(n) + W2_hat(n,:) * S2(n,:)';
    beta_dot(n)     = -gamma3*beta(n) + g0*(sat(n) - u(n));
    rho_dot(n)      = -lambda*(rho(n) - rho_inf) + gamma1*varphi(n);
    
    
    % runge-kutta for signal
    alpha_f(n+1) = rk4('dsc',alpha_f(n),alpha_f_dot(n),ts,n);
    rho(n+1)     = rk4('dsc',rho(n),rho_dot(n),ts,n);
    varphi(n+1)  = rk4('dsc',varphi(n),varphi_dot(n),ts,n);
    zeta(n+1)    = rk4('dsc',zeta(n),zeta_dot(n),ts,n);
    beta(n+1)    = rk4('dsc',beta(n),beta_dot(n),ts,n);
    omega(n+1,:) = rk4('RBFNN',omega(n,:),omega_dot,ts,n);

    % observer
    z_hat_dot(1) = l1*z1_tilda(n) + z_hat(n,2);
    p_dot(n)  = W1_hat * S1_hat(n,:)' + g0 * sat(n) + (l1*l2 + 1)*z1_tilda(n) +l2*e2(n) - z_hat(n,2)*d^2/2;
    
    % runge-kutta
    p(n+1)         = rk4('dsc',p(n),p_dot(n),ts,n);
    z_hat(n+1,1)   = rk4('observer',z_hat(n,1),z_hat_dot(1),ts,n);
    W2_hat(n+1,:)  = rk4('RBFNN',W2_hat(n,:),W2_hat_dot,ts,n);
    W3_hat(n+1,:)  = rk4('RBFNN',W3_hat(n,:),W3_hat_dot,ts,n);
    x(n+1,:)       = rk4('bicycle2',x(n,:),u(n),ts,n);   

end

%-----------------------plot---------------------------%

% figure;
% plot(t*20, y_ref, 'LineWidth', 2);         % trajectory ref
% title('Lane Change Reference Path');
% xlabel('x(m)');
% ylabel('y(m)');
% grid on;
% 
% figure;                             % yaw rate ref
% plot(t, psi_ref_deg, 'LineWidth', 1); 
% title('Yaw Rate of Reference');
% ylabel('angle($^\circ$)','Interpreter','latex');
% xlabel('time(sec)');

% grid on;
% 
% 
% figure;                             % curvature
% plot(t, kappa, 'LineWidth', 2);
% title('Curvature');
% xlabel('Time');
% ylabel('Curvature');
% grid on;
close all;

% 경로 및 각도 시각화
figure;
subplot(3, 1, 1);
% plot(x_ref, y_ref);
plot(t, y_ref);
xlabel('X 좌표');
ylabel('Y 좌표');
title('Reference Path');

subplot(3, 1, 2);
plot(t, psi_ref, 'r', t, psi_ref_dot, 'b');
xlabel('time(s)');
ylabel('angle $\degree$)');
legend('Heading angle', 'Derivative of Heading angle');
title('Heading Angle and its First Derivative');

subplot(3, 1, 3);
plot(t, psi_ref_ddot, 'g');
xlabel('time(s)');
ylabel('second derivative of heading angle');
title('Second Derivative of Heading Angle');

% figure
% 
% % plot(x(1:iter,4), x(1:iter,5),'g','LineWidth', 2);
% plot(x_ref(1:iter), x(1:iter,5),'b','LineWidth', 1);
% hold on
% plot(x_ref(1:iter), y_ref(1:iter),'r--','LineWidth', 1);
% hold off
% 
% xlabel('X(m)');
% ylabel('Y(m)');
% legend('State of Vehicle','Desired Path');
% title('State');
% 
% grid on;

figure

% plot(x(1:iter,4), x(1:iter,5),'g','LineWidth', 2);
plot(t(1:iter) , x(1:iter,5),'g','LineWidth', 1);
hold on
plot(t(1:iter), y_ref(1:iter),'r--','LineWidth', 1);
hold off

xlabel('X(m)');
ylabel('Y(m)');
legend('Y','Y_Ref');
title('Tracking of Lateral');

grid on;

figure

% plot(x(1:iter,4), x(1:iter,5),'g','LineWidth', 2);
plot(t(1:iter) , x(1:iter,6)*180/pi,'b','LineWidth', 1);
hold on
plot(t(1:iter), psi_ref(1:iter)*180/pi,'black--','LineWidth', 1);
hold off

xlabel('angle(m)');
ylabel('time(s)');
legend('psi','psi_Ref');
title('Tracking of Yaw angle');

grid on;

figure

plot(t(1:iter),rho(1:iter),t(1:iter),-rho(1:iter));
hold on;
plot(t(1:iter),z_hat(1:iter,1));
hold off;
xlabel('time(sec)');
ylabel('tracking error(m)');
legend('\rho','-\rho', 'Error');
title('Performance function');
grid on;

figure

plot(t(1:iter),u(1:iter)*180/pi,'g','LineWidth',2);
hold on;
plot(t(1:iter),sat(1:iter)*180/pi,'r--','LineWidth',2);
hold off;

xlabel('time(sec)');
ylabel('steering angle($\degree$)');
legend('u','\Phi');
title('Input and Input saturation');

figure

plot(t(1:iter),z_hat(1:iter,1),t(1:iter),z(1:iter,1)); 
xlabel('time(sec)');
ylabel('estimation');
legend('z\hat{}','z', 'Location', 'best', 'Interpreter', 'latex');
title('Estimation');

figure

plot(t(1:iter),err(1:iter,1),t(1:iter),err(1:iter,3)*180/pi); 
xlabel('time(sec)');
ylabel('estimation');
legend('Error of Lateral','Error of Heading angle', 'Location', 'best', 'Interpreter', 'latex');
title('Error');


