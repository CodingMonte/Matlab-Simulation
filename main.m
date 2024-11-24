%clc;
close all;
clear;

% Generation of Desired Path
% time
t0=0;
ts=0.01;
tf=50;
t=t0:ts:tf;
iter=fix((tf-t0)/ts); %1501

syms a b c d e f

o = 3.7;
t1 = 12;
T = t1 + o;

T2 = 24;
% condition
eq0 = a*t1^5 + b*t1^4 + c*t1^3 + d*t1^2 + e*t1 + f == 0;
eq1 = 5*a*t1^4 + 4*b*t1^3 + 3*c*t1^2 + 2*d*t1 + e == 0; % y'(T) = 0
eq2 = 20*a*t1^3 + 12*b*t1^2 + 6*c*t1 + 2*d == 0; % y''(T) = 0
eq3 = a*T^5 + b*T^4 + c*T^3 + d*T^2 + e*T + f == -3.5; % y(T) = -3.75
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

y_ref(t <= t1) = 0;
y_ref(t1 < t & t < T) = A*t(t1 < t & t < T).^5 + B*t(t1 < t & t < T).^4 + C*t(t1 < t & t < T).^3 + D*t(t1 < t & t < T).^2 + E*t(t1 < t & t < T) + F;

y_ref(T <= t & t <= T2) = -3.5;
y_ref(T2 < t & t <= T2 + o) = -(A*(t(T2 < t & t <= T2 + o)-T2+t1).^5 + B*(t(T2 < t & t <=T2 + o)-T2+t1).^4 ...
    + C*(t(T2 < t & t <=T2 + o)-T2+t1).^3 + D*(t(T2 < t & t <=T2 + o)-T2+t1).^2 + E*(t(T2 < t & t <=T2 + o)-T2+t1) + F)-3.5;
y_ref(T2 + o < t) = 0;

% y_ref = 5*sin(0.5*t);

% Derivative of Desired Path
y_ref_dot = gradient(y_ref, t); % y의 시간에 대한 미분 (첫 번째 미분)

% Derivative of Desired Path
y_ref_ddot = gradient(y_ref_dot, t); % y의 시간에 대한 미분 (첫 번째 미분)

% Desired Heading angle
psi_ref = atan2(y_ref_dot,100*1000/3600); % y의 기울기에 대한 아크탄젠트를 통해 헤딩 각도 계산

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
mu = 0.85;          %            Road friction coefficient 

% Linear tire parameters(Global)
Cf = 38800;        % [N/rad]    Front tire cornering stiffness
Cr = 36530;        % [N/rad]    Rear tire cornering stiffness

% Input saturation

u_bar = 5*pi/180;
u     = 0;
sat   = 0;

% Initial condition

x         = [100*1000/3600 0 0 0 y_ref(1) + 0.4 psi_ref(1)];
z_hat     = [0.4 0];
z1_hat_dot = [0 0];

%     x(1): vx, velocity of vehicle c.g. point in body-fixed x axis [m/s].
%     x(2): vy, velocity of vehicle c.g. point in body-fixed y axis [m/s].
%     x(3): gamma, yaw rate at vehicle c.g. point [rad/s].
%     x(4): X, displacement of vehicle c.g. point in inertial x axis [m].
%     x(5): Y, displacement of vehicle c.g. point in inertial y axis [m].
%     x(6): psi, yaw angle at  vehicle c.g. point [rad].

% Reference of z
zd = zeros(1,iter+1);

% Control gain
k1 = 1.5;
k2 = 60;

% Observer gain
l1 = 20;     
l2 = 30;
p  = 0;     % Output feedback intial condition

% DSC gain
tau     = 0.01;
alpha_f = 0;        % Filltered virtual controller

% Weight function of RBFNN
l       = 5;
sigma   = 1;
Gamma1  = 5;
Gamma2  = 5;
N       = 2;
omega   = zeros(1,l);
W1_hat  = zeros(1,l);
W2_hat  = zeros(1,l);
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
L       =  1;
a1      =  (a23 + L*a43)/L;
a2      =  (a22 + L*a42);
a3      = -(a23+L*a43)/L;
a4      =  (a24 + L*a44 - L*a22 - L^2*a42);
g       =  b1 + L*b2;
h       =  a24 + L*a44 - x(1);
mu0 = 0.85;
g0      = mu0*g;
% d       = (1-mu)*g*u_bar;
d       = 0.7*g*u_bar;
% Initial condition of signals
varphi  = 0;
beta    = 0;
rho     = 3;
lambda  = 1;
rho_inf = 0.05;%%c
gamma1  = 10;
gamma2  = 2;
gamma3  = 5;
rho_dot = 0;
phase1 = 0.25;
phase2 = 0.17;
phase3 = 0.15;

for n = 1:iter
    
    if n*ts == 13
        mu = phase1;
    elseif n*ts == 22
        mu = phase2;
    elseif n*ts == 26
        mu = phase3;
    end

    err(n,1) = x(n,5) - y_ref(n);               % position error
    err(n,2) = x(n,2) + x(n,1)*(x(n,6)-psi_ref(n));
    err(n,3) = x(n,6) - psi_ref(n);             % heading angle error
    err(n,4) = x(n,3) - psi_ref_dot(n);

    z(n,1) = err(n,1) + L*err(n,3);             % modified error
    z(n,2) = err(n,2) + L*err(n,4);

    xi(n)           = (z_hat(n,1) - zd(n) - varphi(n))/rho(n);
    e1(n)           = log((1 + xi(n))/(1 - xi(n))) - zeta(n);
    z1_tilda(n)     = z(n,1) - z_hat(n,1);
    % alpha(n)        = -l1*z1_tilda(n) - (1-xi(n)^2)/2*k1*rho(n)*log((1 + xi(n))/(1 - xi(n)));
    k1(n) = rho(n)/2 + 1;
    alpha(n)        = -l1*z1_tilda(n) - (1-xi(n)^2)/2*k1(n)*log((1 + xi(n))/(1 - xi(n)));
    z_hat(n,2)      = p(n) + l2*z1_tilda(n);
    e2(n)           = z_hat(n,2) - alpha_f(n) - beta(n);
    rho_dot(n)      = -lambda*(rho(n) - rho_inf) + gamma1*(abs(varphi(n)));

    % Layer Matrix
    G1(n,:)   = [z(n,:) psi_ref_dot(n)];
    G_hat(n,:)= [z_hat(n,:) psi_ref_dot(n)];
    G2(n,:)   = [varphi(n) rho_dot(n)*xi(n) e2(n) alpha(n) alpha_f(n)]/(1*xi(n)^2)/rho(n);

    % Gaussian function
    for i = 1:l
        S1(n,i)     = exp(-(G1(n,:)*G1(n,:)')/N^2);
        S1_hat(n,i) = exp(-(G_hat(n,:)*G_hat(n,:)')/N^2);
        S2(n,i)     = exp(-(G2(n,:)*G2(n,:)')/N^2);
    end

    if n >= 2
        S1_dot(n,:) = (S1_hat(n,:) - S1_hat(n-1,:))/ts;
    else
        S1_dot(n,:) = zeros(1,l);
    end

    % Weight vector
    W1_hat       = Gamma1*(S1_hat(n,:) * z1_tilda(n) - omega(n,:));
    W2_hat_dot   = Gamma2*(S2(n,:)*e1(n) - sigma*W2_hat(n));

    NW1(n) = norm(W1_hat);
    NW2(n) = norm(W2_hat(n));

    % Signal
    if n >= 2
        omega_dot = -l1*z1_tilda(n)*S1_hat(n,:) + S1_dot(n,:)*z1_tilda(n) + sigma*W1_hat;
    else 
        omega_dot = zeros(1,l);
    end
    
    alpha_f_dot(n)  = (alpha(n) - alpha_f(n))/tau;    

    % Control input
    u(n) = (-k2*e2(n) - W1_hat * S1_hat(n,:)' - z1_tilda(n) + alpha_f_dot(n) - gamma3*beta(n) + z_hat(n,2)*d^2/2)/g0;

    if u(n) >= u_bar
        sat(n) = u_bar;

    elseif u(n) <= -u_bar
        sat(n) = -u_bar;

    else
        sat(n) = u(n);
    end

    % derivative of signal
    varphi_dot(n)   = -gamma2*varphi(n) + beta(n);
    beta_dot(n)     = -gamma3*beta(n) + g0*(sat(n) - u(n));
    zeta_dot(n)     = -k1(n)*zeta(n)/rho(n) + W2_hat(n,:) * S2(n,:)';
     
    
    % runge-kutta for signal
    alpha_f(n+1) = rk4('dsc',alpha_f(n),alpha_f_dot(n),ts,n);
    rho(n+1)     = rk4('dsc',rho(n),rho_dot(n),ts,n);
    varphi(n+1)  = rk4('dsc',varphi(n),varphi_dot(n),ts,n);
    zeta(n+1)    = rk4('dsc',zeta(n),zeta_dot(n),ts,n);
    beta(n+1)    = rk4('dsc',beta(n),beta_dot(n),ts,n);
    omega(n+1,:) = rk4('RBFNN',omega(n,:),omega_dot,ts,n);

    % observer
    z1_hat_dot(n) = l1*z1_tilda(n) + z_hat(n,2);
    p_dot(n)  = W1_hat * S1_hat(n,:)' + g0 * sat(n) + (l1*l2 + 1)*z1_tilda(n) - z_hat(n,2)*d^2/2 - l2*e2(n);
    % p_dot(n)  = W1_hat * S1_hat(n,:)' + (l1*l2 + 1)*z1_tilda(n);

    % runge-kutta
    p(n+1)         = rk4('dsc',p(n),p_dot(n),ts,n);
    z_hat(n+1,1)   = rk4('observer',z_hat(n,1),z1_hat_dot(n),ts,n);
    W2_hat(n+1,:)  = rk4('RBFNN',W2_hat(n,:),W2_hat_dot,ts,n);
    x(n+1,:)       = rk4('vehicle_full_dynamics',x(n,:),sat(n),ts,n);   

end

%-----------------------plot---------------------------%

figure('Position', [100, 100, 1200, 400]);

plot(x(1:iter,4) , x(1:iter,5),'r','LineWidth', 1);
hold on
plot(x(1:iter,4), y_ref(1:iter),'k:','LineWidth', 1.5);
hold on
plot([x(14/ts,4) x(14/ts,4)] , [-4 -3.6],'k','LineWidth', 2);
hold on
plot([x(22/ts,4) x(22/ts,4)] , [-4 -3.6],'k','LineWidth', 2);
hold on
plot([x(26/ts,4) x(26/ts,4)] , [-4 -3.6],'k','LineWidth', 2);
hold off

xlabel('X(m)');
ylabel('Y(m)');
h_legend = legend('actual','reference','Location', 'southeast','fontsize',12);
legend('box','off');
% legend의 현재 위치 가져오기
legend_pos = get(h_legend, 'Position');

% Y축 위치를 살짝 위로 이동 (Position: [left bottom width height])
legend_pos(2) = legend_pos(2) + 0.05; % bottom 값을 증가시켜 위로 이동

% 변경된 위치 적용
set(h_legend, 'Position', legend_pos);
axis([0 1100 -4 0.5])
text(x(6/ts,4)-5,-3.75,sprintf('\\mu = %.2f', 0.85),'fontsize',13)
text(x(16.5/ts,4)-5,-3.75,sprintf('\\mu = %.2f', phase1),'fontsize',13)
text(x(22.7/ts,4)-5,-3.75,sprintf('\\mu = %.2f', phase2),'fontsize',13)
text(x(31/ts,4)-5,-3.75,sprintf('\\mu = %.2f', phase3),'fontsize',13)

figure('Position', [100, 100, 1200, 400]);

plot(x(1:iter,4), x(1:iter,5),'g','LineWidth', 2);
plot(t(1:iter) , x(1:iter,6)*180/pi,'b','LineWidth', 1);
hold on
plot(t(1:iter), psi_ref(1:iter)*180/pi,'black--','LineWidth', 1);
hold off

xlabel('angle(m)');
ylabel('time(s)');
legend('\psi','\psi_{Ref}');
title('Tracking of Yaw angle');
% 
figure

plot(t(1:10:iter),rho(1:10:iter),'b');
hold on;
plot(t(1:10:iter),-rho(1:10:iter),'r');
hold on;
plot(t(1:10:iter),z_hat(1:10:iter,1),'k');
hold off;
xlabel('time(sec)');
ylabel('tracking error');
legend('$\rho$','-$\rho$', '$z_{1}$','Interpreter','latex','fontsize',14);
legend('box','off');
% axes('position',[0.5 0.16 0.4 0.3])
% box on
% plot(t(1:10:iter),rho(1:10:iter),'b')
% hold on
% plot(t(1:10:iter),-rho(1:10:iter),'r')
% hold on
% plot(t(1:10:iter),z_hat(1:10:iter,1),'k');
% hold off
% axis([11 21 -0.5 0.5])

figure('Position', [100, 100, 1200, 400]);

plot(t(1:iter),u(1:iter)*180/pi,'k:','LineWidth',1.5);
hold on;
plot(t(1:iter),sat(1:iter)*180/pi,'r','LineWidth',1);
hold off;

xlabel('time(sec)');
ylabel('steering angle(\circ)');
legend('\delta','\Phi(\delta)','Location', 'southeast','fontsize',14);
legend('box','off')
axes('position',[0.2 0.25 0.4 0.3])
box on
plot(t(1400:2400),u(1400:2400)*180/pi,'k:','LineWidth',1.5)
hold on
plot(t(1400:2400),sat(1400:2400)*180/pi,'r')
hold off
axis([14 24 -10 10])
plot(t(2300:3400),u(2300:3400)*180/pi,'k:','LineWidth',1.5)
hold on
plot(t(2300:3400),sat(2300:3400)*180/pi,'r')
hold off
axis([23 34 -10 10])

figure

plot(t(1:10:iter),z1_tilda(1:10:iter),'b',t(1:10:iter),z(1:10:iter,2)-z_hat(1:10:iter,2),'k'); 
xlabel('time(sec)');
ylabel('estimation');
legend('$\tilde{z}_{1}$','$\tilde{z}_{2}$', 'Interpreter', 'latex','fontsize',14);
legend('box','off');

figure

plot(t(1:10:iter),NW1(1:10:iter),'b',t(1:10:iter),NW2(1:10:iter),'r:'); 
xlabel('time(sec)');
ylabel('norm of weights of RBFNN');
legend('$\hat{W}_{1}$','$\hat{W}_{2}$', 'Interpreter', 'latex','fontsize',14);
legend('box','off');

figure

plot(t(1:10:iter),err(1:10:iter,1),'k',t(1:10:iter),err(1:10:iter,3),'k:'); 
xlabel('time(sec)');
ylabel('tracking error');
legend('$y$ - $y_{d}$','$\psi$ - $\psi_{d}$','Location', 'northeast', 'Interpreter', 'latex','fontsize',14);
legend('box','off');
% % 
% %save all figures
% hf = get(0,'children');
% hf = flipud(hf);
% for ix=1:length(hf)
%     set(hf(ix),'Units','centimeters');
%     pos = get(hf(ix),'Position');
%     set(hf(ix),'PaperPositionMode','Auto',...
%         'PaperUnits','centimeters','PaperSize',[pos(3), pos(4)]);
%     figurename = ['image_sin\' 'fig' num2str(ix)];
%     print(hf(ix), figurename, '-dpng', '-vector');
% end
