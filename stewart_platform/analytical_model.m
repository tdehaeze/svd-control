%% Clear Workspace and Close figures
clear; close all; clc;

%% Intialize Laplace variable
s = zpk('s');

%% Bode plot options
opts = bodeoptions('cstprefs');
opts.FreqUnits = 'Hz';
opts.MagUnits = 'abs';
opts.MagScale = 'log';
opts.PhaseWrapping = 'on';
opts.xlim = [1 1000];

% Characteristics

L  = 0.055; % Leg length [m]
Zc = 0;     % ?
m  = 0.2;   % Top platform mass [m]
k  = 1e3;   % Total vertical stiffness [N/m]
c  = 2*0.1*sqrt(k*m); % Damping ? [N/(m/s)]

Rx = 0.04; % ?
Rz = 0.04; % ?
Ix = m*Rx^2; % ?
Iy = m*Rx^2; % ?
Iz = m*Rz^2; % ?

% Mass Matrix

M = m*[1   0 0  0         Zc        0;
       0   1 0 -Zc        0         0;
       0   0 1  0         0         0;
       0 -Zc 0  Rx^2+Zc^2 0         0;
       Zc  0 0  0         Rx^2+Zc^2 0;
       0   0 0  0         0         Rz^2];

% Jacobian Matrix

Bj=1/sqrt(6)*[ 1             1          -2          1         1        -2;
               sqrt(3)      -sqrt(3)     0          sqrt(3)  -sqrt(3)   0;
               sqrt(2)       sqrt(2)     sqrt(2)    sqrt(2)   sqrt(2)   sqrt(2);
               0             0           L          L        -L         -L;
               -L*2/sqrt(3) -L*2/sqrt(3) L/sqrt(3)  L/sqrt(3) L/sqrt(3)  L/sqrt(3);
               L*sqrt(2)    -L*sqrt(2)   L*sqrt(2) -L*sqrt(2) L*sqrt(2) -L*sqrt(2)];

% Stifnness and Damping matrices

kv = k/3;     % Vertical Stiffness of the springs [N/m]
kh = 0.5*k/3; % Horizontal Stiffness of the springs [N/m]

K = diag([3*kh, 3*kh, 3*kv, 3*kv*Rx^2/2, 3*kv*Rx^2/2, 3*kh*Rx^2]); % Stiffness Matrix
C = c*K/100000; % Damping Matrix

% State Space System

A  = [ zeros(6) eye(6); ...
      -M\K     -M\C];
Bw = [zeros(6); -eye(6)];
Bu = [zeros(6); M\Bj];

Co = [-M\K -M\C];

D  = [zeros(6) M\Bj];

ST = ss(A,[Bw Bu],Co,D);



% - OUT 1-6: 6 dof
% - IN 1-6 : ground displacement in the directions of the legs
% - IN 7-12: forces in the actuators.

ST.StateName = {'x';'y';'z';'theta_x';'theta_y';'theta_z';...
                'dx';'dy';'dz';'dtheta_x';'dtheta_y';'dtheta_z'};

ST.InputName = {'w1';'w2';'w3';'w4';'w5';'w6';...
                'u1';'u2';'u3';'u4';'u5';'u6'};

ST.OutputName = {'ax';'ay';'az';'atheta_x';'atheta_y';'atheta_z'};

% Transmissibility

TR=ST*[eye(6); zeros(6)];

figure
subplot(231)
bodemag(TR(1,1));
subplot(232)
bodemag(TR(2,2));
subplot(233)
bodemag(TR(3,3));
subplot(234)
bodemag(TR(4,4));
subplot(235)
bodemag(TR(5,5));
subplot(236)
bodemag(TR(6,6));

% Real approximation of $G(j\omega)$ at decoupling frequency

sys1 = ST*[zeros(6); eye(6)]; % take only the forces inputs

dec_fr = 20;
H1 = evalfr(sys1,j*2*pi*dec_fr);
H2 = H1;
D = pinv(real(H2'*H2));
H1 = inv(D*real(H2'*diag(exp(j*angle(diag(H2*D*H2.'))/2)))) ;
[U,S,V] = svd(H1);

wf = logspace(-1,2,1000);
for i  = 1:length(wf)
    H = abs(evalfr(sys1,j*2*pi*wf(i)));
    H_dec = abs(evalfr(U'*sys1*V,j*2*pi*wf(i)));
    for j = 1:size(H,2)
        g_r1(i,j) =  (sum(H(j,:))-H(j,j))/H(j,j);
        g_r2(i,j) =  (sum(H_dec(j,:))-H_dec(j,j))/H_dec(j,j);
        %     keyboard
    end
    g_lim(i) = 0.5;
end

% Coupled and Decoupled Plant "Gershgorin Radii"

figure;
title('Coupled plant')
loglog(wf,g_r1(:,1),wf,g_r1(:,2),wf,g_r1(:,3),wf,g_r1(:,4),wf,g_r1(:,5),wf,g_r1(:,6),wf,g_lim,'--');
legend('$a_x$','$a_y$','$a_z$','$\theta_x$','$\theta_y$','$\theta_z$','Limit');
xlabel('Frequency (Hz)'); ylabel('Gershgorin Radii')



% #+name: fig:gershorin_raddii_coupled_analytical
% #+caption: Gershorin Raddi for the coupled plant
% #+RESULTS:
% [[file:figs/gershorin_raddii_coupled_analytical.png]]


figure;
title('Decoupled plant (10 Hz)')
loglog(wf,g_r2(:,1),wf,g_r2(:,2),wf,g_r2(:,3),wf,g_r2(:,4),wf,g_r2(:,5),wf,g_r2(:,6),wf,g_lim,'--');
legend('$S_1$','$S_2$','$S_3$','$S_4$','$S_5$','$S_6$','Limit');
xlabel('Frequency (Hz)'); ylabel('Gershgorin Radii')

% Decoupled Plant

figure;
bodemag(U'*sys1*V,opts)

% Controller

fc = 2*pi*0.1; % Crossover Frequency [rad/s]
c_gain = 50; %

cont = eye(6)*c_gain/(s+fc);

% Closed Loop System

FEEDIN  = [7:12]; % Input of controller
FEEDOUT = [1:6]; % Output of controller



% Centralized Control

STcen = feedback(ST, inv(Bj)*cont, FEEDIN, FEEDOUT);
TRcen = STcen*[eye(6); zeros(6)];



% SVD Control

STsvd = feedback(ST, pinv(V')*cont*pinv(U), FEEDIN, FEEDOUT);
TRsvd = STsvd*[eye(6); zeros(6)];

% Results

figure
subplot(231)
bodemag(TR(1,1),TRcen(1,1),TRsvd(1,1),opts)
legend('OL','Centralized','SVD')
subplot(232)
bodemag(TR(2,2),TRcen(2,2),TRsvd(2,2),opts)
legend('OL','Centralized','SVD')
subplot(233)
bodemag(TR(3,3),TRcen(3,3),TRsvd(3,3),opts)
legend('OL','Centralized','SVD')
subplot(234)
bodemag(TR(4,4),TRcen(4,4),TRsvd(4,4),opts)
legend('OL','Centralized','SVD')
subplot(235)
bodemag(TR(5,5),TRcen(5,5),TRsvd(5,5),opts)
legend('OL','Centralized','SVD')
subplot(236)
bodemag(TR(6,6),TRcen(6,6),TRsvd(6,6),opts)
legend('OL','Centralized','SVD')
