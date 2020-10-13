%% Clear Workspace and Close figures
clear; close all; clc;

%% Intialize Laplace variable
s = zpk('s');

addpath('STEP');

% Jacobian
% First, the position of the "joints" (points of force application) are estimated and the Jacobian computed.

open('drone_platform_jacobian.slx');

sim('drone_platform_jacobian');

Aa = [a1.Data(1,:);
      a2.Data(1,:);
      a3.Data(1,:);
      a4.Data(1,:);
      a5.Data(1,:);
      a6.Data(1,:)]';

Ab = [b1.Data(1,:);
      b2.Data(1,:);
      b3.Data(1,:);
      b4.Data(1,:);
      b5.Data(1,:);
      b6.Data(1,:)]';

As = (Ab - Aa)./vecnorm(Ab - Aa);

l = vecnorm(Ab - Aa)';

J = [As' , cross(Ab, As)'];

save('./jacobian.mat', 'Aa', 'Ab', 'As', 'l', 'J');

% Simscape Model

open('drone_platform.slx');



% Definition of spring parameters

kx = 0.5*1e3/3; % [N/m]
ky = 0.5*1e3/3;
kz = 1e3/3;

cx = 0.025; % [Nm/rad]
cy = 0.025;
cz = 0.025;



% We load the Jacobian.

load('./jacobian.mat', 'Aa', 'Ab', 'As', 'l', 'J');

% Identification of the plant
% The dynamics is identified from forces applied by each legs to the measured acceleration of the top platform.

%% Name of the Simulink File
mdl = 'drone_platform';

%% Input/Output definition
clear io; io_i = 1;
io(io_i) = linio([mdl, '/Dw'],              1, 'openinput');  io_i = io_i + 1;
io(io_i) = linio([mdl, '/u'],               1, 'openinput');  io_i = io_i + 1;
io(io_i) = linio([mdl, '/Inertial Sensor'], 1, 'openoutput'); io_i = io_i + 1;

G = linearize(mdl, io);
G.InputName  = {'Dwx', 'Dwy', 'Dwz', 'Rwx', 'Rwy', 'Rwz', ...
                'F1', 'F2', 'F3', 'F4', 'F5', 'F6'};
G.OutputName = {'Ax', 'Ay', 'Az', 'Arx', 'Ary', 'Arz'};



% There are 24 states (6dof for the bottom platform + 6dof for the top platform).

size(G)



% #+RESULTS:
% : State-space model with 6 outputs, 12 inputs, and 24 states.


% G = G*blkdiag(inv(J), eye(6));
% G.InputName  = {'Dw1', 'Dw2', 'Dw3', 'Dw4', 'Dw5', 'Dw6', ...
%                 'F1', 'F2', 'F3', 'F4', 'F5', 'F6'};



% Thanks to the Jacobian, we compute the transfer functions in the frame of the legs and in an inertial frame.

Gx = G*blkdiag(eye(6), inv(J'));
Gx.InputName  = {'Dwx', 'Dwy', 'Dwz', 'Rwx', 'Rwy', 'Rwz', ...
                 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'};

Gl = J*G;
Gl.OutputName  = {'A1', 'A2', 'A3', 'A4', 'A5', 'A6'};

% Obtained Dynamics

freqs = logspace(-1, 2, 1000);

figure;

ax1 = subplot(2, 1, 1);
hold on;
plot(freqs, abs(squeeze(freqresp(Gx('Ax', 'Fx'), freqs, 'Hz'))), 'DisplayName', '$A_x/F_x$');
plot(freqs, abs(squeeze(freqresp(Gx('Ay', 'Fy'), freqs, 'Hz'))), 'DisplayName', '$A_y/F_y$');
plot(freqs, abs(squeeze(freqresp(Gx('Az', 'Fz'), freqs, 'Hz'))), 'DisplayName', '$A_z/F_z$');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Amplitude [m/N]'); set(gca, 'XTickLabel',[]);
legend('location', 'southeast');

ax2 = subplot(2, 1, 2);
hold on;
plot(freqs, 180/pi*angle(squeeze(freqresp(Gx('Ax', 'Fx'), freqs, 'Hz'))));
plot(freqs, 180/pi*angle(squeeze(freqresp(Gx('Ay', 'Fy'), freqs, 'Hz'))));
plot(freqs, 180/pi*angle(squeeze(freqresp(Gx('Az', 'Fz'), freqs, 'Hz'))));
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'lin');
ylabel('Phase [deg]'); xlabel('Frequency [Hz]');
ylim([-180, 180]);
yticks([-360:90:360]);

linkaxes([ax1,ax2],'x');



% #+name: fig:stewart_platform_translations
% #+caption: Stewart Platform Plant from forces applied by the legs to the acceleration of the platform
% #+RESULTS:
% [[file:figs/stewart_platform_translations.png]]


freqs = logspace(-1, 2, 1000);

figure;

ax1 = subplot(2, 1, 1);
hold on;
plot(freqs, abs(squeeze(freqresp(Gx('Arx', 'Mx'), freqs, 'Hz'))), 'DisplayName', '$A_{R_x}/M_x$');
plot(freqs, abs(squeeze(freqresp(Gx('Ary', 'My'), freqs, 'Hz'))), 'DisplayName', '$A_{R_y}/M_y$');
plot(freqs, abs(squeeze(freqresp(Gx('Arz', 'Mz'), freqs, 'Hz'))), 'DisplayName', '$A_{R_z}/M_z$');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Amplitude [rad/(Nm)]'); set(gca, 'XTickLabel',[]);
legend('location', 'southeast');

ax2 = subplot(2, 1, 2);
hold on;
plot(freqs, 180/pi*angle(squeeze(freqresp(Gx('Arx', 'Mx'), freqs, 'Hz'))));
plot(freqs, 180/pi*angle(squeeze(freqresp(Gx('Ary', 'My'), freqs, 'Hz'))));
plot(freqs, 180/pi*angle(squeeze(freqresp(Gx('Arz', 'Mz'), freqs, 'Hz'))));
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'lin');
ylabel('Phase [deg]'); xlabel('Frequency [Hz]');
ylim([-180, 180]);
yticks([-360:90:360]);

linkaxes([ax1,ax2],'x');



% #+name: fig:stewart_platform_rotations
% #+caption: Stewart Platform Plant from torques applied by the legs to the angular acceleration of the platform
% #+RESULTS:
% [[file:figs/stewart_platform_rotations.png]]


freqs = logspace(-1, 2, 1000);

figure;

ax1 = subplot(2, 1, 1);
hold on;
for out_i = 1:5
  for in_i = i+1:6
    plot(freqs, abs(squeeze(freqresp(Gl(sprintf('A%i', out_i), sprintf('F%i', in_i)), freqs, 'Hz'))), 'color', [0, 0, 0, 0.2]);
  end
end
for ch_i = 1:6
  plot(freqs, abs(squeeze(freqresp(Gl(sprintf('A%i', ch_i), sprintf('F%i', ch_i)), freqs, 'Hz'))));
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Amplitude [m/N]'); set(gca, 'XTickLabel',[]);

ax2 = subplot(2, 1, 2);
hold on;
for ch_i = 1:6
  plot(freqs, 180/pi*angle(squeeze(freqresp(Gl(sprintf('A%i', ch_i), sprintf('F%i', ch_i)), freqs, 'Hz'))));
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'lin');
ylabel('Phase [deg]'); xlabel('Frequency [Hz]');
ylim([-180, 180]);
yticks([-360:90:360]);

linkaxes([ax1,ax2],'x');



% #+name: fig:stewart_platform_legs
% #+caption: Stewart Platform Plant from forces applied by the legs to displacement of the legs
% #+RESULTS:
% [[file:figs/stewart_platform_legs.png]]


freqs = logspace(-1, 2, 1000);

figure;

ax1 = subplot(2, 1, 1);
hold on;
% plot(freqs, abs(squeeze(freqresp(Gx('Ax', 'Dwx')/s^2, freqs, 'Hz'))), 'DisplayName', '$D_x/D_{w,x}$');
% plot(freqs, abs(squeeze(freqresp(Gx('Ay', 'Dwy')/s^2, freqs, 'Hz'))), 'DisplayName', '$D_y/D_{w,y}$');
% plot(freqs, abs(squeeze(freqresp(Gx('Az', 'Dwz')/s^2, freqs, 'Hz'))), 'DisplayName', '$D_z/D_{w,z}$');
set(gca,'ColorOrderIndex',1)
plot(freqs, abs(squeeze(freqresp(TR(1,1), freqs, 'Hz'))), '--', 'DisplayName', '$D_x/D_{w,x}$');
plot(freqs, abs(squeeze(freqresp(TR(2,2), freqs, 'Hz'))), '--', 'DisplayName', '$D_x/D_{w,x}$');
plot(freqs, abs(squeeze(freqresp(TR(3,3), freqs, 'Hz'))), '--', 'DisplayName', '$D_x/D_{w,x}$');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Transmissibility - Translations');  xlabel('Frequency [Hz]');
legend('location', 'northeast');

ax2 = subplot(2, 1, 2);
hold on;
% plot(freqs, abs(squeeze(freqresp(Gx('Arx', 'Rwx')/s^2, freqs, 'Hz'))), 'DisplayName', '$R_x/R_{w,x}$');
% plot(freqs, abs(squeeze(freqresp(Gx('Ary', 'Rwy')/s^2, freqs, 'Hz'))), 'DisplayName', '$R_y/R_{w,y}$');
% plot(freqs, abs(squeeze(freqresp(Gx('Arz', 'Rwz')/s^2, freqs, 'Hz'))), 'DisplayName', '$R_z/R_{w,z}$');
set(gca,'ColorOrderIndex',1)
plot(freqs, abs(squeeze(freqresp(TR(4,4), freqs, 'Hz'))), '--', 'DisplayName', '$D_x/D_{w,x}$');
plot(freqs, abs(squeeze(freqresp(TR(5,5), freqs, 'Hz'))), '--', 'DisplayName', '$D_x/D_{w,x}$');
plot(freqs, abs(squeeze(freqresp(TR(6,6), freqs, 'Hz'))), '--', 'DisplayName', '$D_x/D_{w,x}$');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Transmissibility - Rotations');  xlabel('Frequency [Hz]');
legend('location', 'northeast');

linkaxes([ax1,ax2],'x');

% Real Approximation of $G$ at the decoupling frequency
% Let's compute a real approximation of the complex matrix $H_1$ which corresponds to the the transfer function $G_c(j\omega_c)$ from forces applied by the actuators to the measured acceleration of the top platform evaluated at the frequency $\omega_c$.

wc = 2*pi*20; % Decoupling frequency [rad/s]

Gc = G({'Ax', 'Ay', 'Az', 'Arx', 'Ary', 'Arz'}, ...
       {'F1', 'F2', 'F3', 'F4', 'F5', 'F6'}); % Transfer function to find a real approximation

H1 = evalfr(Gc, j*wc);



% The real approximation is computed as follows:

D = pinv(real(H1'*H1));
H1 = inv(D*real(H1'*diag(exp(j*angle(diag(H1*D*H1.'))/2))));

% Verification of the decoupling using the "Gershgorin Radii"
% First, the Singular Value Decomposition of $H_1$ is performed:
% \[ H_1 = U \Sigma V^H \]


[U,S,V] = svd(H1);



% Then, the "Gershgorin Radii" is computed for the plant $G_c(s)$ and the "SVD Decoupled Plant" $G_d(s)$:
% \[ G_d(s) = U^T G_c(s) V \]

% This is computed over the following frequencies.

freqs = logspace(-2, 2, 1000); % [Hz]



% Gershgorin Radii for the coupled plant:

Gr_coupled = zeros(length(freqs), size(Gc,2));

H = abs(squeeze(freqresp(Gc, freqs, 'Hz')));
for out_i = 1:size(Gc,2)
    Gr_coupled(:, out_i) = squeeze((sum(H(out_i,:,:)) - H(out_i,out_i,:))./H(out_i, out_i, :));
end



% Gershgorin Radii for the decoupled plant using SVD:

Gd = U'*Gc*V;
Gr_decoupled = zeros(length(freqs), size(Gd,2));

H = abs(squeeze(freqresp(Gd, freqs, 'Hz')));
for out_i = 1:size(Gd,2)
    Gr_decoupled(:, out_i) = squeeze((sum(H(out_i,:,:)) - H(out_i,out_i,:))./H(out_i, out_i, :));
end



% Gershgorin Radii for the decoupled plant using the Jacobian:

Gj = Gc*inv(J');
Gr_jacobian = zeros(length(freqs), size(Gj,2));

H = abs(squeeze(freqresp(Gj, freqs, 'Hz')));

for out_i = 1:size(Gj,2)
    Gr_jacobian(:, out_i) = squeeze((sum(H(out_i,:,:)) - H(out_i,out_i,:))./H(out_i, out_i, :));
end

figure;
hold on;
plot(freqs, Gr_coupled(:,1), 'DisplayName', 'Coupled');
plot(freqs, Gr_decoupled(:,1), 'DisplayName', 'SVD');
plot(freqs, Gr_jacobian(:,1), 'DisplayName', 'Jacobian');
for in_i = 2:6
    set(gca,'ColorOrderIndex',1)
    plot(freqs, Gr_coupled(:,in_i), 'HandleVisibility', 'off');
    set(gca,'ColorOrderIndex',2)
    plot(freqs, Gr_decoupled(:,in_i), 'HandleVisibility', 'off');
    set(gca,'ColorOrderIndex',3)
    plot(freqs, Gr_jacobian(:,in_i), 'HandleVisibility', 'off');
end
plot(freqs, 0.5*ones(size(freqs)), 'k--', 'DisplayName', 'Limit')
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
hold off;
xlabel('Frequency (Hz)'); ylabel('Gershgorin Radii')
legend('location', 'northeast');

% Decoupled Plant
% Let's see the bode plot of the decoupled plant $G_d(s)$.
% \[ G_d(s) = U^T G_c(s) V \]


freqs = logspace(-1, 2, 1000);

figure;
hold on;
for ch_i = 1:6
  plot(freqs, abs(squeeze(freqresp(Gd(ch_i, ch_i), freqs, 'Hz'))), ...
       'DisplayName', sprintf('$G(%i, %i)$', ch_i, ch_i));
end
for in_i = 1:5
  for out_i = in_i+1:6
    plot(freqs, abs(squeeze(freqresp(Gd(out_i, in_i), freqs, 'Hz'))), 'color', [0, 0, 0, 0.2], ...
         'HandleVisibility', 'off');
  end
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Amplitude'); xlabel('Frequency [Hz]');
legend('location', 'southeast');



% #+name: fig:simscape_model_decoupled_plant_svd
% #+caption: Decoupled Plant using SVD
% #+RESULTS:
% [[file:figs/simscape_model_decoupled_plant_svd.png]]


freqs = logspace(-1, 2, 1000);

figure;
hold on;
for ch_i = 1:6
  plot(freqs, abs(squeeze(freqresp(Gj(ch_i, ch_i), freqs, 'Hz'))), ...
       'DisplayName', sprintf('$G(%i, %i)$', ch_i, ch_i));
end
for in_i = 1:5
  for out_i = in_i+1:6
    plot(freqs, abs(squeeze(freqresp(Gj(out_i, in_i), freqs, 'Hz'))), 'color', [0, 0, 0, 0.2], ...
         'HandleVisibility', 'off');
  end
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Amplitude'); xlabel('Frequency [Hz]');
legend('location', 'southeast');

% Diagonal Controller
% The controller $K$ is a diagonal controller consisting a low pass filters with a crossover frequency $\omega_c$ and a DC gain $C_g$.


wc = 2*pi*0.1; % Crossover Frequency [rad/s]
C_g = 50; % DC Gain

K = eye(6)*C_g/(s+wc);



% #+RESULTS:
% [[file:figs/centralized_control.png]]


G_cen = feedback(G, inv(J')*K, [7:12], [1:6]);



% #+RESULTS:
% [[file:figs/svd_control.png]]

% SVD Control

G_svd = feedback(G, pinv(V')*K*pinv(U), [7:12], [1:6]);

% Results
% Let's first verify the stability of the closed-loop systems:

isstable(G_cen)



% #+RESULTS:
% : ans =
% :   logical
% :    1


isstable(G_svd)



% #+RESULTS:
% : ans =
% :   logical
% :    0

% The obtained transmissibility in Open-loop, for the centralized control as well as for the SVD control are shown in Figure [[fig:stewart_platform_simscape_cl_transmissibility]].


freqs = logspace(-3, 3, 1000);

figure

ax1 = subplot(2, 3, 1);
hold on;
plot(freqs, abs(squeeze(freqresp(G(    'Ax', 'Dwx')/s^2, freqs, 'Hz'))), 'DisplayName', 'Open-Loop');
plot(freqs, abs(squeeze(freqresp(G_cen('Ax', 'Dwx')/s^2, freqs, 'Hz'))), 'DisplayName', 'Centralized');
plot(freqs, abs(squeeze(freqresp(G_svd('Ax', 'Dwx')/s^2, freqs, 'Hz'))), 'DisplayName', 'SVD');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Transmissibility - $D_x/D_{w,x}$');  xlabel('Frequency [Hz]');
legend('location', 'southwest');

ax2 = subplot(2, 3, 2);
hold on;
plot(freqs, abs(squeeze(freqresp(G(    'Ay', 'Dwy')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_cen('Ay', 'Dwy')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_svd('Ay', 'Dwy')/s^2, freqs, 'Hz'))));
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Transmissibility - $D_y/D_{w,y}$');  xlabel('Frequency [Hz]');

ax3 = subplot(2, 3, 3);
hold on;
plot(freqs, abs(squeeze(freqresp(G(    'Az', 'Dwz')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_cen('Az', 'Dwz')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_svd('Az', 'Dwz')/s^2, freqs, 'Hz'))));
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Transmissibility - $D_z/D_{w,z}$');  xlabel('Frequency [Hz]');

ax4 = subplot(2, 3, 4);
hold on;
plot(freqs, abs(squeeze(freqresp(G(    'Arx', 'Rwx')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_cen('Arx', 'Rwx')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_svd('Arx', 'Rwx')/s^2, freqs, 'Hz'))));
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Transmissibility - $R_x/R_{w,x}$');  xlabel('Frequency [Hz]');

ax5 = subplot(2, 3, 5);
hold on;
plot(freqs, abs(squeeze(freqresp(G(    'Ary', 'Rwy')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_cen('Ary', 'Rwy')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_svd('Ary', 'Rwy')/s^2, freqs, 'Hz'))));
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Transmissibility - $R_y/R_{w,y}$');  xlabel('Frequency [Hz]');

ax6 = subplot(2, 3, 6);
hold on;
plot(freqs, abs(squeeze(freqresp(G(    'Arz', 'Rwz')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_cen('Arz', 'Rwz')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_svd('Arz', 'Rwz')/s^2, freqs, 'Hz'))));
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Transmissibility - $R_z/R_{w,z}$');  xlabel('Frequency [Hz]');

linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'x');
xlim([freqs(1), freqs(end)]);
