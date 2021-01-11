%% Clear Workspace and Close figures
clear; close all; clc;

%% Intialize Laplace variable
s = zpk('s');

addpath('STEP');

freqs = logspace(-1, 2, 1000);

% Simscape Model - Parameters
% <<sec:stewart_simscape>>

open('drone_platform.slx');



% Definition of spring parameters:

kx = 0.5*1e3/3; % [N/m]
ky = 0.5*1e3/3;
kz = 1e3/3;

cx = 0.025; % [Nm/rad]
cy = 0.025;
cz = 0.025;



% We suppose the sensor is perfectly positioned.

sens_pos_error = zeros(3,1);



% Gravity:

g = 0;



% We load the Jacobian (previously computed from the geometry):

load('jacobian.mat', 'Aa', 'Ab', 'As', 'l', 'J');



% We initialize other parameters:

U = eye(6);
V = eye(6);
Kc = tf(zeros(6));



% #+name: fig:stewart_platform_plant
% #+caption: Considered plant $\bm{G} = \begin{bmatrix}G_d\\G_u\end{bmatrix}$. $D_w$ is the translation/rotation of the support, $\tau$ the actuator forces, $a$ the acceleration/angular acceleration of the top platform
% #+RESULTS:
% [[file:figs/stewart_platform_plant.png]]


%% Name of the Simulink File
mdl = 'drone_platform';

%% Input/Output definition
clear io; io_i = 1;
io(io_i) = linio([mdl, '/Dw'],              1, 'openinput');  io_i = io_i + 1; % Ground Motion
io(io_i) = linio([mdl, '/V-T'],             1, 'openinput');  io_i = io_i + 1; % Actuator Forces
io(io_i) = linio([mdl, '/Inertial Sensor'], 1, 'openoutput'); io_i = io_i + 1; % Top platform acceleration

G = linearize(mdl, io);
G.InputName  = {'Dwx', 'Dwy', 'Dwz', 'Rwx', 'Rwy', 'Rwz', ...
                'F1', 'F2', 'F3', 'F4', 'F5', 'F6'};
G.OutputName = {'Ax', 'Ay', 'Az', 'Arx', 'Ary', 'Arz'};

% Plant
Gu = G(:, {'F1', 'F2', 'F3', 'F4', 'F5', 'F6'});
% Disturbance dynamics
Gd = G(:, {'Dwx', 'Dwy', 'Dwz', 'Rwx', 'Rwy', 'Rwz'});



% There are 24 states (6dof for the bottom platform + 6dof for the top platform).

size(G)



% #+RESULTS:
% : State-space model with 6 outputs, 12 inputs, and 24 states.

% The elements of the transfer matrix $\bm{G}$ corresponding to the transfer function from actuator forces $\tau$ to the measured acceleration $a$ are shown in Figure [[fig:stewart_platform_coupled_plant]].

% One can easily see that the system is strongly coupled.


figure;

% Magnitude
hold on;
for i_in = 1:6
    for i_out = [1:i_in-1, i_in+1:6]
        plot(freqs, abs(squeeze(freqresp(Gu(i_out, i_in), freqs, 'Hz'))), 'color', [0,0,0,0.2], ...
             'HandleVisibility', 'off');
    end
end
plot(freqs, abs(squeeze(freqresp(Gu(i_out, i_in), freqs, 'Hz'))), 'color', [0,0,0,0.2], ...
     'DisplayName', '$G_u(i,j)\ i \neq j$');
set(gca,'ColorOrderIndex',1)
for i_in_out = 1:6
    plot(freqs, abs(squeeze(freqresp(Gu(i_in_out, i_in_out), freqs, 'Hz'))), 'DisplayName', sprintf('$G_u(%d,%d)$', i_in_out, i_in_out));
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
xlabel('Frequency [Hz]'); ylabel('Magnitude');
ylim([1e-2, 1e5]);
legend('location', 'northwest');



% #+name: fig:plant_decouple_jacobian
% #+caption: Decoupled plant $\bm{G}_x$ using the Jacobian matrix $J$
% #+RESULTS:
% [[file:figs/plant_decouple_jacobian.png]]

% We define a new plant:
% \[ G_x(s) = G(s) J^{-T} \]

% $G_x(s)$ correspond to the transfer function from forces and torques applied to the top platform to the absolute acceleration of the top platform.


Gx = Gu*inv(J');
Gx.InputName  = {'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'};

% Decoupling using the SVD
% <<sec:stewart_svd_decoupling>>

% In order to decouple the plant using the SVD, first a real approximation of the plant transfer function matrix as the crossover frequency is required.

% Let's compute a real approximation of the complex matrix $H_1$ which corresponds to the the transfer function $G_u(j\omega_c)$ from forces applied by the actuators to the measured acceleration of the top platform evaluated at the frequency $\omega_c$.

wc = 2*pi*30; % Decoupling frequency [rad/s]

H1 = evalfr(Gu, j*wc);



% The real approximation is computed as follows:

D = pinv(real(H1'*H1));
H1 = inv(D*real(H1'*diag(exp(j*angle(diag(H1*D*H1.'))/2))));



% #+caption: Real part of $G$ at the decoupling frequency $\omega_c$
% #+RESULTS:
% |    4.4 |   -2.1 |   -2.1 |    4.4 |  -2.4 |   -2.4 |
% |   -0.2 |   -3.9 |    3.9 |    0.2 |  -3.8 |    3.8 |
% |    3.4 |    3.4 |    3.4 |    3.4 |   3.4 |    3.4 |
% | -367.1 | -323.8 |  323.8 |  367.1 |  43.3 |  -43.3 |
% | -162.0 | -237.0 | -237.0 | -162.0 | 398.9 |  398.9 |
% |  220.6 | -220.6 |  220.6 | -220.6 | 220.6 | -220.6 |

% Now, the Singular Value Decomposition of $H_1$ is performed:
% \[ H_1 = U \Sigma V^H \]


[U,~,V] = svd(H1);



% #+name: fig:plant_decouple_svd
% #+caption: Decoupled plant $\bm{G}_{SVD}$ using the Singular Value Decomposition
% #+RESULTS:
% [[file:figs/plant_decouple_svd.png]]

% The decoupled plant is then:
% \[ G_{SVD}(s) = U^{-1} G_u(s) V^{-H} \]


Gsvd = inv(U)*Gu*inv(V');

% Verification of the decoupling using the "Gershgorin Radii"
% <<sec:stewart_gershorin_radii>>

% The "Gershgorin Radii" is computed for the coupled plant $G(s)$, for the "Jacobian plant" $G_x(s)$ and the "SVD Decoupled Plant" $G_{SVD}(s)$:

% The "Gershgorin Radii" of a matrix $S$ is defined by:
% \[ \zeta_i(j\omega) = \frac{\sum\limits_{j\neq i}|S_{ij}(j\omega)|}{|S_{ii}(j\omega)|} \]

% This is computed over the following frequencies.

% Gershgorin Radii for the coupled plant:
Gr_coupled = zeros(length(freqs), size(Gu,2));
H = abs(squeeze(freqresp(Gu, freqs, 'Hz')));
for out_i = 1:size(Gu,2)
    Gr_coupled(:, out_i) = squeeze((sum(H(out_i,:,:)) - H(out_i,out_i,:))./H(out_i, out_i, :));
end

% Gershgorin Radii for the decoupled plant using SVD:
Gr_decoupled = zeros(length(freqs), size(Gsvd,2));
H = abs(squeeze(freqresp(Gsvd, freqs, 'Hz')));
for out_i = 1:size(Gsvd,2)
    Gr_decoupled(:, out_i) = squeeze((sum(H(out_i,:,:)) - H(out_i,out_i,:))./H(out_i, out_i, :));
end

% Gershgorin Radii for the decoupled plant using the Jacobian:
Gr_jacobian = zeros(length(freqs), size(Gx,2));
H = abs(squeeze(freqresp(Gx, freqs, 'Hz')));
for out_i = 1:size(Gx,2)
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
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
hold off;
xlabel('Frequency (Hz)'); ylabel('Gershgorin Radii')
legend('location', 'northwest');
ylim([1e-3, 1e3]);

% Verification of the decoupling using the "Relative Gain Array"
% <<sec:stewart_rga>>

% The relative gain array (RGA) is defined as:
% \begin{equation}
%   \Lambda\big(G(s)\big) = G(s) \times \big( G(s)^{-1} \big)^T
% \end{equation}
% where $\times$ denotes an element by element multiplication and $G(s)$ is an $n \times n$ square transfer matrix.

% The obtained RGA elements are shown in Figure [[fig:simscape_model_rga]].


% Relative Gain Array for the coupled plant:
RGA_coupled = zeros(length(freqs), size(Gu,1), size(Gu,2));
Gu_inv = inv(Gu);
for f_i = 1:length(freqs)
    RGA_coupled(f_i, :, :) = abs(evalfr(Gu, j*2*pi*freqs(f_i)).*evalfr(Gu_inv, j*2*pi*freqs(f_i))');
end

% Relative Gain Array for the decoupled plant using SVD:
RGA_svd = zeros(length(freqs), size(Gsvd,1), size(Gsvd,2));
Gsvd_inv = inv(Gsvd);
for f_i = 1:length(freqs)
    RGA_svd(f_i, :, :) = abs(evalfr(Gsvd, j*2*pi*freqs(f_i)).*evalfr(Gsvd_inv, j*2*pi*freqs(f_i))');
end

% Relative Gain Array for the decoupled plant using the Jacobian:
RGA_x = zeros(length(freqs), size(Gx,1), size(Gx,2));
Gx_inv = inv(Gx);
for f_i = 1:length(freqs)
    RGA_x(f_i, :, :) = abs(evalfr(Gx, j*2*pi*freqs(f_i)).*evalfr(Gx_inv, j*2*pi*freqs(f_i))');
end

figure;
tiledlayout(1, 2, 'TileSpacing', 'None', 'Padding', 'None');

ax1 = nexttile;
hold on;
for i_in = 1:6
    for i_out = [1:i_in-1, i_in+1:6]
        plot(freqs, RGA_svd(:, i_out, i_in), '--', 'color', [0 0 0 0.2], ...
             'HandleVisibility', 'off');
    end
end
plot(freqs, RGA_svd(:, 1, 2), '--', 'color', [0 0 0 0.2], ...
     'DisplayName', '$RGA_{SVD}(i,j),\ i \neq j$');

plot(freqs, RGA_svd(:, 1, 1), 'k-', ...
     'DisplayName', '$RGA_{SVD}(i,i)$');
for ch_i = 1:6
    plot(freqs, RGA_svd(:, ch_i, ch_i), 'k-', ...
         'HandleVisibility', 'off');
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Magnitude'); xlabel('Frequency [Hz]');
legend('location', 'southwest');

ax2 = nexttile;
hold on;
for i_in = 1:6
    for i_out = [1:i_in-1, i_in+1:6]
        plot(freqs, RGA_x(:, i_out, i_in), '--', 'color', [0 0 0 0.2], ...
             'HandleVisibility', 'off');
    end
end
plot(freqs, RGA_x(:, 1, 2), '--', 'color', [0 0 0 0.2], ...
     'DisplayName', '$RGA_{X}(i,j),\ i \neq j$');

plot(freqs, RGA_x(:, 1, 1), 'k-', ...
     'DisplayName', '$RGA_{X}(i,i)$');
for ch_i = 1:6
    plot(freqs, RGA_x(:, ch_i, ch_i), 'k-', ...
         'HandleVisibility', 'off');
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
xlabel('Frequency [Hz]'); set(gca, 'YTickLabel',[]);
legend('location', 'southwest');

linkaxes([ax1,ax2],'y');
ylim([1e-5, 1e1]);

% Obtained Decoupled Plants
% <<sec:stewart_decoupled_plant>>

% The bode plot of the diagonal and off-diagonal elements of $G_{SVD}$ are shown in Figure [[fig:simscape_model_decoupled_plant_svd]].


figure;
tiledlayout(3, 1, 'TileSpacing', 'None', 'Padding', 'None');

% Magnitude
ax1 = nexttile([2, 1]);
hold on;
for i_in = 1:6
    for i_out = [1:i_in-1, i_in+1:6]
        plot(freqs, abs(squeeze(freqresp(Gsvd(i_out, i_in), freqs, 'Hz'))), 'color', [0,0,0,0.2], ...
             'HandleVisibility', 'off');
    end
end
plot(freqs, abs(squeeze(freqresp(Gsvd(1, 2), freqs, 'Hz'))), 'color', [0,0,0,0.5], ...
     'DisplayName', '$G_{SVD}(i,j),\ i \neq j$');
set(gca,'ColorOrderIndex',1)
for ch_i = 1:6
    plot(freqs, abs(squeeze(freqresp(Gsvd(ch_i, ch_i), freqs, 'Hz'))), ...
         'DisplayName', sprintf('$G_{SVD}(%i,%i)$', ch_i, ch_i));
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Magnitude'); set(gca, 'XTickLabel',[]);
legend('location', 'northwest');
ylim([1e-1, 1e5])

% Phase
ax2 = nexttile;
hold on;
for ch_i = 1:6
    plot(freqs, 180/pi*angle(squeeze(freqresp(Gsvd(ch_i, ch_i), freqs, 'Hz'))));
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'lin');
ylabel('Phase [deg]'); xlabel('Frequency [Hz]');
ylim([-180, 180]);
yticks([-180:90:360]);

linkaxes([ax1,ax2],'x');



% #+name: fig:simscape_model_decoupled_plant_svd
% #+caption: Decoupled Plant using SVD
% #+RESULTS:
% [[file:figs/simscape_model_decoupled_plant_svd.png]]

% Similarly, the bode plots of the diagonal elements and off-diagonal elements of the decoupled plant $G_x(s)$ using the Jacobian are shown in Figure [[fig:simscape_model_decoupled_plant_jacobian]].


figure;
tiledlayout(3, 1, 'TileSpacing', 'None', 'Padding', 'None');

% Magnitude
ax1 = nexttile([2, 1]);
hold on;
for i_in = 1:6
    for i_out = [1:i_in-1, i_in+1:6]
        plot(freqs, abs(squeeze(freqresp(Gx(i_out, i_in), freqs, 'Hz'))), 'color', [0,0,0,0.2], ...
             'HandleVisibility', 'off');
    end
end
plot(freqs, abs(squeeze(freqresp(Gx(1, 2), freqs, 'Hz'))), 'color', [0,0,0,0.5], ...
     'DisplayName', '$G_x(i,j),\ i \neq j$');
set(gca,'ColorOrderIndex',1)
plot(freqs, abs(squeeze(freqresp(Gx('Ax', 'Fx'), freqs, 'Hz'))), 'DisplayName', '$G_x(1,1) = A_x/F_x$');
plot(freqs, abs(squeeze(freqresp(Gx('Ay', 'Fy'), freqs, 'Hz'))), 'DisplayName', '$G_x(2,2) = A_y/F_y$');
plot(freqs, abs(squeeze(freqresp(Gx('Az', 'Fz'), freqs, 'Hz'))), 'DisplayName', '$G_x(3,3) = A_z/F_z$');
plot(freqs, abs(squeeze(freqresp(Gx('Arx', 'Mx'), freqs, 'Hz'))), 'DisplayName', '$G_x(4,4) = A_{R_x}/M_x$');
plot(freqs, abs(squeeze(freqresp(Gx('Ary', 'My'), freqs, 'Hz'))), 'DisplayName', '$G_x(5,5) = A_{R_y}/M_y$');
plot(freqs, abs(squeeze(freqresp(Gx('Arz', 'Mz'), freqs, 'Hz'))), 'DisplayName', '$G_x(6,6) = A_{R_z}/M_z$');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Magnitude'); set(gca, 'XTickLabel',[]);
legend('location', 'northwest');
ylim([1e-2, 2e6])

% Phase
ax2 = nexttile;
hold on;
plot(freqs, 180/pi*angle(squeeze(freqresp(Gx('Ax', 'Fx'), freqs, 'Hz'))));
plot(freqs, 180/pi*angle(squeeze(freqresp(Gx('Ay', 'Fy'), freqs, 'Hz'))));
plot(freqs, 180/pi*angle(squeeze(freqresp(Gx('Az', 'Fz'), freqs, 'Hz'))));
plot(freqs, 180/pi*angle(squeeze(freqresp(Gx('Arx', 'Mx'), freqs, 'Hz'))));
plot(freqs, 180/pi*angle(squeeze(freqresp(Gx('Ary', 'My'), freqs, 'Hz'))));
plot(freqs, 180/pi*angle(squeeze(freqresp(Gx('Arz', 'Mz'), freqs, 'Hz'))));
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'lin');
ylabel('Phase [deg]'); xlabel('Frequency [Hz]');
ylim([0, 180]);
yticks([0:45:360]);

linkaxes([ax1,ax2],'x');



% #+name: fig:svd_control
% #+caption: Control Diagram for the SVD control
% #+RESULTS:
% [[file:figs/svd_control.png]]


% We choose the controller to be a low pass filter:
% \[ K_c(s) = \frac{G_0}{1 + \frac{s}{\omega_0}} \]

% $G_0$ is tuned such that the crossover frequency corresponding to the diagonal terms of the loop gain is equal to $\omega_c$


wc = 2*pi*80;  % Crossover Frequency [rad/s]
w0 = 2*pi*0.1; % Controller Pole [rad/s]

K_cen = diag(1./diag(abs(evalfr(Gx, j*wc))))*(1/abs(evalfr(1/(1 + s/w0), j*wc)))/(1 + s/w0);
L_cen = K_cen*Gx;
G_cen = feedback(G, pinv(J')*K_cen, [7:12], [1:6]);

K_svd = diag(1./diag(abs(evalfr(Gsvd, j*wc))))*(1/abs(evalfr(1/(1 + s/w0), j*wc)))/(1 + s/w0);
L_svd = K_svd*Gsvd;
G_svd = feedback(G, inv(V')*K_svd*inv(U), [7:12], [1:6]);



% The obtained diagonal elements of the loop gains are shown in Figure [[fig:stewart_comp_loop_gain_diagonal]].


figure;
tiledlayout(3, 1, 'TileSpacing', 'None', 'Padding', 'None');

% Magnitude
ax1 = nexttile([2, 1]);
hold on;
plot(freqs, abs(squeeze(freqresp(L_svd(1, 1), freqs, 'Hz'))), 'DisplayName', '$L_{SVD}(i,i)$');
for i_in_out = 2:6
    set(gca,'ColorOrderIndex',1)
    plot(freqs, abs(squeeze(freqresp(L_svd(i_in_out, i_in_out), freqs, 'Hz'))), 'HandleVisibility', 'off');
end

set(gca,'ColorOrderIndex',2)
plot(freqs, abs(squeeze(freqresp(L_cen(1, 1), freqs, 'Hz'))), ...
     'DisplayName', '$L_{J}(i,i)$');
for i_in_out = 2:6
    set(gca,'ColorOrderIndex',2)
    plot(freqs, abs(squeeze(freqresp(L_cen(i_in_out, i_in_out), freqs, 'Hz'))), 'HandleVisibility', 'off');
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Magnitude'); set(gca, 'XTickLabel',[]);
legend('location', 'northwest');
ylim([5e-2, 2e3])

% Phase
ax2 = nexttile;
hold on;
for i_in_out = 1:6
    set(gca,'ColorOrderIndex',1)
    plot(freqs, 180/pi*angle(squeeze(freqresp(L_svd(i_in_out, i_in_out), freqs, 'Hz'))));
end
set(gca,'ColorOrderIndex',2)
for i_in_out = 1:6
    set(gca,'ColorOrderIndex',2)
    plot(freqs, 180/pi*angle(squeeze(freqresp(L_cen(i_in_out, i_in_out), freqs, 'Hz'))));
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'lin');
ylabel('Phase [deg]'); xlabel('Frequency [Hz]');
ylim([-180, 180]);
yticks([-180:90:360]);

linkaxes([ax1,ax2],'x');

% Closed-Loop system Performances
% <<sec:stewart_closed_loop_results>>

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
% :    1

% The obtained transmissibility in Open-loop, for the centralized control as well as for the SVD control are shown in Figure [[fig:stewart_platform_simscape_cl_transmissibility]].


figure;
tiledlayout(2, 2, 'TileSpacing', 'None', 'Padding', 'None');

ax1 = nexttile;
hold on;
plot(freqs, abs(squeeze(freqresp(G(    'Ax', 'Dwx')/s^2, freqs, 'Hz'))), 'DisplayName', 'Open-Loop');
plot(freqs, abs(squeeze(freqresp(G_cen('Ax', 'Dwx')/s^2, freqs, 'Hz'))), 'DisplayName', 'Centralized');
plot(freqs, abs(squeeze(freqresp(G_svd('Ax', 'Dwx')/s^2, freqs, 'Hz'))), '--', 'DisplayName', 'SVD');
set(gca,'ColorOrderIndex',1)
plot(freqs, abs(squeeze(freqresp(G(    'Ay', 'Dwy')/s^2, freqs, 'Hz'))), 'HandleVisibility', 'off');
plot(freqs, abs(squeeze(freqresp(G_cen('Ay', 'Dwy')/s^2, freqs, 'Hz'))), 'HandleVisibility', 'off');
plot(freqs, abs(squeeze(freqresp(G_svd('Ay', 'Dwy')/s^2, freqs, 'Hz'))), '--', 'HandleVisibility', 'off');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('$D_x/D_{w,x}$, $D_y/D_{w, y}$'); set(gca, 'XTickLabel',[]);
legend('location', 'southwest');

ax2 = nexttile;
hold on;
plot(freqs, abs(squeeze(freqresp(G(    'Az', 'Dwz')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_cen('Az', 'Dwz')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_svd('Az', 'Dwz')/s^2, freqs, 'Hz'))), '--');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('$D_z/D_{w,z}$'); set(gca, 'XTickLabel',[]);

ax3 = nexttile;
hold on;
plot(freqs, abs(squeeze(freqresp(G(    'Arx', 'Rwx')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_cen('Arx', 'Rwx')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_svd('Arx', 'Rwx')/s^2, freqs, 'Hz'))), '--');
set(gca,'ColorOrderIndex',1)
plot(freqs, abs(squeeze(freqresp(G(    'Ary', 'Rwy')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_cen('Ary', 'Rwy')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_svd('Ary', 'Rwy')/s^2, freqs, 'Hz'))), '--');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('$R_x/R_{w,x}$, $R_y/R_{w,y}$');  xlabel('Frequency [Hz]');

ax4 = nexttile;
hold on;
plot(freqs, abs(squeeze(freqresp(G(    'Arz', 'Rwz')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_cen('Arz', 'Rwz')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_svd('Arz', 'Rwz')/s^2, freqs, 'Hz'))), '--');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('$R_z/R_{w,z}$');  xlabel('Frequency [Hz]');

linkaxes([ax1,ax2,ax3,ax4],'xy');
xlim([freqs(1), freqs(end)]);
ylim([1e-3, 1e2]);
