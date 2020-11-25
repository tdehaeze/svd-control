%% Clear Workspace and Close figures
clear; close all; clc;

%% Intialize Laplace variable
s = zpk('s');

freqs = logspace(-1, 2, 1000);

% Gravimeter Model - Parameters
% <<sec:gravimeter_model>>


open('gravimeter.slx')



% The model of the gravimeter is schematically shown in Figure [[fig:gravimeter_model]].

% #+name: fig:gravimeter_model
% #+caption: Model of the gravimeter
% [[file:figs/gravimeter_model.png]]

% The parameters used for the simulation are the following:

l  = 1.0; % Length of the mass [m]
h  = 1.7; % Height of the mass [m]

la = l/2; % Position of Act. [m]
ha = h/2; % Position of Act. [m]

m = 400; % Mass [kg]
I = 115; % Inertia [kg m^2]

k = 15e3; % Actuator Stiffness [N/m]
c = 2e1; % Actuator Damping [N/(m/s)]

deq = 0.2; % Length of the actuators [m]

g = 0; % Gravity [m/s2]

% System Identification
% <<sec:gravimeter_identification>>


%% Name of the Simulink File
mdl = 'gravimeter';

%% Input/Output definition
clear io; io_i = 1;
io(io_i) = linio([mdl, '/F1'], 1, 'openinput');  io_i = io_i + 1;
io(io_i) = linio([mdl, '/F2'], 1, 'openinput');  io_i = io_i + 1;
io(io_i) = linio([mdl, '/F3'], 1, 'openinput');  io_i = io_i + 1;
io(io_i) = linio([mdl, '/Acc_side'], 1, 'openoutput'); io_i = io_i + 1;
io(io_i) = linio([mdl, '/Acc_side'], 2, 'openoutput'); io_i = io_i + 1;
io(io_i) = linio([mdl, '/Acc_top'], 1, 'openoutput'); io_i = io_i + 1;
io(io_i) = linio([mdl, '/Acc_top'], 2, 'openoutput'); io_i = io_i + 1;

G = linearize(mdl, io);
G.InputName  = {'F1', 'F2', 'F3'};
G.OutputName = {'Ax1', 'Az1', 'Ax2', 'Az2'};



% #+name: fig:gravimeter_plant_schematic
% #+caption: Schematic of the gravimeter plant
% #+RESULTS:
% [[file:figs/gravimeter_plant_schematic.png]]

% We can check the poles of the plant:

pole(G)



% #+RESULTS:
% | -0.12243+13.551i   |
% | -0.12243-13.551i   |
% | -0.05+8.6601i      |
% | -0.05-8.6601i      |
% | -0.0088785+3.6493i |
% | -0.0088785-3.6493i |

% As expected, the plant as 6 states (2 translations + 1 rotation)

size(G)



% #+RESULTS:
% : State-space model with 4 outputs, 3 inputs, and 6 states.

% The bode plot of all elements of the plant are shown in Figure [[fig:open_loop_tf]].


figure;
tiledlayout(4, 3, 'TileSpacing', 'None', 'Padding', 'None');

for out_i = 1:4
    for in_i = 1:3
        nexttile;
        plot(freqs, abs(squeeze(freqresp(G(out_i,in_i), freqs, 'Hz'))), '-');
        set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
        xlim([1e-1, 2e1]); ylim([1e-4, 1e0]);

        if in_i == 1
            ylabel('Amplitude [m/N]')
        else
            set(gca, 'YTickLabel',[]);
        end

        if out_i == 4
            xlabel('Frequency [Hz]')
        else
            set(gca, 'XTickLabel',[]);
        end
    end
end



% #+name: fig:gravimeter_decouple_jacobian
% #+caption: Decoupled plant $\bm{G}_x$ using the Jacobian matrix $J$
% #+RESULTS:
% [[file:figs/gravimeter_decouple_jacobian.png]]

% The Jacobian corresponding to the sensors and actuators are defined below:

Ja = [1 0  h/2
      0 1 -l/2
      1 0 -h/2
      0 1  0];

Jt = [1 0  ha
      0 1 -la
      0 1  la];



% And the plant $\bm{G}_x$ is computed:

Gx = pinv(Ja)*G*pinv(Jt');
Gx.InputName  = {'Fx', 'Fz', 'My'};
Gx.OutputName  = {'Dx', 'Dz', 'Ry'};

size(Gx)



% #+RESULTS:
% : size(Gx)
% : State-space model with 3 outputs, 3 inputs, and 6 states.

% The diagonal and off-diagonal elements of $G_x$ are shown in Figure [[fig:gravimeter_jacobian_plant]].


figure;

% Magnitude
hold on;
for i_in = 1:3
    for i_out = [1:i_in-1, i_in+1:3]
        plot(freqs, abs(squeeze(freqresp(Gx(i_out, i_in), freqs, 'Hz'))), 'color', [0,0,0,0.2], ...
             'HandleVisibility', 'off');
    end
end
plot(freqs, abs(squeeze(freqresp(Gx(i_out, i_in), freqs, 'Hz'))), 'color', [0,0,0,0.2], ...
     'DisplayName', '$G_x(i,j)\ i \neq j$');
set(gca,'ColorOrderIndex',1)
for i_in_out = 1:3
  plot(freqs, abs(squeeze(freqresp(Gx(i_in_out, i_in_out), freqs, 'Hz'))), 'DisplayName', sprintf('$G_x(%d,%d)$', i_in_out, i_in_out));
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
xlabel('Frequency [Hz]'); ylabel('Magnitude');
legend('location', 'southeast');
ylim([1e-8, 1e0]);

% Decoupling using the SVD
% <<sec:gravimeter_svd_decoupling>>

% In order to decouple the plant using the SVD, first a real approximation of the plant transfer function matrix as the crossover frequency is required.

% Let's compute a real approximation of the complex matrix $H_1$ which corresponds to the the transfer function $G(j\omega_c)$ from forces applied by the actuators to the measured acceleration of the top platform evaluated at the frequency $\omega_c$.

wc = 2*pi*10; % Decoupling frequency [rad/s]

H1 = evalfr(G, j*wc);



% The real approximation is computed as follows:

D = pinv(real(H1'*H1));
H1 = pinv(D*real(H1'*diag(exp(j*angle(diag(H1*D*H1.'))/2))));



% #+caption: Real approximate of $G$ at the decoupling frequency $\omega_c$
% #+RESULTS:
% |  0.0092 | -0.0039 |  0.0039 |
% | -0.0039 |  0.0048 | 0.00028 |
% |  -0.004 |  0.0038 | -0.0038 |
% | 8.4e-09 |  0.0025 |  0.0025 |


% Now, the Singular Value Decomposition of $H_1$ is performed:
% \[ H_1 = U \Sigma V^H \]


[U,S,V] = svd(H1);



% #+name: fig:gravimeter_decouple_svd
% #+caption: Decoupled plant $\bm{G}_{SVD}$ using the Singular Value Decomposition
% #+RESULTS:
% [[file:figs/gravimeter_decouple_svd.png]]

% The decoupled plant is then:
% \[ \bm{G}_{SVD}(s) = U^{-1} \bm{G}(s) V^{-H} \]


Gsvd = inv(U)*G*inv(V');

size(Gsvd)



% #+RESULTS:
% : size(Gsvd)
% : State-space model with 4 outputs, 3 inputs, and 6 states.

% The 4th output (corresponding to the null singular value) is discarded, and we only keep the $3 \times 3$ plant:

Gsvd = Gsvd(1:3, 1:3);



% The diagonal and off-diagonal elements of the "SVD" plant are shown in Figure [[fig:gravimeter_svd_plant]].

figure;

% Magnitude
hold on;
for i_in = 1:3
    for i_out = [1:i_in-1, i_in+1:3]
        plot(freqs, abs(squeeze(freqresp(Gsvd(i_out, i_in), freqs, 'Hz'))), 'color', [0,0,0,0.2], ...
             'HandleVisibility', 'off');
    end
end
plot(freqs, abs(squeeze(freqresp(Gsvd(i_out, i_in), freqs, 'Hz'))), 'color', [0,0,0,0.2], ...
     'DisplayName', '$G_x(i,j)\ i \neq j$');
set(gca,'ColorOrderIndex',1)
for i_in_out = 1:3
  plot(freqs, abs(squeeze(freqresp(Gsvd(i_in_out, i_in_out), freqs, 'Hz'))), 'DisplayName', sprintf('$G_x(%d,%d)$', i_in_out, i_in_out));
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
xlabel('Frequency [Hz]'); ylabel('Magnitude');
legend('location', 'southwest', 'FontSize', 8);
ylim([1e-8, 1e0]);

% Verification of the decoupling using the "Gershgorin Radii"
% <<sec:gravimeter_gershgorin_radii>>

% The "Gershgorin Radii" is computed for the coupled plant $G(s)$, for the "Jacobian plant" $G_x(s)$ and the "SVD Decoupled Plant" $G_{SVD}(s)$:

% The "Gershgorin Radii" of a matrix $S$ is defined by:
% \[ \zeta_i(j\omega) = \frac{\sum\limits_{j\neq i}|S_{ij}(j\omega)|}{|S_{ii}(j\omega)|} \]


% Gershgorin Radii for the coupled plant:
Gr_coupled = zeros(length(freqs), size(G,2));
H = abs(squeeze(freqresp(G, freqs, 'Hz')));
for out_i = 1:size(G,2)
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
for in_i = 2:3
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
legend('location', 'southwest');
ylim([1e-4, 1e2]);

% Verification of the decoupling using the "Relative Gain Array"
% <<sec:gravimeter_rga>>

% The relative gain array (RGA) is defined as:
% \begin{equation}
%   \Lambda\big(G(s)\big) = G(s) \times \big( G(s)^{-1} \big)^T
% \end{equation}
% where $\times$ denotes an element by element multiplication and $G(s)$ is an $n \times n$ square transfer matrix.

% The obtained RGA elements are shown in Figure [[fig:gravimeter_rga]].


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
for i_in = 1:3
    for i_out = [1:i_in-1, i_in+1:3]
        plot(freqs, RGA_svd(:, i_out, i_in), '--', 'color', [0 0 0 0.2], ...
             'HandleVisibility', 'off');
    end
end
plot(freqs, RGA_svd(:, 1, 2), '--', 'color', [0 0 0 0.2], ...
     'DisplayName', '$RGA_{SVD}(i,j),\ i \neq j$');

plot(freqs, RGA_svd(:, 1, 1), 'k-', ...
     'DisplayName', '$RGA_{SVD}(i,i)$');
for ch_i = 1:3
  plot(freqs, RGA_svd(:, ch_i, ch_i), 'k-', ...
       'HandleVisibility', 'off');
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Magnitude'); xlabel('Frequency [Hz]');
legend('location', 'southwest');

ax2 = nexttile;
hold on;
for i_in = 1:3
    for i_out = [1:i_in-1, i_in+1:3]
        plot(freqs, RGA_x(:, i_out, i_in), '--', 'color', [0 0 0 0.2], ...
             'HandleVisibility', 'off');
    end
end
plot(freqs, RGA_x(:, 1, 2), '--', 'color', [0 0 0 0.2], ...
     'DisplayName', '$RGA_{X}(i,j),\ i \neq j$');

plot(freqs, RGA_x(:, 1, 1), 'k-', ...
     'DisplayName', '$RGA_{X}(i,i)$');
for ch_i = 1:3
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
% <<sec:gravimeter_decoupled_plant>>

% The bode plot of the diagonal and off-diagonal elements of $G_{SVD}$ are shown in Figure [[fig:gravimeter_decoupled_plant_svd]].


figure;
tiledlayout(3, 1, 'TileSpacing', 'None', 'Padding', 'None');

% Magnitude
ax1 = nexttile([2, 1]);
hold on;
for i_in = 1:3
    for i_out = [1:i_in-1, i_in+1:3]
        plot(freqs, abs(squeeze(freqresp(Gsvd(i_out, i_in), freqs, 'Hz'))), 'color', [0,0,0,0.2], ...
             'HandleVisibility', 'off');
    end
end
plot(freqs, abs(squeeze(freqresp(Gsvd(1, 2), freqs, 'Hz'))), 'color', [0,0,0,0.5], ...
     'DisplayName', '$G_{SVD}(i,j),\ i \neq j$');
set(gca,'ColorOrderIndex',1)
for ch_i = 1:3
  plot(freqs, abs(squeeze(freqresp(Gsvd(ch_i, ch_i), freqs, 'Hz'))), ...
       'DisplayName', sprintf('$G_{SVD}(%i,%i)$', ch_i, ch_i));
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Magnitude'); set(gca, 'XTickLabel',[]);
legend('location', 'southwest');
ylim([1e-8, 1e0])

% Phase
ax2 = nexttile;
hold on;
for ch_i = 1:3
  plot(freqs, 180/pi*angle(squeeze(freqresp(Gsvd(ch_i, ch_i), freqs, 'Hz'))));
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'lin');
ylabel('Phase [deg]'); xlabel('Frequency [Hz]');
ylim([-180, 180]);
yticks([-180:90:360]);

linkaxes([ax1,ax2],'x');



% #+name: fig:gravimeter_decoupled_plant_svd
% #+caption: Decoupled Plant using SVD
% #+RESULTS:
% [[file:figs/gravimeter_decoupled_plant_svd.png]]

% Similarly, the bode plots of the diagonal elements and off-diagonal elements of the decoupled plant $G_x(s)$ using the Jacobian are shown in Figure [[fig:gravimeter_decoupled_plant_jacobian]].


figure;
tiledlayout(3, 1, 'TileSpacing', 'None', 'Padding', 'None');

% Magnitude
ax1 = nexttile([2, 1]);
hold on;
for i_in = 1:3
    for i_out = [1:i_in-1, i_in+1:3]
        plot(freqs, abs(squeeze(freqresp(Gx(i_out, i_in), freqs, 'Hz'))), 'color', [0,0,0,0.2], ...
             'HandleVisibility', 'off');
    end
end
plot(freqs, abs(squeeze(freqresp(Gx(1, 2), freqs, 'Hz'))), 'color', [0,0,0,0.5], ...
     'DisplayName', '$G_x(i,j),\ i \neq j$');
set(gca,'ColorOrderIndex',1)
plot(freqs, abs(squeeze(freqresp(Gx(1, 1), freqs, 'Hz'))), 'DisplayName', '$G_x(1,1) = A_x/F_x$');
plot(freqs, abs(squeeze(freqresp(Gx(2, 2), freqs, 'Hz'))), 'DisplayName', '$G_x(2,2) = A_y/F_y$');
plot(freqs, abs(squeeze(freqresp(Gx(3, 3), freqs, 'Hz'))), 'DisplayName', '$G_x(3,3) = R_y/M_y$');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Magnitude'); set(gca, 'XTickLabel',[]);
legend('location', 'southwest');
ylim([1e-8, 1e0])

% Phase
ax2 = nexttile;
hold on;
plot(freqs, 180/pi*angle(squeeze(freqresp(Gx(1, 1), freqs, 'Hz'))));
plot(freqs, 180/pi*angle(squeeze(freqresp(Gx(2, 2), freqs, 'Hz'))));
plot(freqs, 180/pi*angle(squeeze(freqresp(Gx(3, 3), freqs, 'Hz'))));
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'lin');
ylabel('Phase [deg]'); xlabel('Frequency [Hz]');
ylim([-180, 180]);
yticks([0:45:360]);

linkaxes([ax1,ax2],'x');



% #+name: fig:svd_control_gravimeter
% #+caption: Control Diagram for the SVD control
% #+RESULTS:
% [[file:figs/svd_control_gravimeter.png]]


% We choose the controller to be a low pass filter:
% \[ K_c(s) = \frac{G_0}{1 + \frac{s}{\omega_0}} \]

% $G_0$ is tuned such that the crossover frequency corresponding to the diagonal terms of the loop gain is equal to $\omega_c$


wc = 2*pi*10;  % Crossover Frequency [rad/s]
w0 = 2*pi*0.1; % Controller Pole [rad/s]

K_cen = diag(1./diag(abs(evalfr(Gx, j*wc))))*(1/abs(evalfr(1/(1 + s/w0), j*wc)))/(1 + s/w0);
L_cen = K_cen*Gx;
G_cen = feedback(G, pinv(Jt')*K_cen*pinv(Ja));

K_svd = diag(1./diag(abs(evalfr(Gsvd, j*wc))))*(1/abs(evalfr(1/(1 + s/w0), j*wc)))/(1 + s/w0);
L_svd = K_svd*Gsvd;
U_inv = inv(U);
G_svd = feedback(G, inv(V')*K_svd*U_inv(1:3, :));



% The obtained diagonal elements of the loop gains are shown in Figure [[fig:gravimeter_comp_loop_gain_diagonal]].


figure;
tiledlayout(3, 1, 'TileSpacing', 'None', 'Padding', 'None');

% Magnitude
ax1 = nexttile([2, 1]);
hold on;
plot(freqs, abs(squeeze(freqresp(L_svd(1, 1), freqs, 'Hz'))), 'DisplayName', '$L_{SVD}(i,i)$');
for i_in_out = 2:3
  set(gca,'ColorOrderIndex',1)
  plot(freqs, abs(squeeze(freqresp(L_svd(i_in_out, i_in_out), freqs, 'Hz'))), 'HandleVisibility', 'off');
end

set(gca,'ColorOrderIndex',2)
plot(freqs, abs(squeeze(freqresp(L_cen(1, 1), freqs, 'Hz'))), ...
     'DisplayName', '$L_{J}(i,i)$');
for i_in_out = 2:3
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
for i_in_out = 1:3
  set(gca,'ColorOrderIndex',1)
  plot(freqs, 180/pi*angle(squeeze(freqresp(L_svd(i_in_out, i_in_out), freqs, 'Hz'))));
end
set(gca,'ColorOrderIndex',2)
for i_in_out = 1:3
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
% <<sec:gravimeter_closed_loop_results>>

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

% The obtained transmissibility in Open-loop, for the centralized control as well as for the SVD control are shown in Figure [[fig:gravimeter_platform_simscape_cl_transmissibility]].


freqs = logspace(-2, 2, 1000);

figure;
tiledlayout(1, 3, 'TileSpacing', 'None', 'Padding', 'None');

ax1 = nexttile;
hold on;
plot(freqs, abs(squeeze(freqresp(G(    1,1)/s^2, freqs, 'Hz'))), 'DisplayName', 'Open-Loop');
plot(freqs, abs(squeeze(freqresp(G_cen(1,1)/s^2, freqs, 'Hz'))), 'DisplayName', 'Centralized');
plot(freqs, abs(squeeze(freqresp(G_svd(1,1)/s^2, freqs, 'Hz'))), '--', 'DisplayName', 'SVD');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Transmissibility'); xlabel('Frequency [Hz]');
title('$D_x/D_{w,x}$');
legend('location', 'southwest');

ax2 = nexttile;
hold on;
plot(freqs, abs(squeeze(freqresp(G(    2,2)/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_cen(2,2)/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_svd(2,2)/s^2, freqs, 'Hz'))), '--');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
set(gca, 'YTickLabel',[]); xlabel('Frequency [Hz]');
title('$D_z/D_{w,z}$');

ax3 = nexttile;
hold on;
plot(freqs, abs(squeeze(freqresp(G(    3,3)/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_cen(3,3)/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_svd(3,3)/s^2, freqs, 'Hz'))), '--');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
set(gca, 'YTickLabel',[]); xlabel('Frequency [Hz]');
title('$R_y/R_{w,y}$');

linkaxes([ax1,ax2,ax3],'xy');
xlim([freqs(1), freqs(end)]);
xlim([1e-2, 5e1]); ylim([1e-7, 1e-2]);
