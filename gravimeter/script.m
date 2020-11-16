%% Clear Workspace and Close figures
clear; close all; clc;

%% Intialize Laplace variable
s = zpk('s');

% Simscape Model - Parameters

open('gravimeter.slx')



% Parameters

l  = 1.0; % Length of the mass [m]
h  = 1.7; % Height of the mass [m]

la = l/2; % Position of Act. [m]
ha = h/2; % Position of Act. [m]

m = 400; % Mass [kg]
I = 115; % Inertia [kg m^2]

k = 15e3; % Actuator Stiffness [N/m]
c = 0.03; % Actuator Damping [N/(m/s)]

deq = 0.2; % Length of the actuators [m]

g = 0; % Gravity [m/s2]

% System Identification - Without Gravity

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

% \begin{equation}
%   \bm{a} = \begin{bmatrix} a_{1x} \\ a_{1z} \\ a_{2x} \\ a_{2z} \end{bmatrix}
% \end{equation}

% \begin{equation}
%   \bm{\tau} = \begin{bmatrix}\tau_1 \\ \tau_2 \\ \tau_2 \end{bmatrix}
% \end{equation}

% We can check the poles of the plant:


pole(G)



% #+RESULTS:
% #+begin_example
%       -0.000183495485977108 +       13.546056874877i
%       -0.000183495485977108 -       13.546056874877i
%       -7.49842878906757e-05 +      8.65934902322567i
%       -7.49842878906757e-05 -      8.65934902322567i
%       -1.33171230256362e-05 +      3.64924169037897i
%       -1.33171230256362e-05 -      3.64924169037897i
% #+end_example

% The plant as 6 states as expected (2 translations + 1 rotation)

size(G)



% #+RESULTS:
% : State-space model with 4 outputs, 3 inputs, and 6 states.

% The bode plot of all elements of the plant are shown in Figure [[fig:open_loop_tf]].


freqs = logspace(-1, 2, 1000);

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

% The jacobian corresponding to the sensors and actuators are defined below.

Ja = [1 0  h/2
      0 1 -l/2
      1 0 -h/2
      0 1  0];

Jt = [1 0  ha
      0 1 -la
      0 1  la];

Gx = pinv(Ja)*G*pinv(Jt');
Gx.InputName  = {'Fx', 'Fz', 'My'};
Gx.OutputName  = {'Dx', 'Dz', 'Ry'};



% The diagonal and off-diagonal elements of $G_x$ are shown in Figure [[fig:gravimeter_jacobian_plant]].


freqs = logspace(-1, 2, 1000);

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

% Real Approximation of $G$ at the decoupling frequency
% <<sec:gravimeter_real_approx>>

% Let's compute a real approximation of the complex matrix $H_1$ which corresponds to the the transfer function $G_u(j\omega_c)$ from forces applied by the actuators to the measured acceleration of the top platform evaluated at the frequency $\omega_c$.

wc = 2*pi*10; % Decoupling frequency [rad/s]

H1 = evalfr(G, j*wc);



% The real approximation is computed as follows:

D = pinv(real(H1'*H1));
H1 = inv(D*real(H1'*diag(exp(j*angle(diag(H1*D*H1.'))/2))));

% SVD Decoupling
% <<sec:gravimeter_svd_decoupling>>

% First, the Singular Value Decomposition of $H_1$ is performed:
% \[ H_1 = U \Sigma V^H \]


[U,~,V] = svd(H1);



% #+name: fig:gravimeter_decouple_svd
% #+caption: Decoupled plant $\bm{G}_{SVD}$ using the Singular Value Decomposition
% #+RESULTS:
% [[file:figs/gravimeter_decouple_svd.png]]

% The decoupled plant is then:
% \[ G_{SVD}(s) = U^{-1} G_u(s) V^{-H} \]


Gsvd = inv(U)*G*inv(V');



% The diagonal and off-diagonal elements of the "SVD" plant are shown in Figure [[fig:gravimeter_svd_plant]].


freqs = logspace(-1, 2, 1000);

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
legend('location', 'southeast', 'FontSize', 8);
ylim([1e-8, 1e0]);

% TODO Verification of the decoupling using the "Gershgorin Radii"
% <<sec:comp_decoupling>>

% The "Gershgorin Radii" is computed for the coupled plant $G(s)$, for the "Jacobian plant" $G_x(s)$ and the "SVD Decoupled Plant" $G_{SVD}(s)$:

% The "Gershgorin Radii" of a matrix $S$ is defined by:
% \[ \zeta_i(j\omega) = \frac{\sum\limits_{j\neq i}|S_{ij}(j\omega)|}{|S_{ii}(j\omega)|} \]

% This is computed over the following frequencies.

freqs = logspace(-2, 2, 1000); % [Hz]

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
plot(freqs, 0.5*ones(size(freqs)), 'k--', 'DisplayName', 'Limit')
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
hold off;
xlabel('Frequency (Hz)'); ylabel('Gershgorin Radii')
legend('location', 'northwest');
ylim([1e-3, 1e3]);

% TODO Obtained Decoupled Plants
% <<sec:gravimeter_decoupled_plant>>

% The bode plot of the diagonal and off-diagonal elements of $G_{SVD}$ are shown in Figure [[fig:simscape_model_decoupled_plant_svd]].


freqs = logspace(-1, 2, 1000);

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


freqs = logspace(-1, 2, 1000);

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



% The obtained diagonal elements of the loop gains are shown in Figure [[fig:gravimeter_comp_loop_gain_diagonal]].


freqs = logspace(-1, 2, 1000);

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

% TODO Closed-Loop system Performances
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
