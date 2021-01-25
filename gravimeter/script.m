%% Clear Workspace and Close figures
clear; close all; clc;

%% Intialize Laplace variable
s = zpk('s');

freqs = logspace(-1, 3, 1000);

% Gravimeter Model - Parameters
% <<sec:gravimeter_model>>


open('gravimeter.slx')



% The model of the gravimeter is schematically shown in Figure [[fig:gravimeter_model]].

% #+name: fig:gravimeter_model
% #+caption: Model of the gravimeter
% [[file:figs/gravimeter_model.png]]

% #+name: fig:leg_model
% #+caption: Model of the struts
% [[file:figs/leg_model.png]]

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
G.OutputName = {'Ax1', 'Ay1', 'Ax2', 'Ay2'};



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
            ylabel('Amplitude [$\frac{m/s^2}{N}$]')
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

Ja = [1 0 -h/2
      0 1  l/2
      1 0  h/2
      0 1  0];

Jt = [1 0 -ha
      0 1  la
      0 1 -la];



% And the plant $\bm{G}_x$ is computed:

Gx = pinv(Ja)*G*pinv(Jt');
Gx.InputName  = {'Fx', 'Fy', 'Mz'};
Gx.OutputName  = {'Dx', 'Dy', 'Rz'};

size(Gx)



% #+RESULTS:
% : size(Gx)
% : State-space model with 3 outputs, 3 inputs, and 6 states.

% The diagonal and off-diagonal elements of $G_x$ are shown in Figure [[fig:gravimeter_jacobian_plant]].

% It is shown at the system is:
% - decoupled at high frequency thanks to a diagonal mass matrix (the Jacobian being evaluated at the center of mass of the payload)
% - coupled at low frequency due to the non-diagonal terms in the stiffness matrix, especially the term corresponding to a coupling between a force in the x direction to a rotation around z (due to the torque applied by the stiffness 1).

% The choice of the frame in this the Jacobian is evaluated is discussed in Section [[sec:choice_jacobian_reference]].


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



% #+name: fig:gravimeter_rga
% #+caption: Obtained norm of RGA elements for the SVD decoupled plant and the Jacobian decoupled plant
% #+RESULTS:
% [[file:figs/gravimeter_rga.png]]

% The RGA-number is also a measure of diagonal dominance:
% \begin{equation}
%   \text{RGA-number} = \| \Lambda(G) - I \|_\text{sum}
% \end{equation}


% Relative Gain Array for the decoupled plant using SVD:
RGA_svd = zeros(size(Gsvd,1), size(Gsvd,2), length(freqs));
Gsvd_inv = inv(Gsvd);
for f_i = 1:length(freqs)
    RGA_svd(:, :, f_i) = abs(evalfr(Gsvd, j*2*pi*freqs(f_i)).*evalfr(Gsvd_inv, j*2*pi*freqs(f_i))');
end

% Relative Gain Array for the decoupled plant using the Jacobian:
RGA_x = zeros(size(Gx,1), size(Gx,2), length(freqs));
Gx_inv = inv(Gx);
for f_i = 1:length(freqs)
    RGA_x(:, :, f_i) = abs(evalfr(Gx, j*2*pi*freqs(f_i)).*evalfr(Gx_inv, j*2*pi*freqs(f_i))');
end

RGA_num_svd = squeeze(sum(sum(RGA_svd - eye(3))));
RGA_num_x = squeeze(sum(sum(RGA_x - eye(3))));

figure;
hold on;
plot(freqs, RGA_num_svd)
plot(freqs, RGA_num_x)
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
xlabel('Frequency [Hz]'); ylabel('RGA-Number');

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
plot(freqs, abs(squeeze(freqresp(Gx(3, 3), freqs, 'Hz'))), 'DisplayName', '$G_x(3,3) = R_z/M_z$');
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

K_svd = diag(1./diag(abs(evalfr(Gsvd, j*wc))))*(1/abs(evalfr(1/(1 + s/w0), j*wc)))/(1 + s/w0);
L_svd = K_svd*Gsvd;
U_inv = inv(U);



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

% Now the system is identified again with additional inputs and outputs:
% - $x$, $y$ and $R_z$ ground motion
% - $x$, $y$ and $R_z$ acceleration of the payload.


%% Name of the Simulink File
mdl = 'gravimeter';

%% Input/Output definition
clear io; io_i = 1;
io(io_i) = linio([mdl, '/Dx'], 1, 'openinput');  io_i = io_i + 1;
io(io_i) = linio([mdl, '/Dy'], 1, 'openinput');  io_i = io_i + 1;
io(io_i) = linio([mdl, '/Rz'], 1, 'openinput');  io_i = io_i + 1;
io(io_i) = linio([mdl, '/F1'], 1, 'openinput');  io_i = io_i + 1;
io(io_i) = linio([mdl, '/F2'], 1, 'openinput');  io_i = io_i + 1;
io(io_i) = linio([mdl, '/F3'], 1, 'openinput');  io_i = io_i + 1;
io(io_i) = linio([mdl, '/Abs_Motion'], 1, 'openoutput'); io_i = io_i + 1;
io(io_i) = linio([mdl, '/Abs_Motion'], 2, 'openoutput'); io_i = io_i + 1;
io(io_i) = linio([mdl, '/Abs_Motion'], 3, 'openoutput'); io_i = io_i + 1;
io(io_i) = linio([mdl, '/Acc_side'], 1, 'openoutput'); io_i = io_i + 1;
io(io_i) = linio([mdl, '/Acc_side'], 2, 'openoutput'); io_i = io_i + 1;
io(io_i) = linio([mdl, '/Acc_top'], 1, 'openoutput'); io_i = io_i + 1;
io(io_i) = linio([mdl, '/Acc_top'], 2, 'openoutput'); io_i = io_i + 1;

G = linearize(mdl, io);
G.InputName  = {'Dx', 'Dy', 'Rz', 'F1', 'F2', 'F3'};
G.OutputName = {'Ax', 'Ay', 'Arz', 'Ax1', 'Ay1', 'Ax2', 'Ay2'};



% The loop is closed using the developed controllers.

G_cen = lft(G, -pinv(Jt')*K_cen*pinv(Ja));
G_svd = lft(G, -inv(V')*K_svd*U_inv(1:3, :));



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
plot(freqs, abs(squeeze(freqresp(G(    'Ax','Dx')/s^2, freqs, 'Hz'))), 'DisplayName', 'Open-Loop');
plot(freqs, abs(squeeze(freqresp(G_cen('Ax','Dx')/s^2, freqs, 'Hz'))), 'DisplayName', 'Centralized');
plot(freqs, abs(squeeze(freqresp(G_svd('Ax','Dx')/s^2, freqs, 'Hz'))), '--', 'DisplayName', 'SVD');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Transmissibility'); xlabel('Frequency [Hz]');
title('$D_x/D_{w,x}$');
legend('location', 'southwest');

ax2 = nexttile;
hold on;
plot(freqs, abs(squeeze(freqresp(G(    'Ay','Dy')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_cen('Ay','Dy')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_svd('Ay','Dy')/s^2, freqs, 'Hz'))), '--');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
set(gca, 'YTickLabel',[]); xlabel('Frequency [Hz]');
title('$D_y/D_{w,y}$');

ax3 = nexttile;
hold on;
plot(freqs, abs(squeeze(freqresp(G(    'Arz','Rz')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_cen('Arz','Rz')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_svd('Arz','Rz')/s^2, freqs, 'Hz'))), '--');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
set(gca, 'YTickLabel',[]); xlabel('Frequency [Hz]');
title('$R_z/R_{w,z}$');

linkaxes([ax1,ax2,ax3],'xy');
xlim([freqs(1), freqs(end)]);
xlim([1e-2, 5e1]); ylim([1e-2, 1e1]);



% #+name: fig:gravimeter_platform_simscape_cl_transmissibility
% #+caption: Obtained Transmissibility
% #+RESULTS:
% [[file:figs/gravimeter_platform_simscape_cl_transmissibility.png]]


freqs = logspace(-2, 2, 1000);

figure;
hold on;
for out_i = 1:3
    for in_i = out_i+1:3
        set(gca,'ColorOrderIndex',1)
        plot(freqs, abs(squeeze(freqresp(G(    out_i,in_i), freqs, 'Hz'))));
        set(gca,'ColorOrderIndex',2)
        plot(freqs, abs(squeeze(freqresp(G_cen(out_i,in_i), freqs, 'Hz'))));
        set(gca,'ColorOrderIndex',3)
        plot(freqs, abs(squeeze(freqresp(G_svd(out_i,in_i), freqs, 'Hz'))), '--');
    end
end
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Transmissibility'); xlabel('Frequency [Hz]');
ylim([1e-6, 1e3]);

% Robustness to a change of actuator position
% <<sec:robustness_actuator_position>>

% Let say we change the position of the actuators:

la = l/2*0.7; % Position of Act. [m]
ha = h/2*0.7; % Position of Act. [m]

%% Name of the Simulink File
mdl = 'gravimeter';

%% Input/Output definition
clear io; io_i = 1;
io(io_i) = linio([mdl, '/Dx'], 1, 'openinput');  io_i = io_i + 1;
io(io_i) = linio([mdl, '/Dy'], 1, 'openinput');  io_i = io_i + 1;
io(io_i) = linio([mdl, '/Rz'], 1, 'openinput');  io_i = io_i + 1;
io(io_i) = linio([mdl, '/F1'], 1, 'openinput');  io_i = io_i + 1;
io(io_i) = linio([mdl, '/F2'], 1, 'openinput');  io_i = io_i + 1;
io(io_i) = linio([mdl, '/F3'], 1, 'openinput');  io_i = io_i + 1;
io(io_i) = linio([mdl, '/Abs_Motion'], 1, 'openoutput'); io_i = io_i + 1;
io(io_i) = linio([mdl, '/Abs_Motion'], 2, 'openoutput'); io_i = io_i + 1;
io(io_i) = linio([mdl, '/Abs_Motion'], 3, 'openoutput'); io_i = io_i + 1;
io(io_i) = linio([mdl, '/Acc_side'], 1, 'openoutput'); io_i = io_i + 1;
io(io_i) = linio([mdl, '/Acc_side'], 2, 'openoutput'); io_i = io_i + 1;
io(io_i) = linio([mdl, '/Acc_top'], 1, 'openoutput'); io_i = io_i + 1;
io(io_i) = linio([mdl, '/Acc_top'], 2, 'openoutput'); io_i = io_i + 1;

G = linearize(mdl, io);
G.InputName  = {'Dx', 'Dy', 'Rz', 'F1', 'F2', 'F3'};
G.OutputName = {'Ax', 'Ay', 'Arz', 'Ax1', 'Ay1', 'Ax2', 'Ay2'};



% The loop is closed using the developed controllers.

G_cen_b = lft(G, -pinv(Jt')*K_cen*pinv(Ja));
G_svd_b = lft(G, -inv(V')*K_svd*U_inv(1:3, :));



% The new plant is computed, and the centralized and SVD control architectures are applied using the previously computed Jacobian matrices and $U$ and $V$ matrices.

% The closed-loop system are still stable in both cases, and the obtained transmissibility are equivalent as shown in Figure [[fig:gravimeter_transmissibility_offset_act]].


freqs = logspace(-2, 2, 1000);

figure;
tiledlayout(1, 3, 'TileSpacing', 'None', 'Padding', 'None');

ax1 = nexttile;
hold on;
plot(freqs, abs(squeeze(freqresp(G_cen(      'Ax','Dx')/s^2, freqs, 'Hz'))), 'DisplayName', 'Open-Loop');
plot(freqs, abs(squeeze(freqresp(G_cen_b('Ax','Dx')/s^2, freqs, 'Hz'))), 'DisplayName', 'Centralized');
plot(freqs, abs(squeeze(freqresp(G_svd_b('Ax','Dx')/s^2, freqs, 'Hz'))), '--', 'DisplayName', 'SVD');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
ylabel('Transmissibility'); xlabel('Frequency [Hz]');
title('$D_x/D_{w,x}$');
legend('location', 'southwest');

ax2 = nexttile;
hold on;
plot(freqs, abs(squeeze(freqresp(G_cen(      'Ay','Dy')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_cen_b('Ay','Dy')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_svd_b('Ay','Dy')/s^2, freqs, 'Hz'))), '--');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
set(gca, 'YTickLabel',[]); xlabel('Frequency [Hz]');
title('$D_y/D_{w,y}$');

ax3 = nexttile;
hold on;
plot(freqs, abs(squeeze(freqresp(G_cen(      'Arz','Rz')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_cen_b('Arz','Rz')/s^2, freqs, 'Hz'))));
plot(freqs, abs(squeeze(freqresp(G_svd_b('Arz','Rz')/s^2, freqs, 'Hz'))), '--');
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
set(gca, 'YTickLabel',[]); xlabel('Frequency [Hz]');
title('$R_z/R_{w,z}$');

linkaxes([ax1,ax2,ax3],'xy');
xlim([freqs(1), freqs(end)]);
xlim([1e-2, 5e1]); ylim([1e-2, 1e1]);

% Decoupling of the mass matrix

% #+name: fig:gravimeter_model_M
% #+caption: Choice of {O} such that the Mass Matrix is Diagonal
% [[file:figs/gravimeter_model_M.png]]


la = l/2; % Position of Act. [m]
ha = h/2; % Position of Act. [m]

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
G.OutputName = {'Ax1', 'Ay1', 'Ax2', 'Ay2'};



% Decoupling at the CoM (Mass decoupled)

JMa = [1 0 -h/2
       0 1  l/2
       1 0  h/2
       0 1  0];

JMt = [1 0 -ha
       0 1  la
       0 1 -la];

GM = pinv(JMa)*G*pinv(JMt');
GM.InputName  = {'Fx', 'Fy', 'Mz'};
GM.OutputName  = {'Dx', 'Dy', 'Rz'};

figure;

% Magnitude
hold on;
for i_in = 1:3
    for i_out = [1:i_in-1, i_in+1:3]
        plot(freqs, abs(squeeze(freqresp(GM(i_out, i_in), freqs, 'Hz'))), 'color', [0,0,0,0.2], ...
             'HandleVisibility', 'off');
    end
end
plot(freqs, abs(squeeze(freqresp(GM(i_out, i_in), freqs, 'Hz'))), 'color', [0,0,0,0.2], ...
     'DisplayName', '$G_x(i,j)\ i \neq j$');
set(gca,'ColorOrderIndex',1)
for i_in_out = 1:3
    plot(freqs, abs(squeeze(freqresp(GM(i_in_out, i_in_out), freqs, 'Hz'))), 'DisplayName', sprintf('$G_x(%d,%d)$', i_in_out, i_in_out));
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
xlabel('Frequency [Hz]'); ylabel('Magnitude');
legend('location', 'southeast');
ylim([1e-8, 1e0]);

% Decoupling of the stiffness matrix

% #+name: fig:gravimeter_model_K
% #+caption: Choice of {O} such that the Stiffness Matrix is Diagonal
% [[file:figs/gravimeter_model_K.png]]

% Decoupling at the point where K is diagonal (x = 0, y = -h/2 from the schematic {O} frame):

JKa = [1 0  0
       0 1 -l/2
       1 0 -h
       0 1  0];

JKt = [1 0  0
       0 1 -la
       0 1  la];



% And the plant $\bm{G}_x$ is computed:

GK = pinv(JKa)*G*pinv(JKt');
GK.InputName  = {'Fx', 'Fy', 'Mz'};
GK.OutputName  = {'Dx', 'Dy', 'Rz'};

figure;

% Magnitude
hold on;
for i_in = 1:3
    for i_out = [1:i_in-1, i_in+1:3]
        plot(freqs, abs(squeeze(freqresp(GK(i_out, i_in), freqs, 'Hz'))), 'color', [0,0,0,0.2], ...
             'HandleVisibility', 'off');
    end
end
plot(freqs, abs(squeeze(freqresp(GK(i_out, i_in), freqs, 'Hz'))), 'color', [0,0,0,0.2], ...
     'DisplayName', '$G_x(i,j)\ i \neq j$');
set(gca,'ColorOrderIndex',1)
for i_in_out = 1:3
    plot(freqs, abs(squeeze(freqresp(GK(i_in_out, i_in_out), freqs, 'Hz'))), 'DisplayName', sprintf('$G_x(%d,%d)$', i_in_out, i_in_out));
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
xlabel('Frequency [Hz]'); ylabel('Magnitude');
legend('location', 'southeast');
ylim([1e-8, 1e0]);

% Combined decoupling of the mass and stiffness matrices

% #+name: fig:gravimeter_model_KM
% #+caption: Ideal location of the actuators such that both the mass and stiffness matrices are diagonal
% [[file:figs/gravimeter_model_KM.png]]

% To do so, the actuator position should be modified


la = l/2; % Position of Act. [m]
ha = 0; % Position of Act. [m]

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
G.OutputName = {'Ax1', 'Ay1', 'Ax2', 'Ay2'};

JMa = [1 0 -h/2
       0 1  l/2
       1 0  h/2
       0 1  0];

JMt = [1 0 -ha
       0 1  la
       0 1 -la];

GKM = pinv(JMa)*G*pinv(JMt');
GKM.InputName  = {'Fx', 'Fy', 'Mz'};
GKM.OutputName  = {'Dx', 'Dy', 'Rz'};

figure;

% Magnitude
hold on;
for i_in = 1:3
    for i_out = [1:i_in-1, i_in+1:3]
        plot(freqs, abs(squeeze(freqresp(GKM(i_out, i_in), freqs, 'Hz'))), 'color', [0,0,0,0.2], ...
             'HandleVisibility', 'off');
    end
end
plot(freqs, abs(squeeze(freqresp(GKM(i_out, i_in), freqs, 'Hz'))), 'color', [0,0,0,0.2], ...
     'DisplayName', '$G_x(i,j)\ i \neq j$');
set(gca,'ColorOrderIndex',1)
for i_in_out = 1:3
    plot(freqs, abs(squeeze(freqresp(GKM(i_in_out, i_in_out), freqs, 'Hz'))), 'DisplayName', sprintf('$G_x(%d,%d)$', i_in_out, i_in_out));
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
xlabel('Frequency [Hz]'); ylabel('Magnitude');
legend('location', 'southeast');
ylim([1e-8, 1e0]);

% SVD decoupling performances
% <<sec:decoupling_performances>>
% As the SVD is applied on a *real approximation* of the plant dynamics at a frequency $\omega_0$, it is foreseen that the effectiveness of the decoupling depends on the validity of the real approximation.

% Let's do the SVD decoupling on a plant that is mostly real (low damping) and one with a large imaginary part (larger damping).

% Start with small damping, the obtained diagonal and off-diagonal terms are shown in Figure [[fig:gravimeter_svd_low_damping]].

c = 2e1; % Actuator Damping [N/(m/s)]

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
G.OutputName = {'Ax1', 'Ay1', 'Ax2', 'Ay2'};

wc = 2*pi*10; % Decoupling frequency [rad/s]
H1 = evalfr(G, j*wc);
D = pinv(real(H1'*H1));
H1 = pinv(D*real(H1'*diag(exp(j*angle(diag(H1*D*H1.'))/2))));
[U,S,V] = svd(H1);
Gsvd = inv(U)*G*inv(V');

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
     'DisplayName', '$G_{svd}(i,j)\ i \neq j$');
set(gca,'ColorOrderIndex',1)
for i_in_out = 1:3
    plot(freqs, abs(squeeze(freqresp(Gsvd(i_in_out, i_in_out), freqs, 'Hz'))), 'DisplayName', sprintf('$G_{svd}(%d,%d)$', i_in_out, i_in_out));
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
xlabel('Frequency [Hz]'); ylabel('Magnitude');
legend('location', 'northwest');
ylim([1e-8, 1e0]);



% #+name: fig:gravimeter_svd_low_damping
% #+caption: Diagonal and off-diagonal term when decoupling with SVD on the gravimeter with small damping
% #+RESULTS:
% [[file:figs/gravimeter_svd_low_damping.png]]

% Now take a larger damping, the obtained diagonal and off-diagonal terms are shown in Figure [[fig:gravimeter_svd_high_damping]].

c = 5e2; % Actuator Damping [N/(m/s)]

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
G.OutputName = {'Ax1', 'Ay1', 'Ax2', 'Ay2'};

wc = 2*pi*10; % Decoupling frequency [rad/s]
H1 = evalfr(G, j*wc);
D = pinv(real(H1'*H1));
H1 = pinv(D*real(H1'*diag(exp(j*angle(diag(H1*D*H1.'))/2))));
[U,S,V] = svd(H1);
Gsvdd = inv(U)*G*inv(V');

figure;

% Magnitude
hold on;
for i_in = 1:3
    for i_out = [1:i_in-1, i_in+1:3]
        plot(freqs, abs(squeeze(freqresp(Gsvdd(i_out, i_in), freqs, 'Hz'))), 'color', [0,0,0,0.2], ...
             'HandleVisibility', 'off');
    end
end
plot(freqs, abs(squeeze(freqresp(Gsvdd(i_out, i_in), freqs, 'Hz'))), 'color', [0,0,0,0.2], ...
     'DisplayName', '$G_{svd}(i,j)\ i \neq j$');
set(gca,'ColorOrderIndex',1)
for i_in_out = 1:3
    plot(freqs, abs(squeeze(freqresp(Gsvdd(i_in_out, i_in_out), freqs, 'Hz'))), 'DisplayName', sprintf('$G_{svd}(%d,%d)$', i_in_out, i_in_out));
end
hold off;
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
xlabel('Frequency [Hz]'); ylabel('Magnitude');
legend('location', 'northwest');
ylim([1e-8, 1e0]);
