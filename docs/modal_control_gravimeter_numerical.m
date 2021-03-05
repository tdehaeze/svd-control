clc
clear all
close all

%% System properties
g = 100000;
w0 = 2*pi*.5; % MinusK BM1 tablle
l = 0.5; %[m]
la = 1; %[m]
h = 1.7; %[m]
ha = 1.7;% %[m]
m = 400; %[kg]
k = 15e3;%[N/m]
kv = k;
kh = 15e3;
I = 115;%[kg m^2]
dampv = 0.03;
damph = 0.03;
s = tf('s');

%% State-space model
M = [m 0 0
     0 m 0
     0 0 I];

la = l;
ha = h;
kv = k;
kh = k;
 
%Jacobian of the bottom sensor
Js1 = [1 0 h/2
    0 1 -l/2];

%Jacobian of the top sensor
Js2 = [1 0 -h/2
    0 1 0];

%Jacobian of the actuators
Ja = [1 0 ha/2 %Left horizontal actuator
    %1 0 h/2 %Right horizontal actuator
    0 1 -la/2 %Left vertical actuator
    0 1 la/2]; %Right vertical actuator
Jah = [1 0 ha/2];
Jav = [0 1 -la/2 %Left vertical actuator
    0 1 la/2]; %Right vertical actuator
Jta = Ja';
Jtah = Jah';
Jtav = Jav';
K = kv*Jtav*Jav + kh*Jtah*Jah;    
C = dampv*kv*Jtav*Jav+damph*kh*Jtah*Jah;  

E = [1 0 0
    0 1 0
    0 0 1]; %projecting ground motion in the directions of the legs

AA = [zeros(3) eye(3)
     -M\K -M\C];
 
BB = [zeros(3,3)
    M\Jta ];
  
% CC = [[Js1;Js2] zeros(4,3)];
CC = [[Jah;Jav] zeros(3,3)];

% DD = zeros(4,3);
DD = zeros(3);

G = ss(AA,BB,CC,DD);
%% Modal coordinates
[V,D] = eig(M\K);
Mm = V'*M*V; % Modal mass matrix
Dm = V'*C*V; % Modal damping matrix
Km = V'*K*V; % Modal stiffness matrix

Bm = inv(Mm)*V'*Jta;
% Cm = [Js1;Js2]*V;
Cm = [Jah;Jav]*V;


omega = real(sqrt(inv(Mm)*Km));
zeta = real(0.5*inv(Mm)*Dm*inv(omega));

Gm = [1/(s^2+2*zeta(1,1)*omega(1,1)*s+omega(1,1)^2),0,0;
    0,1/(s^2+2*zeta(2,2)*omega(2,2)*s+omega(2,2)^2),0;
    0,0,1/(s^2+2*zeta(3,3)*omega(3,3)*s+omega(3,3)^2)];
figure(1)
bode(G,Cm*Gm*Bm)
figure(2)
bode(G,Gm)

%% Controller
s = tf('s');
w0 = 2*pi*0.1;
Kc = 100/(1+s/w0);
Knet = inv(Bm)*Kc*inv(Cm);
Gc = -lft(G,Knet); 
isstable(Gc)





