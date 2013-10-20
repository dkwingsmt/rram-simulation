%%--------------- Physical constants ---------------
kb = 1.38e-23;  % Boltzmann constant (SI)
%mp = 1.67e-27;  % Proton mass (SI)
me = 9.11e-31;  % Electron mass (SI)
qe = 1.60e-19;  % Electron charge (SI)
hbar = 1.055e-34;   % Reduced Plank constant (SI)

%%--------------- Experiment constants ---------------
n = 41;     % Number of layers of atoms
m = 41;     % 10nm long
h = 0.25e-9;% (m)
T = 300;    % Global temperature (K)

%%--------------- Model constants ---------------
Ea  = 0.66*qe;   % Average migration barrier height of oxygen ions in HfOx (J)
t0 = 10^-13;    % 1/t0 -- Characteristic vibration frequency of oxygen ions (s)
beta = 100;     % Regeneration probability factor (now assumed constant) (1)
Etrap_e = 1.83*qe;  % Trap energy of an empty Vo below the conduction band of HfOx (J)
Etrap_f = 1.97*qe;  % Trap energy of an filled Vo ... (J)
Ef_ud = -1.9*qe;    % Fermi level of the up/downside electrode (to Ec of HfOx) (J)
R0_hop = 1e12;  % Vibration frequency of an electron (Hz)
a0 = 0.33e-9;   % Attenuation length of the electron wave function in a trap (m)
% R_{tunnel}^0 N^{L,R}, the electronic coupling factor between the electrode
% and the dielectric layer (Hz)
k_tunnel_ud = 1e14; 

me_eff = 0.1 * me;  % Electron effective mass for HfO2

gamma=1e-9*qe;  %m
miu=2.5e-9;   %m^2/(V*s)
ratio=100;  %氧空位与绝缘层的电阻率之�?