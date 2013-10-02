close all; clear all; clc;
load distribution

Curren=zeros(1,10);

for k=1:10
vo=history_vo{1,k};

%%--------------- Physical constants ---------------
kb = 1.38e-23;  % Boltzmann constant (SI)
%mp = 1.67e-27;  % Proton mass (SI)
me = 9.11e-31;  % Electron mass (SI)
qe = 1.60e-19;  % Electron charge (SI)
hbar = 1.055e-34;   % Reduced Plank constant (SI)

%%--------------- Model constants ---------------
Ea  = 1.0*qe;   % Average migration barrier height of oxygen ions in HfOx (J)
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

%%--------------- Experiment constants ---------------
n = 41;     % Number of layers of atoms
m = 101;    % 50nm long (??)
h = 0.25e-9;% (m)
v = 1+0.1*k;      % Voltage of anode (higher y) to cathode (y=0) (V)
T = 300;    % Global temperature (K)
deltat = 0.00005;   % Initial step of time (s)
time = 0.008;      % Total experiment time (s)

nh = n * h;
mh = m * h;
kbT = kb * T;
y_atoms = (1:n - 0.5) * h;  % Y-position of every layer of atoms
VH = y_atoms / (n * h) * v; % Homogenrous component of electrical potential (V)

%%--------------- Calculate current ---------------

%% Calculate R_vo, the hopping rate between two cavancies

% Where vacancies at, in the unit of index of the matrix vo
[vo_i, vo_j] = find(vo);    % column vectors; vo([y, x] in [vo_y, vo_x]) = 1
% Where vacancies at, the absolute position (nm)
vo_y = vo_i * h;
vo_x = vo_j * h;
vo_ymat = repmat(vo_y, [1, length(vo_y)]);
vo_xmat = repmat(vo_x, [1, length(vo_x)]);

% Distances between each two vacancies (nm)
vo_distance = sqrt((vo_ymat-vo_ymat').^2 + (vo_xmat-vo_xmat').^2);
% Voltage between vacacies, vo_VH(m, n) = VH(m) - VH(n) (V)
vo_VH = (vo_ymat - vo_ymat') * v / nh;
% Hopping rate between two vacancies
R_vo = R0_hop * exp(-vo_distance/a0 + qe*vo_VH'/kbT);
% You can not hop to yourself, set diag(R_vo) = 0
R_vo = R_vo - diag(diag(R_vo));


%% Calculate Rin/out_u/d, a.k.a R^{i/o L/R}_n in the paper
% the hopping rate between a vacancy and an electrode

% Number of filled states in the electrode above Ev+
% For given E1, E2, the Fermi-Dirac integral
% F_above(E1, E2) = \int_{E1}^{+\infty} 1/(1+exp( (E-E2)/kT )dE
%  = kT * \int_0^{+\infty} 1/(1+exp( t-(E2-E1)/kT )dt
%  = kT * ln( 1 + exp((E2-E1)/kT) )
% F_below(E1, E2) = \int_{-\infty}^{E1} 1/(1+exp( (E2-E)/kT )dE
%  = \int_{-E1}^{+\infty} 1/(1+exp( (E+E2)/kT )dE
%  = kT * ln( 1 + exp(-(E2-E1)/kT) )
% Here VH is used for V. Note that upside is anode and downside is cathode
Fin_u = kbT * log(1 + exp((Ef_ud+Etrap_e + qe*(vo_VH-v)) /kbT));
Fin_d = kbT * log(1 + exp((Ef_ud+Etrap_e + qe*(vo_VH-0)) /kbT));
Fout_u = kbT * log(1 + exp(-(Ef_ud+Etrap_f + qe*(vo_VH-v)) /kbT));
Fout_d = kbT * log(1 + exp(-(Ef_ud+Etrap_f + qe*(vo_VH-0)) /kbT));

% Tunneling probability from electrode into the trap, note VH(x) = x*v/nh
% Assume T = exp(\int k_ sqrt(b_ - y) dy)
k_ = -2 / hbar * sqrt(2 * me_eff * (qe * v / nh));
b_ = Etrap_e / (qe * v / nh) + vo_y;
Tin_d_e = exp( k_ * (2/3) * ((b_ - vo_y).^(3/2) - (b_ - 0).^(3/2)) );
Tin_u_e = exp( k_ * (2/3) * ((b_ - nh).^(3/2) - (b_ - vo_y).^(3/2)) );
b_ = Etrap_f / (qe * v / nh) + vo_y;
Tout_d_f = exp( k_ * (2/3) * ((b_ - vo_y).^(3/2) - (b_ - 0).^(3/2)) );
Tout_u_f = exp( k_ * (2/3) * ((b_ - nh).^(3/2) - (b_ - vo_y).^(3/2)) );

% Hopping rates between a trap and an electrode
Rin_d = k_tunnel_ud * Fin_d * Tin_d_e;
Rin_u = k_tunnel_ud * Fin_u * Tin_u_e;
Rout_d = k_tunnel_ud * Fout_d * Tout_d_f;
Rout_u = k_tunnel_ud * Fout_u * Tout_u_f;

%R_vo
%Rin_d
%Rin_u
%Rout_d
%Rout_u

func = @(f) (1-f).*(R_vo'*f) - f.*(R_vo*(1-f)) + (Rin_u+Rin_d).*(1-f) - (Rout_d+Rout_u).*f;
%f1 = fsolve(@(f) (1-f).*(R_vo'*f) - f.*(R_vo*(1-f)) + (Rin_u+Rin_d).*(1-f) - (Rout_d+Rout_u).*f, ones(length(vo_i), 1)*0.5);
%Il = -qe * ((1-f1)' * Rin_u - f1' * Rout_u)

%% === Solve f in the equation in order to solve current
% First calculate the Jacobian matrix of the function
% which is in the form of j = jc + jf * f
% where jc is the constant part, a [nvo, nvo] matrix (nvo is amount of vo)
% and jf is [jf1, jf2, ... jf_nvo], a group of matrix of length numvo
% Note that during the calculation, for convenience,
% jf will be stretched to [nvo^2, nvo]

% func = (1-f).*(R_vo'*f) - f.*(R_vo*(1-f)) + (Rin_u+Rin_d).*(1-f) - (Rout_d+Rout_u).*f

nvo = length(R_vo);
nvo2 = nvo^2;   

%# The following equations used Eistein summation notation
%# (1-f).*(R_vo'*f) part
% jc1 = R_vo'
% jf1_ij f_j = - R_ji f_i - delta_ij (R_kireshape(1:9, 3, 3); f_k)
%# - f.*(R_vo*(1-f)) part
% jc2_ii = - R_ij    // Sum j
% jf2_ij f_j = R_ij f_i + delta_ij (R_ik f_k)
%# (Rin_u+Rin_d).*(1-f) - (Rout_d+Rout_u).*f part
% jc3 = diag(- Rin_u - Rin_d - Rout_d - Rout_u)


jc = R_vo' + diag(- sum(R_vo, 2) - Rin_u - Rin_d - Rout_d - Rout_u);
jf = sparse(nvo2, nvo);
%jf = sym(zeros(nvo2, nvo));
% Matrices to be put into jf in two methods
% 1 1 1 
% 2 2 2
% 3 3 3
mat_vert = - R_vo' + R_vo;
% 123
%     123
%         123
mat_diag = - R_vo + R_vo';
for i = 1:nvo
    m = diag(mat_vert(:, i));
    m(i, :) = m(i, :) + mat_diag(i, :);
    jf(((i-1)*nvo+1) : (i*nvo), :) = m;
end

%% === Now solve ===
f0 = ones(nvo, 1) * 0.5;

fs = zeros(nvo, 5);
resids = zeros(nvo, 5);

for loop = 1:5
    ftry = (1-f0).*(R_vo'*f0)-f0.*(R_vo*(1-f0))+(Rin_u+Rin_d).*(1-f0)-(Rout_d+Rout_u).*f0;
    jac = jc + reshape(jf * f0, [nvo, nvo]);
    resid = double(jac\ftry);
    f0 = max(0, min(1, f0 - resid));

    fs(:, loop) = f0;
    resids(:, loop) = resid;
end
fs;
resids;

%% The current
I = -qe * (Rin_u'*(1-f0) - Rout_u'*f0);
Curren(1,k)=I;
end