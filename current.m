close all; clear all; clc;

load distribution
load_consts

Curren=zeros(1,10);

for k=1:20
vo=history_vo{1,k};
v = 0.4+0.1*k;      % Voltage of anode (higher y) to cathode (y=0) (V)

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

% Calculate jf
% Matrices to be put into jf in two methods
% 1 1 1 
% 2 2 2
% 3 3 3
mat_vert = - R_vo' + R_vo;
% 123
%     123
%         123
jfi = zeros(3*nvo2, 1);
jfj = zeros(3*nvo2, 1);
jfs = zeros(3*nvo2, 1);
jfnum = 0;
mat_diag = - R_vo + R_vo';
for i = 1:nvo    
    m = sparse(1:nvo, 1:nvo, mat_vert(:, i));
    m(i, :) = m(i, :) + mat_diag(i, :);
    
    [nowjfi, nowjfj, nowjfs] = find(m);
    nowjfend = jfnum + length(nowjfi);    
    jfi((jfnum+1):nowjfend) = nowjfi + (i-1)*nvo;
    jfj((jfnum+1):nowjfend) = nowjfj;
    jfs((jfnum+1):nowjfend) = nowjfs;
    jfnum = nowjfend;
end
jfi = jfi(jfi ~= 0);
jfj = jfj(jfi ~= 0);
jfs = jfs(jfi ~= 0);
jf = sparse(jfi, jfj, jfs);

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
