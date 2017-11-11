%% BMS Demo: 2D Plane Strain
% ======================================================================= %
% This code computes the k(w) band structure for a 2D model of a square
% unit cell with 23040 DOFs (after enforcing periodicity). It does so using
% the full model (Warning this calculations is time consuming) and with the
% BMS model.

% Dimitri Krattiger
% Oct. 4 2017

clear;
clc;
close all

%% Add subfolders with dependent libraries to Matlab path
% ======================================================================= %

% add current libraries folder and all subfolders to path
addpath(genpath('libraries'))

%% Solution Options
% ======================================================================= %

% number of dispersion branches to compute
n_curves = 8;

%% Circular Inclusion model dimensions
% ======================================================================= %

% unit cell dimensions
Lx = 0.05;
Ly = 0.05;

% lattice vectors
r1 = [Lx;0];
r2 = [0;Ly];
R = [r1,r2];

%% Create frequency vector
% ======================================================================= %

omega_max = 10e4*2*pi;

% note, the full solution takes a long time (even with dynamic reduction).
n_om = 21;
omega_full = linspace(0,omega_max,n_om);

n_om_BMS = 101;
omega_BMS = linspace(0,omega_max,n_om_BMS);

%% Load Model
% ======================================================================= %

% Load model (mass and stiffness matrices and coordinates)
model_savestring = '8LobeInclusion_23040DOF_2ndOrderEles_FreeModel';
load(model_savestring);

% Load mesh plotting details (only necessary if unit-cell geometry is plotted)
model_plot_savestring = '8LobeInclusion_23040DOF_2ndOrderEles_FreeModel_Plot';
load(model_plot_savestring)

%% Display unit cell geometry
% ======================================================================= %

figure(1);clf
h_patch = plot_FEM_model(coordinates,patchfaces,C,fedges);
axis equal
title('Unit Cell')
axis off
drawnow

%% k(w) Dispersion Solution
% ======================================================================= %

% compute node sets for overall unit cell
[node_sets] = find_node_sets(coordinates,R);
dof_sets = node2dof(node_sets,2);

% Full FE model w(k) Dispersion Solution
% Note this takes a long time (about 200s per freq. pt on my laptop)
tstart_full = tic;
C_free = [];
options.n_curves = n_curves;
% options.full_eig           = false;
% options.dynamicReduction   = false;

profile clear
profile on
[kappa_full,~,t_wloop_full] = dispersion_solver_k_w(omega_full,K_free,C_free,M_free,dof_sets,R,options);

profile viewer

t_full = toc(tstart_full);

%% Perform BMS Reduction
% =================================================================== %

% BMS reduction parameters (stored in options structure)
options_BMS.InteriorMethod        = 'CB+';
options_BMS.BoundaryMethod        = 'exact';
options_BMS.n_FI                  = 30;
options_BMS.n_CC                  = 12;

% perform BMS reduction 
[K_BMS,M_BMS,dof_sets_BMS,info_BMS,T_BMS] = BMS(K_free,M_free,coordinates,R,options_BMS);

% solve for BMS dispersion
C_BMS = [];
clear options
options.n_curves = n_curves;

% compute BMS solution
[kappa_BMS,~,t_wloop_BMS] = dispersion_solver_k_w(omega_BMS,K_BMS,C_BMS,M_BMS,dof_sets_BMS,R,options);

%% plot BMS k(w) solution
% ======================================================================= %
curve_plot = 1:n_curves;

figure(5);clf;hold on;view(3)
h1 = plot3(real(kappa_BMS(curve_plot,:)),imag(kappa_BMS(curve_plot,:)),omega_BMS/(2*pi),'g.');hold on
h2 = plot3(real( kappa_full(curve_plot,:)),imag(kappa_full(curve_plot,:)),omega_full/(2*pi),'ko');hold on
legend([h1(1),h2(1)],'BMS','Full')

%% plot dispersion
% ======================================================================= %

figure(4);clf

% specify plot options (case sensitive)
options_plot.Markers = {'o','.'};
options_plot.LineStyles = {'none','-'};
options_plot.ThreeD = true;
options_plot.legendstrings = {'full','BMS'};

% list out plot vectors
omegas = {omega_full,omega_BMS};
kappas = {kappa_full(curve_plot,:),kappa_BMS(curve_plot,:)};

% call plotting routine
h = dispersion_plot_k_w(omegas,kappas,options_plot);

%% Remove subfolders from Matlab path
% ======================================================================= %
% rmpath(genpath('libraries'))