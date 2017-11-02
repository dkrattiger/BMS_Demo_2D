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

% find path for current file
functionPathName = mfilename('fullpath');
istringcut = find(functionPathName=='/',1,'last');
functionPathName = functionPathName(1:istringcut);

% add current libraries folder and all subfolders to path
addpath(genpath([functionPathName,'libraries']))

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

n_om = 2;
omega_max = 3e4*2*pi;
omega = linspace(0,omega_max,n_om);

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
tstart_full = tic;
C_free = [];
options.n_curves = n_curves;
[kappa_full,~,t_wloop_full] = dispersion_solver_k_w(omega,K_free,C_free,M_free,dof_sets,R,options);

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
options.n_curves = n_curves;

% compute BMS solution
n_om_BMS = 100;
omega2 = linspace(0,omega_max,n_om_BMS);
[kappa_BMS,~,t_wloop_BMS] = dispersion_solver_k_w(omega2,K_BMS,C_BMS,M_BMS,dof_sets_BMS,R,options);

%% plot BMS k(w) solution
% ======================================================================= %
curve_plot = 1:n_curves;
figure(5);clf;hold on;view(3)
h1 = plot3(real(kappa_BMS(curve_plot,:)),imag(kappa_BMS(curve_plot,:)),omega/(2*pi),'g.');hold on
h2 = plot3(real( kappa_full(curve_plot,:)),imag(kappa_full(curve_plot,:)),omega/(2*pi),'ko');hold on
legend([h1(1),h2(1)],'BMS','Full')

options_plot.MarkerOrder = {'.','o'};
options_plot.ThreeD = true;
om_plot = 1:length(omega);
omega_plot = omega(om_plot);
kappas_plot = {kappa_BMS(curve_plot,om_plot),...
               kappa_full(curve_plot,om_plot)};
figure(4);clf
h1 = dispersion_plot_k_w(omega_plot,kappas_plot,options_plot);
h2 = dispersion_plot_k_w(omega_plot,kappas_plot,options_plot);



%% plot dispersion
% ======================================================================= %

figure(2);clf
h1 = plot(kappa_plot,f_full,'k-','linewidth',1.5);hold on
h2 = plot(kappa_plot,f_BMS,'g--','linewidth',1.5);
set(gca,'xtick',linspace(0,kappa_plot(end),length(sym_pts)))
set(gca,'xticklabels',sym_pts)
xlabel('Wave number')
ylabel('Frequency (Hz)')

%% Remove subfolders from Matlab path
% ======================================================================= %
rmpath(genpath([functionPathName,'libraries']))
