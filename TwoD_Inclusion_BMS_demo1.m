%% BMS Demo: 2D Plane Strain
% ======================================================================= %
% This code computes the full band structure for a 2D model of a square
% unit cell with 23040 DOFs (after enforcing periodicity). It also computes
% the BMS band structure and plots a band-structure comparison as well as a
% mode shape comparison

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
n_curves = 10;

%% Circular Inclusion model dimensions
% ======================================================================= %

% unit cell dimensions
Lx = 0.05;
Ly = 0.05;

% lattice vectors
r1 = [Lx;0];
r2 = [0;Ly];
R = [r1,r2];

%% Create Wave vector
% ======================================================================= %

sym_pts = {'\Gamma','X','M','\Gamma'};
kap_ref = 3; % wave vector refinement parameter (increments will double the refinement)
[kappa,kappa_plot] = wave_vector(sym_pts,kap_ref,R);
n_kap = size(kappa,2);

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


%% w(k) Dispersion Solution
% ======================================================================= %

% compute node sets for overall unit cell
[node_sets] = find_node_sets(coordinates,R);
dof_sets = node2dof(node_sets,2);

% Full FE model w(k) Dispersion Solution
tstart_full = tic;
[omega_full,PHI_full,t_kloop_full] = dispersion_solver_w_k(kappa,K_free,M_free,dof_sets,R,n_curves);
t_full = toc(tstart_full);

% Convert to Hz
f_full = omega_full/(2*pi);
        
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
[w_BMS,PHI_BMS,t_kloop_BMS] = dispersion_solver_w_k(kappa,K_BMS,M_BMS,dof_sets_BMS,R,n_curves);

% Convert to Hz
f_BMS = w_BMS/(2*pi);

% evaluate and display frequency errors 
% (ignore zero-frequencies in error calcuation)
tol = 1e-5;
e_freq = 100*abs(f_BMS-f_full)./f_full;
e_freq_max = max(max(e_freq(f_full>tol*max(max(f_full)))));

% BMS timing
t_BMS = sum(t_kloop_BMS) + info_BMS.t_up_front;

% print performance summary
fprintf('\nMaximum frequency error: %5.2e%%',e_freq_max);
fprintf('\nBMS model band-structure calculation time: %5.2fs',t_BMS);
fprintf('\nfull model band-structure calculation time: %5.2fs\n',t_full);
% disp(['BMS band-structure calculation time: ',num2str(t_BMS),'s']);

%% plot dispersion
% ======================================================================= %

figure(2);clf
h1 = plot(kappa_plot,f_full,'k-','linewidth',1.5);hold on
h2 = plot(kappa_plot,f_BMS,'g--','linewidth',1.5);
set(gca,'xtick',linspace(0,kappa_plot(end),length(sym_pts)))
set(gca,'xticklabels',sym_pts)
xlabel('Wave number')
ylabel('Frequency (Hz)')

drawnow


%% plot full model mode shape
% ======================================================================= %

% form periodicity transformation matrices
T_per = Periodic_Boundary_Conditions(dof_sets);

% choose mode to plot
i_phi_plot = 7;
k_sel = 9;

% mode shape magnitude
mag_PHI = 0.07;

% get full mode shape
PHI_plot_full = PHI_full(:,i_phi_plot,k_sel);

% normalize full mode shape
[~,i_max] = max(imag(PHI_plot_full));
PHI_plot_full = PHI_plot_full/PHI_plot_full(i_max);
PHI_plot_full = PHI_plot_full/max(abs(PHI_plot_full));
PHI_plot_full = (Lx*mag_PHI)*PHI_plot_full;

% compute periodicity transformation for selected k point
kvec = kappa(:,k_sel);
lam  = exp(-1i*R*kvec);    
lam(end+1:3) = 0;
T_per_k = T_per.s0 + T_per.s1*lam(1) + T_per.s2*lam(2) + T_per.s3*lam(3)...
        + T_per.s12*lam(1)*lam(2) + T_per.s23*lam(2)*lam(3) + T_per.s13*lam(1)*lam(3)...
        + T_per.s123*lam(1)*lam(2)*lam(3);

% Expand full mode shape to include periodically-constrained DOFs
PHI_plot_full = T_per_k*PHI_plot_full;

% Add mode shape deflections to original coordinate locations
coordinates_phi_full = coordinates;
coordinates_phi_full(:,1) = coordinates_phi_full(:,1) + real(PHI_plot_full(1:2:end));
coordinates_phi_full(:,2) = coordinates_phi_full(:,2) + real(PHI_plot_full(2:2:end));

% plot full model mode
figure(3);clf
[hpatch,hline] = plot_FEM_model(coordinates_phi_full,patchfaces,C,fedges);
axis equal

%% plot BMS model mode shape
% ======================================================================= %

% Expand BMS mode shape to include periodically-constrained DOFs
PHI_plot_BMS = BMS_Plus_Mode_Expansion(PHI_BMS(:,i_phi_plot,k_sel),...
    dof_sets_BMS,kvec,R,T_BMS,K_BMS,M_BMS);

% Align phase of BMS mode shape with full mode shape phase
PHI_plot_BMS = PHI_plot_BMS*(PHI_plot_BMS'*PHI_plot_full);

% normalize BMS mode shape
PHI_plot_BMS = (Lx*mag_PHI)*PHI_plot_BMS/[diag(max(abs(PHI_plot_BMS)))]';

% Add mode shape deflections to original coordinate locations
coordinates_phi_BMS = coordinates;
coordinates_phi_BMS(:,1) = coordinates_phi_BMS(:,1) + real(PHI_plot_BMS(1:2:end));
coordinates_phi_BMS(:,2) = coordinates_phi_BMS(:,2) + real(PHI_plot_BMS(2:2:end));

% plot BMS model mode over top of full model mode (just edge lines)
[hpatchBMS,hlineBMS] = plot_FEM_model(coordinates_phi_BMS,patchfaces,C,fedges);
axis equal
delete(hpatchBMS)
set(hlineBMS,'linestyle','--','color','green')

% BMS mode error
e_phi = 1-abs(PHI_plot_full'*PHI_plot_BMS)/(sqrt(PHI_plot_full'*PHI_plot_full)*sqrt(PHI_plot_BMS'*PHI_plot_BMS));

% add title to plot
title(['Mode Comparison: branch ',num2str(i_phi_plot),', k-point ',num2str(k_sel),' of ',num2str(n_kap),', error=',num2str(e_phi)])

% add legend to plot
legend([hline(1),hlineBMS(1)],'full','BMS');

%% Remove subfolders from Matlab path
% ======================================================================= %
rmpath(genpath([functionPathName,'libraries']))
