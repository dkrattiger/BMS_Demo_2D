function [xlocs,ylocs,elenodes,pattern,patchfaces,C,fedges] = ...
     mesh_square_resonator(n,n_ele_side,n_ele_incl,L)
         

% clear;clc;close all
% 
% % element size
% n = 3;
% 
% mrl = 12;
% n_ele_ss = 2*mrl;
% n_ele_c = 1*mrl;
% n_ele_layer = 0*mrl;
% n_ele_bs = 1*mrl;
% 
% % % mesh size
% % n_ele_ss = 16;
% % n_ele_c = 4;
% % n_ele_layer = 0;
% % n_ele_bs = 8;
% 
% % model dimensions
% r_ss = 0.1;
% r_bs = 0.5;
% r_c = 0.25;
% r_l = 0.25;
% 
% % % distortion parameters cube
% cir_dist = 0.5;
% sqr_dist = 0;%0.1;
% n_petals = 8;
% % sqr_bow = 0.15;
% theta_offset = 0*pi/8;

% % number of elements and nodes in different model sections
% n_ele_cube = n_ele_ss^2;
% n_ele_side = n_ele_ss*(n_ele_c+n_ele_layer+n_ele_bs);
% n_nodes_square = (n_ele_ss*(n-1)+1)*(n_ele_ss*(n-1)+1);
% n_nodes_side = (n_ele_ss*(n-1)+1)*((n_ele_c+n_ele_layer+n_ele_bs)*(n-1)+1);

n_nodes_side = n_ele_side*(n-1)+1;
n_nodes_incl = n_ele_incl*(n-1)+1;
n_ele = n_ele_side^2;
n_nodes = n_nodes_side^2;
%% cube nodes;
% ======================================================================= %
vec = linspace(0,L,n_nodes_side);
[xgrid,ygrid] = meshgrid(vec,vec);

% % distort cube
% theta_vec = atan2(vec,r_ss);
% rad_dist_vec = (1+sqr_dist*cos(n_petals*(theta_vec-theta_offset)));
% rad_dist_mat = (rad_dist_vec'*rad_dist_vec((n-1)*n_ele_ss+1:-1:1))/rad_dist_vec(end);
% 
% 
% thetas = atan2(ygrid,xgrid);
% radius = sqrt(xgrid.^2+ygrid.^2).*rad_dist_mat;
% bowmat = (1+abs(sqr_bow*cos(2*thetas)));
% xgrid = radius.*cos(thetas).*bowmat;
% ygrid = radius.*sin(thetas).*bowmat;


% x_s = r_c*cos(theta)*(1+distortion*cos(n_petals*theta));
% figure(123)
% bar3(reshape(testmat,[(n-1)*n_ele_ss+1,(n-1)*n_ele_ss+1]))
% % surf(reshape(thetas,[(n-1)*n_ele_ss+1,(n-1)*n_ele_ss+1]))
% plot(xgrid,ygrid,'.k')

xlocs = xgrid(:);
ylocs = ygrid(:);


%% define element connectivity matrix
% ======================================================================= %
elenodes = zeros(n_ele,n^2);
nodes_square = reshape(1:n_nodes,n_nodes_side,n_nodes_side);


% sides = [(n+1:n:n^2-2*n+1),(n^2-n+2:1:n^2-1),(n^2-n:-n:n+n),(n-1:-1:2)];
corners = [1,n,n^2,n^2-n+1];
sides = [2:n-1,2*n:n:n^2-n,(n^2-1:-1:n^2-n+2),(n^2-2*n+1:-n:n+1)];
interior = 1:n^2; interior([corners,sides]) = [];
elenode_order = [corners,sides,interior];

for i = 1:n_ele_side
    for j = 1:n_ele_side
        
        y_bl = (i-1)*(n-1)+1;
        x_bl = (j-1)*(n-1)+1;

        elenodes_temp = nodes_square([x_bl:x_bl+(n-1)],[y_bl:y_bl+(n-1)])';
        elenodes((i-1)*n_ele_side+j,:) = ...
            elenodes_temp(elenode_order);
    end
end

%% Define pattern (identifying element materials)
% ======================================================================= %
pattern = ones(n_ele_side,n_ele_side);
incl_start = floor((n_ele_side-n_ele_incl)/2)+1;
pattern(incl_start:(incl_start+n_ele_incl-1),incl_start:(incl_start+n_ele_incl-1)) = 0;
% pattern
% pause
pattern = pattern(:);

%% Create patches for element faces and plot
% ======================================================================= %
ele_patch_edge = [1,5+0*(n-2):4+1*(n-2),...
                  2,5+1*(n-2):4+2*(n-2),...
                  3,5+2*(n-2):4+3*(n-2),...
                  4,5+3*(n-2):4+4*(n-2)];
ele_patch_full = [ele_patch_edge,4*(n-1)+1:n^2];
      
% six_patches = [bottom;top;side12;side23;side34;side41];
patchfaces = zeros(n_ele,4*(n-1));
patchfaces2 = zeros(n_ele,n^2);
[~,elenode_order_reverse] = sort(elenode_order);
for i = 1:n_ele
    patchfaces(i,:) = elenodes(i,ele_patch_edge);
    patchfaces2(i,:) = elenodes(i,elenode_order_reverse);
end

% coloring vector
C = ([pattern,pattern,pattern]+1)/3;

% node coordinate array
xylocs = [xlocs,ylocs];

% 

%% Compute Edge lines
% ======================================================================= %
fedges = [];

% loop through different materials
% for i = 1:size(patchfaces,2)-2
% i_patch2tria = [size(patchfaces,2)*ones(size(patchfaces,2)-2,1),...
%                 [1:(size(patchfaces,2)-2)]',[2:(size(patchfaces,2)-1)]'];
            
n_tria_patch = 2*(n-1)^2;
i_patch2tria = zeros(n_tria_patch,3);
% i_patch2tria(1,:) = [(size(patchfaces,2)),1,2];
index2=0;
for i = 1:(n-1)
    for j = 1:(n-1)
        index = (i-1)*(n)+j;
        index2 = index2+1;
        i_patch2tria(index2*2-1,:) = ([index,index+1,index+n]);
        i_patch2tria(index2*2,:) = ([index+1,index+n+1,index+n]);
    
    end
end
%         i_patch2tria
% ele_patch_full
    
for i = 1:max(pattern)+1
    
    % convert to triangular mesh
    i_mat = find(pattern == (i-1));
    n_patch_mat = size(patchfaces(i_mat),1);
    trifaces = zeros(n_tria_patch*n_patch_mat,3);
    for j = 1:n_patch_mat
        patch_node_index = patchfaces2(i_mat(j),:);
        
        trifaces((j-1)*n_tria_patch+1:j*n_tria_patch,:) = ...
            patch_node_index(i_patch2tria);
    end
%         trifaces([2:2:end],:) = patchfaces(i_mat,[3,4,1]);
    
        
    % compute and plot model edges
    tria = triangulation(trifaces,xylocs);
    fedges_i = featureEdges(tria,10*pi/180);
    fedges = [fedges;fedges_i];
end

% find unique segments
[~,i_unique,~] = unique(sort(fedges,2),'rows','stable');
fedges = fedges(i_unique,:);

% connect segments
count = 1;
fedgecell{count,1} = fedges(1,:)';
fedges(1,:) = nan;
while any(~isnan(fedges(:)))
    
    % see if any segments connect from beginning
    [i_r1,i_c1] =  find(fedges == fedgecell{count,1}(1));
    
    if ~isempty(i_r1);
        fedgecell{count,1} = [fedges(i_r1(1),(rem(i_c1(1),2))+1); fedgecell{count,1}];
        fedges(i_r1,:) = nan;
    end
    
    % see if any segments connect from end
    [i_r2,i_c2] =  find(fedges == fedgecell{count,1}(end));
    
    if ~isempty(i_r2)        
        fedgecell{count,1} = [fedgecell{count,1};fedges(i_r2(1),(rem(i_c2(1),2))+1)];
        fedges(i_r2(1),:) = nan;        
    end
    
    
    if isempty(i_r1) && isempty(i_r2)
        count = count + 1;
        [i_r,~] = find(~isnan(fedges));
        if ~isempty(i_r)
            fedgecell{count,1} = fedges(i_r(1),:)';
            fedges(i_r(1),:) = nan;
        end
    end
end

fedges = fedgecell;

% %% plot patches
% % ======================================================================= %
% 
% figure(2);clf
% plot_unit_cell(xylocs,patchfaces,C,fedges)