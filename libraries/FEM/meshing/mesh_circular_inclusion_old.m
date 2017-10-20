function [xloc_new,yloc_new,elenodes,pattern,patchfaces,C,fedges] = ...
     mesh_circular_inclusion(n,n_ele_ss,n_ele_c,n_ele_layer,...
         n_ele_bs,r_ss,r_bs,r_c,r_l,cir_dist,n_petals,...
         theta_offset)
         
% number of elements and nodes in different model sections
n_ele_square = n_ele_ss^2;
n_ele_side = n_ele_ss*(n_ele_c+n_ele_layer+n_ele_bs);
n_nodes_square = (n_ele_ss*(n-1)+1)*(n_ele_ss*(n-1)+1);
n_nodes_side = (n_ele_ss*(n-1)+1)*((n_ele_c+n_ele_layer+n_ele_bs)*(n-1)+1);


%% cube nodes;
% ======================================================================= %
vec = linspace(-r_ss,r_ss,(n-1)*n_ele_ss+1);
[xgrid,ygrid] = meshgrid(vec,vec);
% [xgrid,ygrid] = ndgrid(vec,vec);

xloc_cube = xgrid(:);
yloc_cube = ygrid(:);

%% 1/4 circle nodes
% ======================================================================= %
thetas = linspace(-pi/4,pi/4,(n-1)*n_ele_ss+1);

xgrid = zeros((n-1)*n_ele_ss+1,(n-1)*(n_ele_c+n_ele_layer+n_ele_bs)+1);
ygrid = zeros((n-1)*n_ele_ss+1,(n-1)*(n_ele_c+n_ele_layer+n_ele_bs)+1);

for i = 1:(n-1)*n_ele_ss+1
    theta = thetas(i);
    
    % point on cube from which to start radial points
    x_sc = r_ss;
    y_sc = r_ss*(-1+2*(i-1)/((n-1)*n_ele_ss)); 
    
    
    % point on sphere on which to end radial points
    x_s = r_c*cos(theta);
    y_s = r_c*sin(theta);

    % point on layer on which to end radial points
    x_l = r_l*cos(theta);
    y_l = r_l*sin(theta);
    

    % point on layer on which to end radial points
    x_bc = r_bs;
    y_bc = r_bs*tan(theta);  

    xvec = [linspace(x_sc,x_s,(n-1)*n_ele_c+1),...
            linspace(x_s,x_l,(n-1)*n_ele_layer+1),...
            linspace(x_l,x_bc,(n-1)*n_ele_bs+1)];
    xvec([(n-1)*n_ele_c+2,(n-1)*(n_ele_c+n_ele_layer)+3]) = [];


    yvec = [linspace(y_sc,y_s,(n-1)*n_ele_c+1),...
            linspace(y_s,y_l,(n-1)*n_ele_layer+1),...
            linspace(y_l,y_bc,(n-1)*n_ele_bs+1)];
    yvec([(n-1)*n_ele_c+2,(n-1)*(n_ele_c+n_ele_layer)+3]) = [];

    xgrid(i,:) = xvec;
    ygrid(i,:) = yvec;
end

xloc_sphr_sixth = xgrid(:);
yloc_sphr_sixth = ygrid(:);

% combine six sphere segments
xloc = [xloc_cube;   xloc_sphr_sixth;   -yloc_sphr_sixth;  -xloc_sphr_sixth;  yloc_sphr_sixth];
yloc = [yloc_cube;   yloc_sphr_sixth;   xloc_sphr_sixth;   -yloc_sphr_sixth;  -xloc_sphr_sixth];

node_cube = reshape(1:n_nodes_square,(n_ele_ss*(n-1)+1),(n_ele_ss*(n-1)+1));
i_l1 = squeeze(node_cube(:,1))';
i_c1 = squeeze(node_cube(:,(n_ele_ss*(n-1)+1)))';
i_d1 = squeeze(node_cube(1,:));
i_t1 = (squeeze(node_cube((n_ele_ss*(n-1)+1),:)));

node_rect = reshape(1:n_nodes_side,(n_ele_ss*(n-1)+1),((n_ele_c+n_ele_layer+n_ele_bs)*(n-1)+1));
i_l2 = flipud(squeeze(node_rect(:,1))');
i_c2 = flipud(squeeze(node_rect(:,(n_ele_c+n_ele_layer+n_ele_bs)*(n-1)+1))');
i_d2 = flipud(squeeze(node_rect(1,:))');
i_t2 = flipud(squeeze(node_rect((n_ele_ss*(n-1)+1),:))');

i_remap2 = zeros(n_nodes_side,1);
i_remap3 = zeros(n_nodes_side,1);
i_remap4 = zeros(n_nodes_side,1);
i_remap5 = zeros(n_nodes_side,1);

% central cube;
i_remap1 = (1:n_nodes_square)';
n_nodes = max(i_remap1);

% positive x section
i_remap2(i_l2) = i_remap1(i_c1);
i_remap2(~i_remap2) = n_nodes+(1:sum(~i_remap2));
n_nodes = max(i_remap2);

% positive y section
i_remap3(i_l2) = i_remap1(i_t1((n_ele_ss*(n-1)+1):-1:1));
i_remap3(i_d2) = i_remap2(i_t2);
i_remap3(~i_remap3) = n_nodes+(1:sum(~i_remap3));
n_nodes = max(i_remap3);

% negative x section
i_remap4(i_l2) = i_remap1(i_l1((n_ele_ss*(n-1)+1):-1:1));
i_remap4(i_d2) = i_remap3(i_t2);
i_remap4(~i_remap4) = n_nodes+(1:sum(~i_remap4));
n_nodes = max(i_remap4);

% negative y section
i_remap5(i_l2) = i_remap1(i_d1);
i_remap5(i_d2) = i_remap4(i_t2);
i_remap5(i_t2) = i_remap2(i_d2);
i_remap5(~i_remap5) = n_nodes+(1:sum(~i_remap5));
n_nodes = max(i_remap5);

i_remap = [i_remap1;i_remap2;i_remap3;i_remap4;i_remap5];

xloc_new = nan(n_nodes,1);
yloc_new = nan(n_nodes,1);

xloc_new(i_remap) = xloc;
yloc_new(i_remap) = yloc;

% i_remap
% 
% figure(888);clf
% plot(xloc_new,yloc_new,'k.');hold on
% for i = 1:length(xloc_new)
%     text(xloc_new(i),yloc_new(i),num2str(i));
% end
% pause

%% Distort Nodes
% ======================================================================= %

thetas = atan2(yloc_new,xloc_new);
radius = sqrt(xloc_new.^2+yloc_new.^2);

% compute radial distortion
r_end = abs(r_bs./cos(thetas));
r_end(abs(xloc_new)<abs(yloc_new)) = abs(r_bs./sin(thetas(abs(xloc_new)<abs(yloc_new))));
rad_dist = radius/r_c;
rad_dist(radius>r_c) = -(radius(radius>r_c)-r_c)./(r_end(radius>r_c)-r_c)+1;

% was using this squared version but for very strong elliptic distortion it
% caused 
% rad_dist = rad_dist.^2;
% rad_dist = rad_dist.^(1/2);

% r_end = abs(r_bs./cos(thetas));
% r_end(abs(xloc_new)<abs(yloc_new)) = abs(r_bs./sin(thetas(abs(xloc_new)<abs(yloc_new))));
% rad_dist = radius/r_c;
% rad_dist = -(radius-r_c)./(r_end-r_c)+1;
% rad_dist = rad_dist.^2;

% choose model type
model_select = 1;
switch model_select
    case 1 % helicoid catenoid ( use cir_dist = 0.500, n_petals=8, theta_offset=0) 
        theta_dist = (1-cir_dist^2)*(1+(cir_dist)*1./(1+(cir_dist)*cos(n_petals*(thetas-theta_offset))))-1;
    case 2 % sinusoidal variation
        theta_dist = (cir_dist*cos(n_petals*(thetas-theta_offset)));
    case 3 % elliptical
        a = 1+cir_dist; 
        b = 1-cir_dist;        
        theta_dist = ((a*b)./sqrt(b^2*cos(floor(n_petals/2)*(thetas-theta_offset)).^2 + a^2*sin(floor(n_petals/2)*(thetas-theta_offset)).^2))-1;
    case 4 % potato ( use cir_dist = 0.3500?, n_petals=8, theta_offset=0) 
        rand_vals3 = randn(n_petals,2);
        rand_vals3 = [ 0.87132       1.2884
                      -2.2002       1.5935
                     -0.26545       2.3494
                     -0.36676     -0.66763
                     -0.66521      0.67268
                       1.2276      0.50229
                      0.77345       1.0332
                      0.17259       1.7022];
        rand_vals = rand_vals3;
        theta_dist = zeros(size(thetas));
        for i = 1:n_petals
            amplitude = rand_vals(i,1)+1i*rand_vals(i,2);
            theta_dist = theta_dist + (cir_dist/n_petals)*real(amplitude*exp(1i*(i)*(thetas-theta_offset)));
        end
end


% for delta = linspace(0,1)
%     radiustemp = radius.*(1+delta*rad_dist.*theta_dist);
% 
%     xloc_newtemp = radiustemp.*cos(thetas);
%     yloc_newtemp = radiustemp.*sin(thetas);
% 
%     figure(989)
%     plot(xloc_newtemp,yloc_newtemp,'k.')
%     axis equal
%     pause(0.1);
% end


% % radius = radius.*(1+rad_dist).*(1+theta_dist);
% radius = radius.*(1+rad_dist.*theta_dist);
% xloc_new = radius.*cos(thetas);
% yloc_new = radius.*sin(thetas);
% 
% figure(898)
% plot(xloc_new,yloc_new,'k.')
% axis equal
% pause

radius = radius.*(1+rad_dist.*theta_dist);

xloc_new = radius.*cos(thetas);
yloc_new = radius.*sin(thetas);

%% define element connectivity matrix
% ======================================================================= %
elenodes_square = zeros(n_ele_square,n^2);
nodes_square = reshape(1:n_nodes_square,n_ele_ss*(n-1)+1,n_ele_ss*(n-1)+1);
pattern_square = zeros(n_ele_square,1);

corners = [1,n,n^2,n^2-n+1];

% sides = [(n+1:n:n^2-2*n+1),(n^2-n+2:1:n^2-1),(n^2-n:-n:n+n),(n-1:-1:2)];
sides = [2:n-1,2*n:n:n^2-n,(n^2-1:-1:n^2-n+2),(n^2-2*n+1:-n:n+1)];
interior = 1:n^2; interior([corners,sides]) = [];
elenode_order = [corners,sides,interior];
elenode_order = 1:n^2;

for i = 1:n_ele_ss
    for j = 1:n_ele_ss
        
        y_bl = (i-1)*(n-1)+1;
        x_bl = (j-1)*(n-1)+1;

        elenode_cube = nodes_square([x_bl:x_bl+(n-1)],[y_bl:y_bl+(n-1)])';
        elenodes_square((i-1)*n_ele_ss+j,:) = ...
            elenode_cube(elenode_order);
    end
end

elenodes_side = zeros(n_ele_side,n^2);
nodes_side = reshape(1:n_nodes_side,n_ele_ss*(n-1)+1,(n_ele_c+n_ele_layer+n_ele_bs)*(n-1)+1);
pattern_side = zeros(n_ele_side,1);

for i = 1:n_ele_ss
    for j = 1:(n_ele_c+n_ele_layer+n_ele_bs)
        x_bl = (i-1)*(n-1)+1; 
        y_bl = (j-1)*(n-1)+1; 

        elenode_cube = nodes_side([x_bl:x_bl+(n-1)],[y_bl:y_bl+(n-1)])';

        elenodes_side((i-1)*(n_ele_c+n_ele_layer+n_ele_bs)+j,:) = ...
            elenode_cube(elenode_order);

        if j > n_ele_c

            % pattern materials: 0 = lead, 1 = epoxy, 2 = rubber
            if j > n_ele_c+n_ele_layer
                % specify outer cube material
                pattern_side((i-1)*(n_ele_c+n_ele_layer+...
                    n_ele_bs)+j) = 1;
            else
                % specify spherical coating material
                pattern_side((i-1)*(n_ele_c+n_ele_layer+...
                    n_ele_bs)+j) = 2;
            end
        end
    end
end
n_elements = n_ele_square+4*n_ele_side;
elenodes = [elenodes_square;...
            elenodes_side+1*n_nodes_square+0*n_nodes_side;...
            elenodes_side+1*n_nodes_square+1*n_nodes_side;...
            elenodes_side+1*n_nodes_square+2*n_nodes_side;...
            elenodes_side+1*n_nodes_square+3*n_nodes_side];
pattern = [pattern_square;pattern_side;pattern_side;pattern_side;pattern_side];

elenodes = i_remap(elenodes);

%% Create patches for element faces and plot
% ======================================================================= %
% ele_patch_edge = [1,5+0*(n-2):4+1*(n-2),...
%                   2,5+1*(n-2):4+2*(n-2),...
%                   3,5+2*(n-2):4+3*(n-2),...
%                   4,5+3*(n-2):4+4*(n-2)];
ele_patch_edge = [1:n-1,...
                  n:n:(n^2-n),...
                  n^2:-1:(n^2-n+2),...
                  (n^2-n+1):-n:(n+1)];
% 1:n-1
% n:n:(n^2-n)
% n^2:-1:(n^2-n+2)
% (n^2-n+1):-n:(n+1)
% pause
ele_patch_full = [ele_patch_edge,4*(n-1)+1:n^2];
      
% six_patches = [bottom;top;side12;side23;side34;side41];
patchfaces = zeros(n_elements,4*(n-1));
patchfaces2 = zeros(n_elements,n^2);
[~,elenode_order_reverse] = sort(elenode_order);
for i = 1:n_elements
    patchfaces(i,:) = elenodes(i,ele_patch_edge);
    patchfaces2(i,:) = elenodes(i,elenode_order_reverse);
end

% coloring vector
C = ([pattern,pattern,pattern]+1)/3;



% node coordinate array
xylocs = [xloc_new,yloc_new];

%% Define Boundary sets
% ======================================================================= %

node_rect = reshape(1:n_nodes_side,(n_ele_ss*(n-1)+1),((n_ele_c+n_ele_layer+n_ele_bs)*(n-1)+1));

i_r3 = squeeze(node_rect(2:n_ele_ss*(n-1),(n_ele_c+n_ele_layer+n_ele_bs)*(n-1)+1));
i_rd3 = squeeze(node_rect(1,(n_ele_c+n_ele_layer+n_ele_bs)*(n-1)+1));
i_rt3 = squeeze(node_rect(n_ele_ss*(n-1)+1,(n_ele_c+n_ele_layer+n_ele_bs)*(n-1)+1));

i_l4 = i_remap4(i_r3((n_ele_ss*(n-1)-1):-1:1));i_l4 = i_l4(:);
i_r4 = i_remap2(i_r3);i_r4 = i_r4(:);
i_d4 = i_remap5(i_r3);i_d4 = i_d4(:);
i_t4 = i_remap3(i_r3((n_ele_ss*(n-1)-1):-1:1));i_t4 = i_t4(:);

i_ld4 = i_remap4(i_rt3);
i_rd4 = i_remap2(i_rd3);
i_rt4 = i_remap2(i_rt3);
i_lt4 = i_remap4(i_rd3);


%% Compute Edge lines
% ======================================================================= %
fedges = [];
n_tria_patch = 2*(n-1)^2;
i_patch2tria = zeros(n_tria_patch,3);
index2=0;
for i = 1:(n-1)
    for j = 1:(n-1)
        index = (i-1)*(n)+j;
        index2 = index2+1;
        i_patch2tria(index2*2-1,:) = ([index,index+1,index+n]);
        i_patch2tria(index2*2,:) = ([index+1,index+n+1,index+n]);
    
    end
end
    
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