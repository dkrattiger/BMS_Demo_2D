function [xloc_new,yloc_new,zloc_new,elenodes,pattern,patchfaces_lite,C_lite,fedges] = ...
     mesh_spherical_resonator(n,n_ele_sc,n_ele_sph,n_ele_layer,...
         n_ele_bc,r_sc,r_bc,r_s,r_l,one_patch_per_face)

if nargin<10
    one_patch_per_face = true;
end
     
% clear;clc;close all
% 
% % element size
% n = 3;
% 
% % mesh size
% n_ele_sc = 2;
% n_ele_sph = 1;
% n_ele_layer = 1;
% n_ele_bc = 1;
% r_sc = 2;
% r_bc = 5;
% r_s = 4;
% r_l = 4.5;

% 
% n_ele_sc = 4;
% n_ele_sph = 2;
% n_ele_layer = 2;
% n_ele_bc = 2;

% n_ele_sph = n_ele_sph;
% n_ele_sc = n_ele_sc;
% n_ele_sc = n_ele_sc;
n_ele_cube = n_ele_sc^3;
n_ele_side = n_ele_sc*(n_ele_sph+n_ele_layer+n_ele_bc)*n_ele_sc;
n_nodes_cube = (n_ele_sc*(n-1)+1)*(n_ele_sc*(n-1)+1)*(n_ele_sc*(n-1)+1);
n_nodes_side = (n_ele_sc*(n-1)+1)*((n_ele_sph+n_ele_layer+n_ele_bc)*(n-1)+1)*(n_ele_sc*(n-1)+1);

% r_sc = 1;
% r_bc = 8;
% r_s = 4;
% r_l = 4.5;

%% cube nodes;
% ======================================================================= %
vec = linspace(-r_sc,r_sc,(n-1)*n_ele_sc+1);
% [xgrid,ygrid,zgrid] = meshgrid(vec,vec,vec);
[xgrid,ygrid,zgrid] = ndgrid(vec,vec,vec);

xloc_cube = xgrid(:);
yloc_cube = ygrid(:);
zloc_cube = zgrid(:);

%% 1/6 sphere nodes
% ======================================================================= %
thetas = linspace(-pi/4,pi/4,(n-1)*n_ele_sc+1);

xgrid = zeros((n-1)*n_ele_sc+1,(n-1)*(n_ele_sph+n_ele_layer+n_ele_bc)+1,(n-1)*n_ele_sc+1);
ygrid = zeros((n-1)*n_ele_sc+1,(n-1)*(n_ele_sph+n_ele_layer+n_ele_bc)+1,(n-1)*n_ele_sc+1);
zgrid = zeros((n-1)*n_ele_sc+1,(n-1)*(n_ele_sph+n_ele_layer+n_ele_bc)+1,(n-1)*n_ele_sc+1);

for i = 1:(n-1)*n_ele_sc+1
    theta = thetas(i);
    phi_lim = atan(cos(theta));
    phis = linspace(-phi_lim,phi_lim,(n-1)*n_ele_sc+1);
    for j = 1:(n-1)*n_ele_sc+1
        
        phi = phis(j);
        
        % point on cube from which to start radial points
        x_sc = r_sc;
        y_sc = r_sc*(-1+2*(i-1)/((n-1)*n_ele_sc));
        z_sc = r_sc*(-1+2*(j-1)/((n-1)*n_ele_sc));     
        
        % point on sphere on which to end radial points
        x_s = r_s*cos(phi)*cos(theta);
        y_s = r_s*cos(phi)*sin(theta);
        z_s = r_s*sin(phi);
        
        % point on layer on which to end radial points
        x_l = r_l*cos(phi)*cos(theta);
        y_l = r_l*cos(phi)*sin(theta);
        z_l = r_l*sin(phi);
        
        % point on layer on which to end radial points
        x_bc = r_bc;
%         y_bc = r_bc*(-1+2*(i-1)/((n-1)*n_ele_sc));
%         z_bc = r_bc*(-1+2*(j-1)/((n-1)*n_ele_sc));  
        y_bc = r_bc*tan(theta);
        z_bc = r_bc*tan(phi)/cos(theta);   
            
        xvec = [linspace(x_sc,x_s,(n-1)*n_ele_sph+1),...
            linspace(x_s,x_l,(n-1)*n_ele_layer+1),...
            linspace(x_l,x_bc,(n-1)*n_ele_bc+1)];
        xvec([(n-1)*n_ele_sph+2,(n-1)*(n_ele_sph+n_ele_layer)+3]) = [];
        
        
        yvec = [linspace(y_sc,y_s,(n-1)*n_ele_sph+1),...
            linspace(y_s,y_l,(n-1)*n_ele_layer+1),...
            linspace(y_l,y_bc,(n-1)*n_ele_bc+1)];
        yvec([(n-1)*n_ele_sph+2,(n-1)*(n_ele_sph+n_ele_layer)+3]) = [];
        
        zvec = [linspace(z_sc,z_s,(n-1)*n_ele_sph+1),...
            linspace(z_s,z_l,(n-1)*n_ele_layer+1),...
            linspace(z_l,z_bc,(n-1)*n_ele_bc+1)];
        zvec([(n-1)*n_ele_sph+2,(n-1)*(n_ele_sph+n_ele_layer)+3]) = [];

        xgrid(i,:,j) = xvec;
        ygrid(i,:,j) = yvec;
        zgrid(i,:,j) = zvec;
    end
end

xloc_sphr_sixth = xgrid(:);
yloc_sphr_sixth = ygrid(:);
zloc_sphr_sixth = zgrid(:);

plot3(xloc_sphr_sixth,yloc_sphr_sixth,zloc_sphr_sixth,'k.');
for i = 1:n_nodes_side
    text(xloc_sphr_sixth(i),yloc_sphr_sixth(i),zloc_sphr_sixth(i),num2str(i))
end

% pause


% combine six sphere segments
xloc = [xloc_cube;   xloc_sphr_sixth;   -yloc_sphr_sixth;  -xloc_sphr_sixth;  yloc_sphr_sixth;  -zloc_sphr_sixth;   zloc_sphr_sixth];
yloc = [yloc_cube;   yloc_sphr_sixth;   xloc_sphr_sixth;   -yloc_sphr_sixth;  -xloc_sphr_sixth;  yloc_sphr_sixth;   yloc_sphr_sixth];
zloc = [zloc_cube;   zloc_sphr_sixth;   zloc_sphr_sixth;   zloc_sphr_sixth;   zloc_sphr_sixth;   xloc_sphr_sixth;   -xloc_sphr_sixth];
          
node_cube = reshape(1:n_nodes_cube,(n_ele_sc*(n-1)+1),(n_ele_sc*(n-1)+1),(n_ele_sc*(n-1)+1));
i_l1 = flipud(squeeze(node_cube(:,1,:))');
i_r1 = flipud(squeeze(node_cube(:,(n_ele_sc*(n-1)+1),:))');
i_f1 = flipud(squeeze(node_cube(1,:,:))');
i_b1 = flipud(squeeze(node_cube((n_ele_sc*(n-1)+1),:,:))');
i_d1 = flipud(squeeze(node_cube(:,:,1)));
i_t1 = flipud(squeeze(node_cube(:,:,(n_ele_sc*(n-1)+1))));

node_rect = reshape(1:n_nodes_side,(n_ele_sc*(n-1)+1),((n_ele_sph+n_ele_layer+n_ele_bc)*(n-1)+1),(n_ele_sc*(n-1)+1));
i_l2 = flipud(squeeze(node_rect(:,1,:))');
i_r2 = flipud(squeeze(node_rect(:,(n_ele_sph+n_ele_layer+n_ele_bc)*(n-1)+1,:))');
i_f2 = flipud(squeeze(node_rect(1,:,:))');
i_b2 = flipud(squeeze(node_rect((n_ele_sc*(n-1)+1),:,:))');
i_d2 = flipud(squeeze(node_rect(:,:,1)));
i_t2 = flipud(squeeze(node_rect(:,:,(n_ele_sc*(n-1)+1))));

% i_remap1 = zeros(n_nodes,1);
i_remap2 = zeros(n_nodes_side,1);
i_remap3 = zeros(n_nodes_side,1);
i_remap4 = zeros(n_nodes_side,1);
i_remap5 = zeros(n_nodes_side,1);
i_remap6 = zeros(n_nodes_side,1);
i_remap7 = zeros(n_nodes_side,1);

% central cube;
i_remap1 = (1:n_nodes_cube)'; 
n_nodes = max(i_remap1);

% positive x section
i_remap2(i_l2) = i_remap1(i_r1);
i_remap2(~i_remap2) = n_nodes+(1:sum(~i_remap2));
n_nodes = max(i_remap2);

% positive y section
i_remap3(i_l2) = i_remap1(fliplr(i_b1));
i_remap3(i_f2) = i_remap2(i_b2);
i_remap3(~i_remap3) = n_nodes+(1:sum(~i_remap3));
n_nodes = max(i_remap3);

% negative x section
i_remap4(i_l2) = i_remap1(fliplr(i_l1));
i_remap4(i_f2) = i_remap3(i_b2);
i_remap4(~i_remap4) = n_nodes+(1:sum(~i_remap4));
n_nodes = max(i_remap4);

% negative y section
i_remap5(i_l2) = i_remap1(i_f1);
i_remap5(i_f2) = i_remap4(i_b2);
i_remap5(i_b2) = i_remap2(i_f2);
i_remap5(~i_remap5) = n_nodes+(1:sum(~i_remap5));
n_nodes = max(i_remap5);

% positive z section
i_remap6(i_l2) = i_remap1(rot90(i_t1,-1));
i_remap6(i_d2) = i_remap2(i_t2);
i_remap6(i_b2) = i_remap3(i_t2);
i_remap6(i_t2) = i_remap4(flipud(i_t2));
i_remap6(i_f2) = i_remap5(flipud(i_t2));
i_remap6(~i_remap6) = n_nodes+(1:sum(~i_remap6));
n_nodes = max(i_remap6);

% positive z section
i_remap7(i_l2) = i_remap1(fliplr(rot90(i_d1,1)));
i_remap7(i_t2) = i_remap2(i_d2);
i_remap7(i_b2) = i_remap3(flipud(i_d2));
i_remap7(i_d2) = i_remap4(flipud(i_d2));
i_remap7(i_f2) = i_remap5((i_d2));
i_remap7(~i_remap7) = n_nodes+(1:sum(~i_remap7));
n_nodes = max(i_remap7);

i_remap = [i_remap1;i_remap2;i_remap3;i_remap4;i_remap5;i_remap6;i_remap7];

xloc_new = nan(n_nodes,1);
yloc_new = nan(n_nodes,1);
zloc_new = nan(n_nodes,1);

xloc_new(i_remap) = xloc;
yloc_new(i_remap) = yloc;
zloc_new(i_remap) = zloc;

%% define element connectivity matrix
% ======================================================================= %
elenodes_cube = zeros(n_ele_cube,n^3);
nodes_cube = reshape(1:n_nodes_cube,n_ele_sc*(n-1)+1,n_ele_sc*(n-1)+1,n_ele_sc*(n-1)+1);
% pattern_cube = zeros(n_ele_cube,1);
pattern_cube = ones(n_ele_cube,1);

% corners = [1,n^2-n+1,n^2,n];
% sides = [(n+1:n:n^2-2*n+1),(n^2-n+2:1:n^2-1),(n^2-n:-n:n+n),(n-1:-1:2)];
% interior = 1:n^2; interior([corners,sides]) = [];
% flat_elenode_order = [corners,sides,interior];
% elenode_order = [];
% for i = 1:n
%     elenode_order = [elenode_order,flat_elenode_order+(n^2)*(i-1)];
% end
% elenode_order = 1:n^3;


for i = 1:n_ele_sc
    for j = 1:n_ele_sc
        for k = 1:n_ele_sc
            x_bl = (i-1)*(n-1)+1;
            y_bl = (j-1)*(n-1)+1;
            z_bl = (k-1)*(n-1)+1;
            
            elenode_cube = zeros(n,n,n);
            for q = 1:n
                elenode_cube(:,:,q) = nodes_cube([x_bl:x_bl+(n-1)],[y_bl:y_bl+(n-1)],z_bl+q-1);
            end
            
            elenodes_cube((i-1)*n_ele_sc*n_ele_sc+(j-1)*n_ele_sc+k,:) = ...
                elenode_cube(:);
        end
    end
end

elenodes_side = zeros(n_ele_side,n^3);
nodes_side = reshape(1:n_nodes_side,n_ele_sc*(n-1)+1,(n_ele_sph+n_ele_layer+n_ele_bc)*(n-1)+1,n_ele_sc*(n-1)+1);
% pattern_side = zeros(n_ele_side,1);
pattern_side = ones(n_ele_side,1);

for i = 1:n_ele_sc
    for j = 1:(n_ele_sph+n_ele_layer+n_ele_bc)
        for k = 1:n_ele_sc
%             x_bl = (j-1)*(n-1)+1; 
%             y_bl = (i-1)*(n-1)+1;           
%             z_bl = (k-1)*(n-1)+1;
            
            x_bl = (i-1)*(n-1)+1;           
            y_bl = (j-1)*(n-1)+1; 
            z_bl = (k-1)*(n-1)+1;
            
            elenode_cube = zeros(n,n,n);
            for q = 1:n
                elenode_cube(:,:,q) = nodes_side([x_bl:x_bl+(n-1)],[y_bl:y_bl+(n-1)],z_bl+q-1);
            end
            
            elenodes_side((i-1)*(n_ele_sph+n_ele_layer+n_ele_bc)*n_ele_sc+(j-1)*n_ele_sc+k,:) = ...
                elenode_cube(:);
            
            if j > n_ele_sph
                
                % pattern materials: 0 = lead, 1 = epoxy, 2 = rubber
                % pattern materials: 1 = lead, 2 = epoxy, 3 = rubber
                if j > n_ele_sph+n_ele_layer
                    % specify outer cube material
%                     pattern_side((i-1)*(n_ele_sph+n_ele_layer+...
%                         n_ele_bc)*n_ele_sc+(j-1)*n_ele_sc+k) = 1;
                    pattern_side((i-1)*(n_ele_sph+n_ele_layer+...
                        n_ele_bc)*n_ele_sc+(j-1)*n_ele_sc+k) = 2;
                else
                    % specify spherical coating material
%                     pattern_side((i-1)*(n_ele_sph+n_ele_layer+...
%                         n_ele_bc)*n_ele_sc+(j-1)*n_ele_sc+k) = 2;
                    pattern_side((i-1)*(n_ele_sph+n_ele_layer+...
                        n_ele_bc)*n_ele_sc+(j-1)*n_ele_sc+k) = 3;
                end
            end
        end
    end
end
n_elements = n_ele_cube+6*n_ele_side;
elenodes = [elenodes_cube;...
            elenodes_side+1*n_nodes_cube+0*n_nodes_side;...
            elenodes_side+1*n_nodes_cube+1*n_nodes_side;...
            elenodes_side+1*n_nodes_cube+2*n_nodes_side;...
            elenodes_side+1*n_nodes_cube+3*n_nodes_side;...
            elenodes_side+1*n_nodes_cube+4*n_nodes_side;...
            elenodes_side+1*n_nodes_cube+5*n_nodes_side];
pattern = [pattern_cube;pattern_side;pattern_side;pattern_side;...
           pattern_side;pattern_side;pattern_side];

elenodes = i_remap(elenodes);

% for i = 1:size(elenodes,1)
%     figure(999);clf;
%     plot3(xloc_new(elenodes(i,:)),yloc_new(elenodes(i,:)),zloc_new(elenodes(i,:)),'k.')
    

%% Create patches for element faces
% ======================================================================= %
% one_patch_per_face = false;
if one_patch_per_face
%     bottom = [1,5+0*(n-2):4+1*(n-2),...
%               2,5+1*(n-2):4+2*(n-2),...
%               3,5+2*(n-2):4+3*(n-2),...
%               4,5+3*(n-2):4+4*(n-2)];
%     top = (n-1)*n^2+bottom;
%     side12 = [1:n^2:(n-1)*n^2+1,...
%              (n-1)*n^2+5:(n-1)*n^2+4+(n-2),...
%              (n-1)*n^2+2:-n^2:2,...
%              4+(n-2):-1:5];
%     side23 = [2:n^2:(n-1)*n^2+2,...
%               (n-1)*n^2+5+1*(n-2):(n-1)*n^2+4+2*(n-2),...
%               (n-1)*n^2+3:-n^2:3,...
%               4+2*(n-2):-1:5+1*(n-2)];
%     side34 = [3:n^2:(n-1)*n^2+3,...
%               (n-1)*n^2+5+2*(n-2):(n-1)*n^2+4+3*(n-2),...
%               (n-1)*n^2+4:-n^2:4,...
%               4+3*(n-2):-1:5+2*(n-2)];
%     side41 = [4:n^2:(n-1)*n^2+4,...
%               (n-1)*n^2+5+3*(n-2):(n-1)*n^2+4+4*(n-2),...
%               (n-1)*n^2+1:-n^2:1,...
%               4+4*(n-2):-1:5+3*(n-2)];
% 
%     % ensure that normals point outwards
%     bottom = fliplr(bottom);
%     side12 = fliplr(side12);
%     side23 = fliplr(side23);
%     side34 = fliplr(side34);
%     side41  = fliplr(side41);
    
    bottom = [1,2:1:(n-1),...
              n:n:(n^2-n),...
              n^2:-1:(n^2-n+2),...
              (n^2-n+1):-n:(n+1)];
    top = (n-1)*n^2+fliplr(bottom);
    side12 = [1:n^2:((n-1)*n^2+1),...
              (n-1)*n^2+(2:1:n),...
              n + (((n-2)*n^2):-n^2:0),...
              (n-1):-1:2];
    side23 = [n + (0:n^2:n^2*(n-1)),...
             (n^2*(n-1)+2*n):n:(n^3),...
             (n^2*(n-1)):-n^2:n^2,...
             (n^2-n):-n:2*n];
    side34 = [n^2:n^2:n^3,...
              (n^3-1):-1:(n^3-n+1),...
              (n^3-n+1-n^2):-n^2:(n^2-n+1),...
              (n^2-n+2):1:(n^2-1)];
    side41 = [(n^2-n+1):n^2:(n^3-n+1),...
              (n^3-2*n+1):-n:(n^2*(n-1)+1),...
              (n^2*(n-2)+1):-n^2:1,...
              (n+1):n:(n^2-2*n+1)];
    % ensure that normals point outwards
    % bottom = fliplr(bottom);
    top = fliplr(top);
    side12 = fliplr(side12);
    side23 = fliplr(side23);
    side34 = fliplr(side34);
    side41  = fliplr(side41);
    
else
    elenode_cube = zeros(n,n,n);
    elenode_cube(elenode_order) = 1:n^3;

    for i = 1:n-1
        for j = 1:n-1
            bottom((i-1)*(n-1)+j,:) = [elenode_cube(i,j,1),elenode_cube(i,j+1,1),...
                      elenode_cube(i+1,j+1,1),elenode_cube(i+1,j,1)];
            top((i-1)*(n-1)+j,:) = [elenode_cube(i,j,n),elenode_cube(i,j+1,n),...
                      elenode_cube(i+1,j+1,n),elenode_cube(i+1,j,n)];
            side12((i-1)*(n-1)+j,:) = [elenode_cube(1,i,j),elenode_cube(1,i,j+1),...
                      elenode_cube(1,i+1,j+1),elenode_cube(1,i+1,j)];
            side23((i-1)*(n-1)+j,:) = [elenode_cube(i,n,j),elenode_cube(i,n,j+1),...
                      elenode_cube(i+1,n,j+1),elenode_cube(i+1,n,j)];
            side34((i-1)*(n-1)+j,:) = [elenode_cube(n,i,j),elenode_cube(n,i,j+1),...
                      elenode_cube(n,i+1,j+1),elenode_cube(n,i+1,j)];
            side41((i-1)*(n-1)+j,:) = [elenode_cube(i,1,j),elenode_cube(i,1,j+1),...
                      elenode_cube(i+1,1,j+1),elenode_cube(i+1,1,j)];
        end
    end
    
    % ensure that normals point outwards
    bottom = fliplr(bottom);
    side12 = fliplr(side12);
    side23 = fliplr(side23);
    
end
      
six_patches = [bottom;top;side12;side23;side34;side41];
[n_patch_per_ele,n_patch_nodes] = size(six_patches);
patchfaces = zeros(n_elements*n_patch_per_ele,n_patch_nodes);
for i = 1:n_elements
    for j = 1:size(six_patches,1);
        index = (i-1)*size(six_patches,1) + j;
        patchfaces(index,:) = elenodes(i,six_patches(j,:));
     end
end

% coloring vector
pattern6vec = repmat(pattern,[1,n_patch_per_ele]).';
pattern6vec = pattern6vec(:);
C = ([pattern6vec,pattern6vec,pattern6vec])/3;

%% remove faces to create section view
% ======================================================================= %
i_node_keep0 = true(size(xloc_new));
i_node_keep1 = yloc_new>-1e-6;
i_node_keep2 = yloc_new>-1e-6;
% i_node_keep1 = (xloc_new>-1e-6 | yloc_new>-1e-6);
% i_node_keep2 = (xloc_new>-1e-6 | yloc_new>-1e-6);

i_ele_keep0 = all(i_node_keep0(elenodes),2) & pattern(:) == 1;
i_ele_keep1 = all(i_node_keep1(elenodes),2) & pattern(:) == 2;
i_ele_keep2 = all(i_node_keep2(elenodes),2) & pattern(:) == 3;

i_ele_keep = i_ele_keep0 | i_ele_keep1 | i_ele_keep2;
i_patch_keep = repmat(i_ele_keep,1,n_patch_per_ele)';i_patch_keep = i_patch_keep(:);


patchfaces_keep = patchfaces(i_patch_keep,:);
C_keep = C(i_patch_keep,:);
pattern6vec_keep = pattern6vec(i_patch_keep);

%% remove doubled up patches (so that only exterior surface remains
% ======================================================================= %

patchfaces_sorted = sort(patchfaces_keep,2);
i_rem_patch = false(6*n_elements,1);

[~,~,i_recon] = unique(patchfaces_sorted,'rows','stable');
for i = 1:length(i_recon)
    if sum(i_recon == i_recon(i)) > 1
        i_rem_patch(i) = true;
    end
end

% remove all but surface patches
patchfaces_lite = patchfaces_keep;
patchfaces_lite(i_rem_patch,:) = [];
n_patches = size(patchfaces_lite,1);


% % determine which nodes are not in any patches
% i_patch_nodes_lite = unique(patchfaces_lite(:));
% i_node_map = zeros(n_nodes,1);
% i_node_map(i_patch_nodes_lite) = [1:length(i_patch_nodes_lite)];

% renumber patche nodes so they correspond to nodelist once un-involved 
% nodes have been removed
% patchfaces_lite = i_node_map(patchfaces_lite);
coordinates = [xloc_new,yloc_new,zloc_new];
% xyzlocs_lite = xyzlocs(i_patch_nodes_lite,:);
% xyzlocs_lite = xyzlocs;

% Patch Coloring Array
C_lite = C_keep;
C_lite(i_rem_patch,:) = [];

pattern6vec_lite = pattern6vec_keep;
pattern6vec_lite(i_rem_patch) = [];

remove_back_patch = false;
if remove_back_patch
    
    % camera view vector     
    az = -42.5*pi/180;
    el = 22*pi/180;
    camvec = -[cos(el)*sin(az);-cos(el)*cos(az);sin(el)];
    
    face_normals = zeros(n_patches,3);
    for j = 1:n_patches
        
        i = 2;
        vec_dot = -1;
        % if vectors are 180 degrees apart move to next set of nodes in element
        while vec_dot<-0.99

            vec1 = coordinates(patchfaces_lite(j,i-1),:)-coordinates(patchfaces_lite(j,i),:);
            vec2 = coordinates(patchfaces_lite(j,i+1),:)-coordinates(patchfaces_lite(j,i),:);
            vec1 = vec1/norm(vec1);
            vec2 = vec2/norm(vec2);
            vec_dot = dot(vec1,vec2);
            i = i+1;
        end
        
        % surface normal to the current face (normalized to unit length)
        face_normals(j,:) = cross(vec2,vec1);
        face_normals(j,:) = face_normals(j,:)/norm(face_normals(j,:));
    end
     
    % keep patches that are showing some part of their face to the camera
    patch_keep = face_normals*camvec<0.05;
    patch_keep = face_normals(:,1)<-1e-3 | face_normals(:,2)<-1e-3 | face_normals(:,3)>1e-3;
    patchfaces_lite = patchfaces_lite(patch_keep,:);
    C_lite = C_lite(patch_keep,:);
    pattern6vec_lite = pattern6vec_lite(patch_keep,:);
end


%% Compute Edge lines
% ======================================================================= %
fedges = [];

% loop through different materials
for i = 1:max(pattern6vec_lite)
    
    % convert to triangular mesh
    i_mat = find(pattern6vec_lite == (i));
    n_pts = size(patchfaces_lite,2);
    
    tri_faces = zeros(length(i_mat)*(n_pts-2),3);
    
    for j = 1:length(i_mat)
        node_ind = 1:n_pts;
        
        for k = 1:n_pts-2
            
            n_pts_left = n_pts-k+1;
            
            dxdydz =  coordinates(patchfaces_lite(i_mat(j),node_ind([2:n_pts_left,1])),:)-...
                      coordinates(patchfaces_lite(i_mat(j),node_ind),:);            
            dnorm = sqrt(sum(dxdydz.*dxdydz,2));

            sharp = sum(dxdydz.*dxdydz([n_pts_left,1:(n_pts_left-1)],:),2)./(dnorm.*dnorm([n_pts_left,1:(n_pts_left-1)]));
            
            
            [~,i_max] = min(abs(sharp));
            tri_faces((j-1)*(n_pts-2)+k,:) = patchfaces_lite(i_mat(j),node_ind(mod(i_max+[-2:0],n_pts_left)+1));
            
            node_ind(i_max) = [];
            
        end
    end
        % compute and plot model edges
        tria = triangulation(tri_faces,coordinates);
        fedges_i = featureEdges(tria,40*pi/180);
        fedges = [fedges;fedges_i];
    
end



% find unique segments
[~,i_unique,~] = unique(sort(fedges,2),'rows','stable');
fedges = fedges(i_unique,:);

% % connect segments
% connect_segs = true;
% if connect_segs
%     count = 1;
%     fedgecell{count,1} = fedges(1,:)';
%     fedges(1,:) = nan;
%     while any(~isnan(fedges(:)))
% 
%         % see if any segments connect from beginning
%         [i_r1,i_c1] =  find(fedges == fedgecell{count,1}(1));
% 
%         if ~isempty(i_r1);
%             fedgecell{count,1} = [fedges(i_r1(1),(rem(i_c1(1),2))+1); fedgecell{count,1}];
%             fedges(i_r1(1),:) = nan;
%         end
% 
%         % see if any segments connect from end
%         [i_r2,i_c2] =  find(fedges == fedgecell{count,1}(end));
% 
%         if ~isempty(i_r2)        
%             fedgecell{count,1} = [fedgecell{count,1};fedges(i_r2(1),(rem(i_c2(1),2))+1)];
%             fedges(i_r2(1),:) = nan;        
%         end
% 
% 
%         if isempty(i_r1) && isempty(i_r2)
%             count = count + 1;
%             [i_r,~] = find(~isnan(fedges));
%             if ~isempty(i_r)
%                 fedgecell{count,1} = fedges(i_r(1),:)';
%                 fedges(i_r(1),:) = nan;
%             end
%         end
%     end
% 
%     fedges = fedgecell;
% end

% %% plot patches
% % ======================================================================= %
% 
% figure(2);clf
% view(3)
% h1 = patch('faces',patchfaces_lite,...
%     'vertices',xyzlocs,'facecolor','w');
% set(h1,'edgecolor',[1,1,1]*0.85)
% set(h1,'FaceColor','flat')
% set(h1,'FaceVertexCData',C_lite)
% % set(h1,'facealpha',0.5)
% 
% % set(h1,'FaceColor','interp')
% % set(h2,'FaceVertexCData',color_vals)
% set(h1,'CDataMapping','direct')
% 
% % pause
% 
% % plot model edges
% hold on
% hline1 = plot3(xloc_new(fedges)',yloc_new(fedges)',zloc_new(fedges)');
% set(hline1,'color','k')
% set(hline1,'linestyle','-')
% set(hline1,'linewidth',2)
% 
% xlabel('x (m)');ylabel('y (m)');title('Unit Cell');
% axis equal
% 
% set(gcf,'renderer','opengl')