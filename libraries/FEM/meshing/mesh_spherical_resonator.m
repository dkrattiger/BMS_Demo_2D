function [xloc_new,yloc_new,zloc_new,elenodes,pattern,patchfaces_lite,C_lite,fedges] = ...
     mesh_spherical_resonator(n,n_ele_sc,n_ele_sph,n_ele_layer,...
         n_ele_bc,r_sc,r_bc,r_s,r_l,one_patch_per_face)

% The mesh is defined with an internal cube surrounded on all sides by six
% pyramid-esque cubes. The overall mesh then constitutes a larger cube.
     
if nargin<10
    one_patch_per_face = true;
end 

%% Example values to un-comment when debugging
% ======================================================================= %
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
% one_patch_per_face = true;
% 
% r_sc = 2;
% r_bc = 5;
% r_s = 4;
% r_l = 4.5;

%% number of elements along various mesh dimensions
% ======================================================================= %
n_ele_cube = n_ele_sc^3;
n_ele_side = n_ele_sc*(n_ele_sph+n_ele_layer+n_ele_bc)*n_ele_sc;
n_nodes_cube = (n_ele_sc*(n-1)+1)*(n_ele_sc*(n-1)+1)*(n_ele_sc*(n-1)+1);
n_nodes_side = (n_ele_sc*(n-1)+1)*((n_ele_sph+n_ele_layer+n_ele_bc)*(n-1)+1)*(n_ele_sc*(n-1)+1);

%% cube nodes;
% ======================================================================= %
vec = linspace(-r_sc,r_sc,(n-1)*n_ele_sc+1);
[xgrid,ygrid,zgrid] = ndgrid(vec,vec,vec);

xloc_cube = xgrid(:);
yloc_cube = ygrid(:);
zloc_cube = zgrid(:);

%% 1/6 sphere nodes
% ======================================================================= %

xgrid = zeros((n-1)*(n_ele_sph+n_ele_layer+n_ele_bc)+1,(n-1)*n_ele_sc+1,(n-1)*n_ele_sc+1);
ygrid = zeros((n-1)*(n_ele_sph+n_ele_layer+n_ele_bc)+1,(n-1)*n_ele_sc+1,(n-1)*n_ele_sc+1);
zgrid = zeros((n-1)*(n_ele_sph+n_ele_layer+n_ele_bc)+1,(n-1)*n_ele_sc+1,(n-1)*n_ele_sc+1);

% loop through a 2d grid of points and define radial mesh lines which will
% be discretized into a set of nodes
for i = 1:(n-1)*n_ele_sc+1
    for j = 1:(n-1)*n_ele_sc+1
        
        % point on small cube (internal to sphere) from which to start radial points
        x_sc = r_sc;
        y_sc = r_sc*(-1+2*(i-1)/((n-1)*n_ele_sc));
        z_sc = r_sc*(-1+2*(j-1)/((n-1)*n_ele_sc));            
        X_sc = [x_sc;y_sc;z_sc];
        
        % point on external unit cell cube from which to start radial points
        x_bc = r_bc;
        y_bc = r_bc*(-1+2*(i-1)/((n-1)*n_ele_sc));
        z_bc = r_bc*(-1+2*(j-1)/((n-1)*n_ele_sc));            
        X_bc = [x_bc;y_bc;z_bc];
        
        % point on sphere on which to end radial points
        a_s = (r_s-norm(X_bc))/(norm(X_sc)-norm(X_bc));
        X_s = a_s*X_sc + (1-a_s)*X_bc;
        x_s = X_s(1);y_s = X_s(2);z_s = X_s(3);
        
        % point on spherical layer on which to end radial points
        a_l = (r_l-norm(X_bc))/(norm(X_sc)-norm(X_bc));
        X_l = a_l*X_sc + (1-a_l)*X_bc;
        x_l = X_l(1);y_l = X_l(2);z_l = X_l(3);

        % define list of nodal coordinates
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

        xgrid(:,i,j) = xvec;
        ygrid(:,i,j) = yvec;
        zgrid(:,i,j) = zvec;
    end
end

% vectorize the 
xloc_sphr_sixth = xgrid(:);
yloc_sphr_sixth = ygrid(:);
zloc_sphr_sixth = zgrid(:);

% combine six sphere segments
xloc = [xloc_cube;   xloc_sphr_sixth;   -yloc_sphr_sixth;  -xloc_sphr_sixth;  yloc_sphr_sixth;  -zloc_sphr_sixth;   zloc_sphr_sixth];
yloc = [yloc_cube;   yloc_sphr_sixth;   xloc_sphr_sixth;   -yloc_sphr_sixth;  -xloc_sphr_sixth;  yloc_sphr_sixth;   yloc_sphr_sixth];
zloc = [zloc_cube;   zloc_sphr_sixth;   zloc_sphr_sixth;   zloc_sphr_sixth;   zloc_sphr_sixth;   xloc_sphr_sixth;   -xloc_sphr_sixth];

%% identify nodes that overlap using "unique" function (on rounded coordinates)
% ======================================================================= %
coords = [xloc,yloc,zloc];


% round coordinates so that overlapping points are numerically equivalent
% (should probably fix this at some point to account for the scale of the
% coordinate system.)
coordsr = round(10^8*coords)/10^8;

% find unique coordinates
[~,i_reduce,i_expand] = unique(coordsr,'rows');


% stitch together separate meshes to give a unified mesh with no degenerate
% nodes
coords_unique = coords(i_reduce,:);
xloc_new = coords_unique(:,1);
yloc_new = coords_unique(:,2);
zloc_new = coords_unique(:,3);

%% define element connectivity matrix
% ======================================================================= %
elenodes_cube = zeros(n_ele_cube,n^3);
nodes_cube = reshape(1:n_nodes_cube,n_ele_sc*(n-1)+1,n_ele_sc*(n-1)+1,n_ele_sc*(n-1)+1);
% pattern_cube = zeros(n_ele_cube,1);
pattern_cube = ones(n_ele_cube,1);

for i = 1:n_ele_sc
    for j = 1:n_ele_sc
        for k = 1:n_ele_sc
            x_bl = (i-1)*(n-1)+1;
            y_bl = (j-1)*(n-1)+1;
            z_bl = (k-1)*(n-1)+1;
            
            elenode_cube = zeros(n,n,n);
            for q = 1:n                
                elenode_cube(:,:,q) = nodes_cube([x_bl:(x_bl+(n-1))],[y_bl:(y_bl+(n-1))],z_bl+q-1);
            end
            
            elenodes_cube((i-1)*n_ele_sc*n_ele_sc+(j-1)*n_ele_sc+k,:) = ...
                elenode_cube(:);
        end
    end
end

elenodes_side = zeros(n_ele_side,n^3);
nodes_side = reshape(1:n_nodes_side,(n_ele_sph+n_ele_layer+n_ele_bc)*(n-1)+1,n_ele_sc*(n-1)+1,n_ele_sc*(n-1)+1);
pattern_side = ones(n_ele_side,1);

for i = 1:(n_ele_sph+n_ele_layer+n_ele_bc)
    for j = 1:n_ele_sc
        for k = 1:n_ele_sc
            
            x_bl = (i-1)*(n-1)+1;           
            y_bl = (j-1)*(n-1)+1; 
            z_bl = (k-1)*(n-1)+1;
            
            elenode_cube = zeros(n,n,n);
            for q = 1:n
                elenode_cube(:,:,q) = nodes_side([x_bl:x_bl+(n-1)],[y_bl:y_bl+(n-1)],z_bl+q-1);
            end
            
            elenodes_side((i-1)*n_ele_sc*n_ele_sc+(j-1)*n_ele_sc+k,:) = ...
                elenode_cube(:);
            
            if i > n_ele_sph
                
                % pattern materials: 1 = lead, 2 = epoxy, 3 = rubber
                if i > n_ele_sph+n_ele_layer
                    % specify outer cube material
                    pattern_side((i-1)*n_ele_sc*n_ele_sc+(j-1)*n_ele_sc+k) = 2;
                else
                    % specify spherical coating material
                    pattern_side((i-1)*n_ele_sc*n_ele_sc+(j-1)*n_ele_sc+k) = 3;
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

% map the uncoupled sections to the unique set of nodes found earlier
elenodes = i_expand(elenodes);

%% Create patches for element faces
% ======================================================================= %
% one_patch_per_face = false;

% this section keeps just the perimiter nodes for each element. This
% is useful for plotting element edges and simple mesh diagrams, but
% discards some information when plotting interpolated values such as
% stress over the element faces.
%
% For stress, it makes more sense to  break every face up into a set of 
% smaller faces (
if one_patch_per_face
    
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
    elenode_cube(:) = 1:n^3;

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
section_view = true;
if section_view
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
else
    i_patch_keep = true(size(patchfaces,1),1);

end

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
%     size(tri_faces)
%     size(coordinates)
    
    if ~isempty(tri_faces)
        % compute and plot model edges
        tria = triangulation(tri_faces,coordinates);
        fedges_i = featureEdges(tria,40*pi/180);
        fedges = [fedges;fedges_i];
    end
end



% % find unique segments
% [~,i_unique,~] = unique(sort(fedges,2),'rows','stable');
% fedges = fedges(i_unique,:);


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