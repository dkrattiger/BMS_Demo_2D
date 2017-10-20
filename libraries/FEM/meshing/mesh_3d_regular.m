function [xloc,yloc,zloc,elenodes,patchfaces_lite,C_lite,fedges] = ...
     mesh_3d_regular(n,n_ele_x,n_ele_y,n_ele_z,Lx,Ly,Lz,pattern)

% clear;clc;close all
% 
% % element size
% n = 2;
% 
% n_ele_x = 2;
% n_ele_y = 1;
% n_ele_z = 1;
% 
% Lx = 1;
% Ly = 1;
% Lz = 1;


n_ele = n_ele_x*n_ele_y*n_ele_z;
n_nodes = (n_ele_x*(n-1)+1)*(n_ele_y*(n-1)+1)*(n_ele_z*(n-1)+1);

%% cube nodes;
% ======================================================================= %
vecx = linspace(0,Lx,(n-1)*n_ele_x+1);
vecy = linspace(0,Ly,(n-1)*n_ele_y+1);
vecz = linspace(0,Lz,(n-1)*n_ele_z+1);
[xgrid,ygrid,zgrid] = meshgrid(vecx,vecy,vecz);

xloc = xgrid(:);
yloc = ygrid(:);
zloc = zgrid(:);

%% define element connectivity matrix
% ======================================================================= %
elenodes = zeros(n_ele,n^3);
nodes = reshape(1:n_nodes,n_ele_y*(n-1)+1,n_ele_x*(n-1)+1,n_ele_z*(n-1)+1);
pattern = pattern(:);
% pattern = zeros(n_ele,1);

% corners = [1,n^2-n+1,n^2,n];
% sides = [(n+1:n:n^2-2*n+1),(n^2-n+2:1:n^2-1),(n^2-n:-n:n+n),(n-1:-1:2)];
% interior = 1:n^2; interior([corners,sides]) = [];
% flat_elenode_order = [corners,sides,interior];
% elenode_order = [];
% for i = 1:n
%     elenode_order = [elenode_order,flat_elenode_order+(n^2)*(i-1)];
% end

elenode_order = 1:n^3;

for i = 1:n_ele_x
    for j = 1:n_ele_y
        for k = 1:n_ele_z
            
            x_bl = (i-1)*(n-1)+1;
            y_bl = (j-1)*(n-1)+1;
            z_bl = (k-1)*(n-1)+1;
            
            elenode = zeros(n,n,n);
            for q = 1:n
                elenode(:,:,q) = nodes([y_bl:y_bl+(n-1)],[x_bl:x_bl+(n-1)],z_bl+q-1);
            end
            
            elenodes((i-1)*n_ele_y*n_ele_z+(j-1)*n_ele_z+k,:) = ...
                elenode(elenode_order);
                        
            
        end
    end
end

%% Create patches for element faces and plot
% ======================================================================= %
% bottom = [1,5+0*(n-2):4+1*(n-2),...
%           2,5+1*(n-2):4+2*(n-2),...
%           3,5+2*(n-2):4+3*(n-2),...
%           4,5+3*(n-2):4+4*(n-2)];
% top = (n-1)*n^2+bottom;
% side12 = [1:n^2:(n-1)*n^2+1,...
%          (n-1)*n^2+5:(n-1)*n^2+4+(n-2),...
%          (n-1)*n^2+2:-n^2:2,...
%          4+(n-2):-1:5];
% side23 = [2:n^2:(n-1)*n^2+2,...
%           (n-1)*n^2+5+1*(n-2):(n-1)*n^2+4+2*(n-2),...
%           (n-1)*n^2+3:-n^2:3,...
%           4+2*(n-2):-1:5+1*(n-2)];
% side34 = [3:n^2:(n-1)*n^2+3,...
%           (n-1)*n^2+5+2*(n-2):(n-1)*n^2+4+3*(n-2),...
%           (n-1)*n^2+4:-n^2:4,...
%           4+3*(n-2):-1:5+2*(n-2)];
% side41 = [4:n^2:(n-1)*n^2+4,...
%           (n-1)*n^2+5+3*(n-2):(n-1)*n^2+4+4*(n-2),...
%           (n-1)*n^2+1:-n^2:1,...
%           4+4*(n-2):-1:5+3*(n-2)];
      
% figure(121);clf;hold on
% plot3(xloc,yloc,zloc,'.')
% for i = 1:length(xloc)
%     text(xloc(i),yloc(i),zloc(i),num2str(i))
% end
% plot3(xloc(n),yloc(n),zloc(n),'ro')
% view(3)
% pause

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
          
         
    
%          (n-1)*n^2+5:(n-1)*n^2+4+(n-2),...
%          (n-1)*n^2+2:-n^2:2,...
%          4+(n-2):-1:5];
% side23 = [2:n^2:(n-1)*n^2+2,...
%           (n-1)*n^2+5+1*(n-2):(n-1)*n^2+4+2*(n-2),...
%           (n-1)*n^2+3:-n^2:3,...
%           4+2*(n-2):-1:5+1*(n-2)];
% side34 = [3:n^2:(n-1)*n^2+3,...
%           (n-1)*n^2+5+2*(n-2):(n-1)*n^2+4+3*(n-2),...
%           (n-1)*n^2+4:-n^2:4,...
%           4+3*(n-2):-1:5+2*(n-2)];
% side41 = [4:n^2:(n-1)*n^2+4,...
%           (n-1)*n^2+5+3*(n-2):(n-1)*n^2+4+4*(n-2),...
%           (n-1)*n^2+1:-n^2:1,...
%           4+4*(n-2):-1:5+3*(n-2)];
      
% ensure that normals point outwards
% bottom = fliplr(bottom);
top = fliplr(top);
side12 = fliplr(side12);
side23 = fliplr(side23);
side34 = fliplr(side34);
side41  = fliplr(side41);
      
six_patches = [bottom;top;side12;side23;side34;side41];
patchfaces = zeros(n_ele*6,4*(n-1));
for i = 1:n_ele
    for j = 1:size(six_patches,1);
        index = (i-1)*size(six_patches,1) + j;
        patchfaces(index,:) = elenodes(i,six_patches(j,:));
    end
end

% coloring vector
pattern6vec = [pattern,pattern,pattern,pattern,pattern,pattern].';
pattern6vec = pattern6vec(:);
C = ([pattern6vec,pattern6vec,pattern6vec]+1)/3;

%% remove doubled up patches (so that only exterior surface remains
% ======================================================================= %

patchfaces_sorted = sort(patchfaces,2);
i_rem_patch = false(6*n_ele,1);

[~,~,i_recon] = unique(patchfaces_sorted,'rows');
for i = 1:length(i_recon)
    if sum(i_recon == i_recon(i)) > 1
        i_rem_patch(i) = true;
    end
end

% remove all but surface patches
patchfaces_lite = patchfaces;
patchfaces_lite(i_rem_patch,:) = [];

xyzlocs = [xloc,yloc,zloc];

% Patch Coloring Array
C_lite = C;
C_lite(i_rem_patch,:) = [];

pattern6vec_lite = pattern6vec;
pattern6vec_lite(i_rem_patch) = [];

%% Compute Edge lines
% ======================================================================= %
fedges = [];

% loop through different materials
for i = 1:max(pattern6vec_lite)+1
    
    % convert to triangular mesh
    i_mat = find(pattern6vec_lite == (i-1));
    n_pts = size(patchfaces_lite,2);
    
    tri_faces = zeros(length(i_mat)*(n_pts-2),3);
    
    for j = 1:length(i_mat)
        node_ind = 1:n_pts;
        
        for k = 1:n_pts-2
            
            n_pts_left = n_pts-k+1;
            
            dxdydz =  xyzlocs(patchfaces_lite(i_mat(j),node_ind([2:n_pts_left,1])),:)-...
                      xyzlocs(patchfaces_lite(i_mat(j),node_ind),:);            
            dnorm = sqrt(sum(dxdydz.*dxdydz,2));

            sharp = sum(dxdydz.*dxdydz([n_pts_left,1:(n_pts_left-1)],:),2)./(dnorm.*dnorm([n_pts_left,1:(n_pts_left-1)]));
            
            
            [~,i_max] = min(abs(sharp));
            tri_faces((j-1)*(n_pts-2)+k,:) = patchfaces_lite(i_mat(j),node_ind(mod(i_max+[-2:0],n_pts_left)+1));
            
            node_ind(i_max) = [];
            
        end
    end
    
    % compute model edges
    if ~isempty(tri_faces)        
        tria = triangulation(tri_faces,xyzlocs);
        fedges_i = featureEdges(tria,40*pi/180);
        fedges = [fedges;fedges_i];
    end
    
end

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
% hline1 = plot3(xloc(fedges)',yloc(fedges)',zloc(fedges)');
% set(hline1,'color','k')
% set(hline1,'linestyle','-')
% set(hline1,'linewidth',2)
% 
% xlabel('x (m)');ylabel('y (m)');title('Unit Cell');
% axis equal
% 
% set(gcf,'renderer','opengl')
% 
% 
% figure(2323)
% plot3(xloc,yloc,zloc,'k.')
% for i = 1:length(xloc)
%     text(xloc(i),yloc(i),zloc(i),num2str(i))
% end

% break