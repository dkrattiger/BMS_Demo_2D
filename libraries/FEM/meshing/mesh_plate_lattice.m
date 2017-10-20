function [coordinates,elenodes,elenode_order,pattern_vec,patchfaces,C,f_edges] = ...
                    mesh_plate_lattice(n,n_ps,Lx,Ly,n_cellx,n_celly,double_up)

%% Specify Rectangular Model configuration
% ======================================================================= %

% "pattern" is a 2D matrix. Every term in "pattern" corresponds to a
% finite element. The value of each term decides which material the element
% will be made of. 0 = fiber (dark material), 1 = matrix (light material)

% model_select = 'big_plus';
model_select = 'square_and_plus_holes';
% model_select = 'square_and_plus_holes_b';
% model_select = 'square_and_plus_holes_with_resonator';
% model_select = 'square_and_plus_holes_with_resonator_b';
switch model_select
    case 'big_plus'
        
        % Dynamics of Lattice Materials Model
        % Square-Plus Hole Mesh
        pattern = nan*ones(n_ps);
        pattern(:,n_ps*(7/16)+1:(9/16)*n_ps) = 1;
        pattern(n_ps*(7/16)+1:(9/16)*n_ps,:) = 1;
    case 'square_and_plus_holes'
        
        % Dynamics of Lattice Materials Model
        % Square-Plus Hole Mesh
        pattern = nan*ones(n_ps);
        pattern(:,n_ps*(7/16)+1:(9/16)*n_ps) = 1;
        pattern(n_ps*(7/16)+1:(9/16)*n_ps,:) = 1;
        pattern(n_ps*(3/16)+1:(13/16)*n_ps,n_ps*(3/16)+1:(13/16)*n_ps) = 1;
        pattern(n_ps*(5/16)+1:(11/16)*n_ps,n_ps*(5/16)+1:(11/16)*n_ps) = nan;
        
    case 'square_and_plus_holes_b'
        
        % Dynamics of Lattice Materials Model
        % Square-Plus Hole Mesh
        pattern = nan*ones(n_ps);
        pattern(:,n_ps*(3/16)+1:(13/16)*n_ps) = 1;
        pattern(n_ps*(3/16)+1:(13/16)*n_ps,:) = 1;
        pattern((1/16)*n_ps+1:(15/16)*n_ps,n_ps*(5/16)+1:(11/16)*n_ps) = nan;
        pattern(n_ps*(5/16)+1:(11/16)*n_ps,(1/16)*n_ps+1:(15/16)*n_ps) = nan;
        
%     case 'square_and_plus_holes_with_resonator'
%         
%         % Dynamics of Lattice Materials Model
%         % Square-Plus Hole with resonator Mesh
%         pattern = nan*ones(n_ps);
%         pattern(:,n_ps*(7/16)+1:(9/16)*n_ps) = 1;
%         pattern(n_ps*(7/16)+1:(9/16)*n_ps,:) = 1;
%         pattern(n_ps*(3/16)+1:(13/16)*n_ps,n_ps*(3/16)+1:(13/16)*n_ps) = 1;
%         pattern(n_ps*(5/16)+1:(11/16)*n_ps,n_ps*(5/16)+1:(11/16)*n_ps) = nan;
%         pattern(n_ps*(6/16)+1:(10/16)*n_ps,n_ps*(6/16)+1:(10/16)*n_ps) = 1;
%         pattern(n_ps*(7/16)+1:(9/16)*n_ps,n_ps*(5/16)+1:n_ps*(6/16)) = 1;
    
    case 'square_and_plus_holes_with_resonator'
        
        % Dynamics of Lattice Materials Model
        % Square-Plus Hole Mesh
        pattern = nan*ones(n_ps);
        pattern(:,n_ps*(3/16)+1:(13/16)*n_ps) = 1;
        pattern(n_ps*(3/16)+1:n_ps*(13/16),:) = 1;
        pattern(n_ps*(1/16)+1:n_ps*(15/16),n_ps*(5/16)+1:(11/16)*n_ps) = nan;
        pattern(n_ps*(5/16)+1:n_ps*(11/16),(1/16)*n_ps+1:(15/16)*n_ps) = nan;
        pattern(n_ps*(2/16)+1:n_ps*(14/16),n_ps*(6/16)+1:(10/16)*n_ps) = 1;
        pattern(n_ps*(6/16)+1:n_ps*(10/16),(2/16)*n_ps+1:(14/16)*n_ps) = 1;
        pattern(n_ps*(6/16)+1:n_ps*(10/16),(1/16)*n_ps+1:(2/16)*n_ps) = 1;
        

    case 'square_and_plus_holes_with_resonator_b'
        
        % Dynamics of Lattice Materials Model
        % Square-Plus Hole with resonator Mesh
        pattern = nan*ones(n_ps);
        pattern(:,n_ps*(7/16)+1:(9/16)*n_ps) = 1;
        pattern(n_ps*(7/16)+1:(9/16)*n_ps,:) = 1;
        pattern(n_ps*(3/16)+1:(13/16)*n_ps,n_ps*(3/16)+1:(13/16)*n_ps) = 1;
        pattern(n_ps*(5/16)+1:(11/16)*n_ps,n_ps*(5/16)+1:(11/16)*n_ps) = nan;
        pattern([n_ps*(0/16)+1:(2/16)*n_ps,n_ps*(14/16)+1:(16/16)*n_ps],[n_ps*(0/16)+1:(6/16)*n_ps,n_ps*(8/16)+1:(16/16)*n_ps]) = 1;
        pattern([n_ps*(0/16)+1:(6/16)*n_ps,n_ps*(10/16)+1:(16/16)*n_ps],[n_ps*(0/16)+1:(2/16)*n_ps,n_ps*(14/16)+1:(16/16)*n_ps]) = 1;
% 
%     case 'square_and_plus_holes_with_resonator_d'
%         
%         % Dynamics of Lattice Materials Model
%         % Square-Plus Hole with resonator Mesh
%         pattern = nan*ones(n_ps);
%         pattern(:,n_ps*(7/16)+1:(9/16)*n_ps) = 1;
%         pattern(n_ps*(7/16)+1:(9/16)*n_ps,:) = 1;
%         pattern(n_ps*(3/16)+1:(13/16)*n_ps,n_ps*(3/16)+1:(13/16)*n_ps) = 1;
%         pattern(n_ps*(5/16)+1:(11/16)*n_ps,n_ps*(5/16)+1:(11/16)*n_ps) = nan;
%         pattern([n_ps*(0/16)+1:(2/16)*n_ps,n_ps*(14/16)+1:(16/16)*n_ps],[n_ps*(3/16)+1:(7/16)*n_ps,n_ps*(8/16)+1:(13/16)*n_ps]) = 1;
%         pattern([n_ps*(3/16)+1:(7/16)*n_ps,n_ps*(8/16)+1:(13/16)*n_ps],[n_ps*(0/16)+1:(2/16)*n_ps,n_ps*(14/16)+1:(16/16)*n_ps]) = 1;

end

%% Create Mesh
% ======================================================================= %

% create a supercell by repeating unit cell
pattern = repmat(pattern,[n_celly,n_cellx]);

% number of elements in each direction (void elements will be removed
% later)
n_ele_x = size(pattern,2);
n_ele_y = size(pattern,1);

% number of nodes (including void elements
n_nodes = (n_ele_y*(n-1)+1)*(n_ele_x*(n-1)+1);  

% number of elements (ignoring air)
n_eles  = sum(sum(~isnan(pattern)));

% 2d grid
xvec = linspace(0,Lx,(n_ele_x*(n-1)+1));
yvec = linspace(0,Ly,(n_ele_y*(n-1)+1));
[xgrid,ygrid] = meshgrid(xvec,yvec);

% vectorize grid
xlocs = reshape(xgrid,n_nodes,1);
ylocs = reshape(ygrid,n_nodes,1);
pattern_vec = reshape(pattern,numel(pattern),1);
pattern_vec(isnan(pattern_vec)) = [];
node_vec = 1:n_nodes;
% nodes = reshape(node_vec,n_ele_y*(n-1)+1,n_ele_x*(n-1)+1);
nodes = reshape(node_vec,n_ele_x*(n-1)+1,n_ele_y*(n-1)+1);

% define element connectivity matrix
elenodes = zeros(n_eles,n^2);

% corners = [1,n,n^2,n^2-n+1];
% sides = [(2:n-1),(n+n:n:n^2-n),(n^2-1:-1:n^2-n+2),(n^2-2*n+1:-n:n+1)];
% interior = 1:n^2; interior([corners,sides]) = [];
% 
corners = [1,n^2-n+1,n^2,n];
sides = [(n^2-2*n+1:-n:n+1),(n^2-1:-1:n^2-n+2),(n+n:n:n^2-n),(2:n-1)];
interior = 1:n^2; interior([corners,sides]) = [];

elenode_order = [corners,sides,interior];
% elenode_order = 1:n^2;




% patchnode_order = [corners,sides];
% patchnode_order = [(1:n-1),(n:n:n^2-n),(n^2:-1:n^2-n+2),(n^2-n+1:-n:n+1)];
% [~,patchnode_order] = sort([corners,sides]);
patchnode_order = [1,5:n+2,2,n+3:2*n,3,2*n+1:3*n-2,4,3*n-1:4*n-4];
% pause


count = 0;
for i = 1:n_ele_x
    for j = 1:n_ele_y
        
        if ~isnan(pattern(j,i))
            count = count + 1;
            y_bl = (j-1)*(n-1)+1;
            x_bl = (i-1)*(n-1)+1;
%             elenode_square = nodes([y_bl:y_bl+(n-1)],[x_bl:x_bl+(n-1)])';
            elenode_square = nodes([x_bl:x_bl+(n-1)],[y_bl:y_bl+(n-1)]);       
            elenodes(count,:) = elenode_square(elenode_order);
        end
    end
end

%% remove void nodes
% ======================================================================= %
% remove void nodes from "xlocs" and "ylocs" and remap numbers in 
% "elenodes"
i_keep = unique(elenodes(:));
i_transfer = zeros(size(xlocs));
i_transfer(i_keep) = 1:length(i_keep);
xlocs = xlocs(i_keep);
ylocs = ylocs(i_keep);
coordinates = [xlocs,ylocs];
n_nodes = length(xlocs);

% renumber elenodes to reflect updated node numbering
elenodes = i_transfer(elenodes);

% list of nodes on element edge (tracing counterclockwise around element)
% boundary = [1,5:n+2,2,n+3:2*n,3,2*n+1:3*n-2,4,3*n-1:4*n-4];

patchfaces = elenodes(:,patchnode_order);

% colormap array
C = [pattern_vec,pattern_vec,pattern_vec];

%% Find Model Edges
% ======================================================================= %

coordinates = [xlocs,ylocs];
% compute triangulation and then find feature edges
tria = triangulation([patchfaces(:,[1,2,3]);patchfaces(:,[3,4,1])],coordinates);
f_edges = featureEdges(tria,60*pi/180);

%% Create a second layer of patches and create side walls to give illusion of thickness
% ======================================================================= %

if double_up
    h = Lx/50;
    xlocs2 = [xlocs;xlocs];
    ylocs2 = [ylocs;ylocs];
    zlocs2 = [ones(n_nodes,1)*h/2;-ones(n_nodes,1)*h/2];
    faces{1} = elenodes(:,patchnode_order);
    faces{2} = elenodes(:,fliplr(patchnode_order))+n_nodes;
    faces{3} = [elenodes(:,patchnode_order(n:-1:1)),...
                elenodes(:,patchnode_order(1:n))+n_nodes];
    faces{4} = [elenodes(:,patchnode_order((2*n-1):-1:n))...
                elenodes(:,patchnode_order(n:(2*n-1)))+n_nodes,];
    faces{5} = [elenodes(:,patchnode_order((3*n-2):-1:(2*n-1))),...
                elenodes(:,patchnode_order((2*n-1):(3*n-2)))+n_nodes];
    faces{6} = [elenodes(:,patchnode_order([1,(4*n-4):-1:(3*n-2)])),...
                elenodes(:,patchnode_order([(3*n-2):(4*n-4),1]))+n_nodes];


    coordinates = [xlocs2,ylocs2,zlocs2];
    C = repmat(C,[6,1]);
    patchfaces = zeros(6*n_eles,length(patchnode_order));
    for i = 1:6
        
        % pad faces with -1 values where lengths are too short (change to
        % nan later)
        if size(faces{i},2)<size(faces{1},2)
            faces{i} = [faces{i},-1*ones(n_eles,size(faces{1},2)-size(faces{i},2))];
        end        
        patchfaces(i:6:end,:) = faces{i};
    end
else
%     faces = elenodes(:,boundary);
    coordinates = [xlocs,ylocs];
    patchfaces = elenodes(:,patchnode_order);
end



%% Remove Doubled up Patches (so only surface remains)
% ======================================================================= %
patchfaces_sorted = sort(patchfaces,2);
i_rem_patch = false(6*n_eles,1);

[~,~,i_recon] = unique(patchfaces_sorted,'rows','stable');
for i = 1:length(i_recon)
    if sum(i_recon == i_recon(i)) > 1
        i_rem_patch(i) = true;
    end
end

% remove all but surface patches
patchfaces(i_rem_patch,:) = [];
patchfaces(patchfaces == -1) = nan;
C(i_rem_patch,:) = [];
n_patches = size(patchfaces,1);

%% remove patches whose facenormals are entirely negative
% ======================================================================= %

remove_back_patch = true;
% remove_back_patch = false;

if remove_back_patch && double_up
    
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

            vec1 = coordinates(patchfaces(j,i-1),:)-coordinates(patchfaces(j,i),:);
            vec2 = coordinates(patchfaces(j,i+1),:)-coordinates(patchfaces(j,i),:);
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
    patch_keep = face_normals*camvec<0.0;
    patchfaces = patchfaces(patch_keep,:);
    C = C(patch_keep,:);
end


