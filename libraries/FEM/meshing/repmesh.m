function [coordinates,elenodes,patches,pattern,C,R,fedges,Phi] = ...
    repmesh(coordinates,elenodes,patches,pattern,C,R,fedges,Phi,reps)

% figure(1);clf
% plot(coordinates(:,1),coordinates(:,2),'r.');hold on

% original mesh size
[n_nodes,n_dim] = size(coordinates);

% element properties
[n_elements,nodes_per_ele] = size(elenodes);

% reps can at most have as many dimensions as their are coordinates
reps = reps(1:n_dim);

% repeat pattern
pattern = repmat(pattern,[prod(reps),1]);

% repeat pattern
Phi = repmat(Phi,[prod(reps),1]);


% repeat color matrix
C = repmat(C,[prod(reps),1]);

% repeat node coordinates
n_nodes_temp = n_nodes;

% loop through repeat directions
for i = 1:length(reps)
    
    coordinates_temp = zeros(reps(i)*n_nodes_temp,n_dim);
    % repeat coordinates some number of times
    for j = 1:reps(i)
        for k = 1:n_dim
            coordinates_temp((j-1)*n_nodes_temp+1:j*n_nodes_temp,k) = ...
                coordinates(:,k) + (j-1)*R(k,i);
            
        end
    end
    coordinates = coordinates_temp;
    n_nodes_temp = size(coordinates,1);
end

% Find nodes that aren't unique so they can be connected
% ======================================================================= %
coordinates_r = coordinates;
for i = 1:n_dim
    coordinates_r(:,i) = coordinates_r(:,i)/(max(coordinates_r(:,i))-min(coordinates_r(:,i)));
end
coordinates_r = round(coordinates_r*1e10)/1e10;
[~,i_reduce,i_expand] = unique(coordinates_r,'rows','stable');
i_reduce = i_reduce';
i_expand = i_expand';
coordinates = coordinates(i_reduce,:);

%% repeat elenodes array
% ======================================================================= %
% element properties
[n_elements,nodes_per_ele] = size(elenodes);

elenodes_temp = zeros(prod(reps)*n_elements,nodes_per_ele);
for i = 1:prod(reps)
    elenodes_temp((i-1)*n_elements+1:i*n_elements,:) = elenodes + (i-1)*n_nodes;
end
elenodes = elenodes_temp;    
elenodes = i_expand(elenodes);

%% repeat elenodes array
% ======================================================================= %

if iscell(fedges)
    n_edges = length(fedges);
    fedges_temp = cell(prod(reps)*n_edges,1);
    for i = 1:prod(reps)
        for j = 1:n_edges
            fedges_temp{(i-1)*n_edges+j} = fedges{j} + (i-1)*n_nodes;
        end
    end
    fedges = fedges_temp;
    for i = 1:length(fedges)
        fedges{i} = i_expand(fedges{i});
    end
    
else
    [n_edges,nodes_per_edge] = size(fedges);
    fedges_temp = zeros(prod(reps)*n_edges,nodes_per_edge);
    for i = 1:prod(reps)
        fedges_temp((i-1)*n_edges+1:i*n_edges,:) = fedges + (i-1)*n_nodes;
    end
%     for i = 1:length(fedges)
    fedges = i_expand(fedges);
%     end
end
% elenodes = i_expand(elenodes);

%% repeat patches array
% ======================================================================= %

% element properties
[n_patches,nodes_per_patch] = size(patches);

patches_temp = zeros(prod(reps)*n_patches,nodes_per_patch);
for i = 1:prod(reps)
    patches_temp((i-1)*n_patches+1:i*n_patches,:) = patches + (i-1)*n_nodes;
end
patches = patches_temp;
patches = i_expand(patches);


for i = 1:n_dim
    R(:,i) = R(:,i)*reps(i);
end
    
%% Find nodes that aren't unique so they can be connected
% ======================================================================= %
% coordinates_r = coordinates;
% for i = 1:n_dim
%     coordinates_r(:,i) = coordinates_r(:,i)/(max(coordinates_r(:,i))-min(coordinates_r(:,i)));
% end
% coordinates_r = round(coordinates_r*1e10)/1e10;

% [~,i_reduce,i_expand] = unique(coordinates_r,'rows','stable');
% i_reduce = i_reduce';
% i_expand = i_expand';

% coordinates = coordinates(i_reduce,:);
% elenodes = i_expand(elenodes);
% patches = i_expand(patches);

if ~isempty(Phi)
    i_reduce2 = node2dof(i_reduce,2);
    max(max(i_reduce2))
    size(Phi)
    Phi = Phi(i_reduce2,:,:);
end