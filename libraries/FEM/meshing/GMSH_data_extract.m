% The unified code to import GMSH mesh files
% clear all
% clc
%%%%%%%%%%%%%%%%%%%%%%%%%   Control Box   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,elenodes] = ...
    GMSH_data_extract(filename)
% filename = 'LR.msh';
% tolerance = 0.00001;  % Choose a threshold for "close enough to zero"


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Openning the file to import the nodes
text = '$Nodes';
fid = fopen(filename);                       %opent he file.
tline = fgetl(fid);                          %Read line from file, removing newline characters

% first, find nodes
while ischar(tline)                          %Determine whether item is character array.  returns 0 if worng 1 if right.
    bol = strcmp(tline, text);               %Compare strings with case sensitivity.  If the same then 1 if not then 0.
    if (bol)                                %If(0) then wrong.  If(1) then right and execute.
%         disp(tline)
%         trigger = 1;
        
        % first line of "Nodes" contains the number of nodes
        tline = fgetl(fid); %Getting the first values of the matrix.
        n_nodes = str2num(tline);
        
        tline = fgetl(fid); %Get to next line to skip the node numbers.
        %space that is right after the text line
        %tline = fgetl(fid);  Not essential now.  it was for the original
        %file becasue of the presence of a spcae right after the text.
        
        % determine number of dimensions 
        n_dim = length(str2num(tline))-1;
        
        % preallocate coordinate array and node index vector
        X = zeros(n_nodes,n_dim);
        nodes = zeros(n_nodes,1);
        for i = 1:n_nodes
            
            % extract values stored on current line of file
            linevals = str2num(tline);
            
            % first value is the node number
            nodes(i) = linevals(1);
            
            % second through last values are node coordinates
            X(i,:) = linevals(2:end);
            tline = fgetl(fid);
        end        
    end
    tline = fgetl(fid);
end
fclose(fid);
% Nodes=MatfrixN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ELEMENTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
text = '$Elements';
end_text = '$EndElements';

fid = fopen(filename);                       %opent he file.
tline = fgetl(fid);                          %Read line from file, removing newline characters
% cntr = 0; cntr2 = 0; cntr3 = 0; cntr4 = 0;

while ischar(tline)                          %Determine whether item is character array.  returns 0 if worng 1 if right.
    bol = strcmp(tline, text);               %Compare strings with case sensitivity.  If the same then 1 if not then 0.
    if (bol)                                 %If(0) then wrong.  If(1) then right and execute.
%         disp(tline)
%         trigger = 1;
        
        % first line in elements section gives number of elements
        tline = fgetl(fid); %Getting the first values of the matrix.
        n_ele = str2num(tline);
        
        % get next line
        tline = fgetl(fid);
        
        % preallocate element node table
        elenodes = [];
        max_nodes_in_ele = 0;
        elements = zeros(n_ele,1);
        
        for i = 1:n_ele
            linevals = str2num(tline);
            
            node_start = linevals(3)+4;
            nodes_in_ele = length(linevals)+1-node_start;
            
            % expand element node table if a larger element is encountered
            if nodes_in_ele > max_nodes_in_ele
%                 nodes_in_ele
%                 max_nodes_in_ele
                
                elenodes = [elenodes,nan(n_ele,nodes_in_ele-max_nodes_in_ele)];
                max_nodes_in_ele = nodes_in_ele;
            end
            
            % store element nodes
%             linevals(5:end)
%             node_start
%             size(elenodes)
            elenodes(i,1:nodes_in_ele) = linevals(node_start:end);
            
            % add element number to ele
            elements(i) = linevals(1);
            tline = fgetl(fid);
        end
    end
    tline = fgetl(fid);
end

fclose(fid);