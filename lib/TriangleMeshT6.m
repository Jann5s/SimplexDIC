function [newcoor, newconn] = TriangleMeshT6(coor,conn)
% [coor, conn] = TriangleMeshT6(coor,conn), converts a T3 mesh to a T6 mesh
% by creating new nodes at the centers of all element edges.
%
% coor are the nodal coordinates and conn the connectivity, to plot the
% mesh use for example:
% patch('Vertices',coor,'Faces',conn(:,1:3));
% or
% patch('Vertices',coor,'Faces',conn(:,[1,4,2,5,3,6]));


Ne = size(conn,1);
Nn = size(coor,1);

% a list of triangle edges (node pairs)
edges = [1,2;2,3;3,1];

% build the edge-element matrix
% ===========================


% matrix of size [Nn, Nn]

I = zeros(6*Ne,1);
J = zeros(6*Ne,1);

% for each element
for ke = 1:Ne

    k1 = (1:3) + 6*(ke-1);
    k2 = (4:6) + 6*(ke-1);
    
    % elemental connectivity
    con = conn(ke,:);
    
    % forward edges
    I(k1,1) = con(edges(:,1));
    J(k1,1) = con(edges(:,2));

    % and the reverse edges
    I(k2,1) = con(edges(:,2));
    J(k2,1) = con(edges(:,1));

end

E = sparse(I,J,ones(6*Ne,1),Nn,Nn);

% squash this matrix to get the edges
[I, J] = find(tril(E));
edgelist = [I(:), J(:)];

% a map from the unused node to the new node position
emap = [5, 6, 4];

% create the new nodes
% ===========================

Nedge = size(edgelist,1);

newcoor = zeros(Nn+Nedge,2);
newconn = zeros(Ne,6);
newcoor(1:Nn,:) = coor;
newconn(:,1:3) = conn;
for kn = 1:Nedge
    
    edge = edgelist(kn,:);
    
    % get the original elements
    els = find(sum(ismember(conn,edge),2) == 2);
        
    % get the connectivity of one
    con = conn(els(1),:);
    
    % get the coordinates of one of these elements
    pos = coor(con,:);
    
    % find which coordinates are of the edge we have
    pcon = ismember(con,edge);
        
    % create the new coordinate
    newcoor(kn+Nn,:) = mean(pos(pcon(1,:),:),1);
    
    % add the connectivity
    for k = 1:numel(els)
        i = emap(~ismember(conn(els(k),:),edge));
        newconn(els(k),i) = kn+Nn;
    end
end

