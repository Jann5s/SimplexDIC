function nodes = TriangleContour(coor,conn)
% nodes = TriangleContour(coor,conn), returns all nodes that are on the
% contour of the mesh

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
[I, ~] = find(tril(E)==1);

% the contour nodes:
nodes = I;
