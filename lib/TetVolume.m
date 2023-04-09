function V = TetVolume(coor,conn)
%  V = TetVolume(coor,conn), compute the signed volume of a tetrahedron
%  using the shoelace formula, the volume will be positive if the nodes of
%  the tet are ordered correctly.

X1 = coor(conn(:,1),:);
X2 = coor(conn(:,2),:);
X3 = coor(conn(:,3),:);
X4 = coor(conn(:,4),:);

% use vertex 4 as the origin
A = X1 - X4;
B = X2 - X4;
C = X3 - X4;

% volume is det(|a,b,c|)/6
V = (1/6) * (...
    A(:,1) .* B(:,2) .* C(:,3) + ...
    B(:,1) .* C(:,2) .* A(:,3) + ...
    C(:,1) .* A(:,2) .* B(:,3) - ...
    A(:,3) .* B(:,2) .* C(:,1) - ...
    B(:,3) .* C(:,2) .* A(:,1) - ...
    C(:,3) .* A(:,2) .* B(:,1) );
