function P = TriangleMeshT3mapT6(coor,conn)
% P = TriangleMeshT3mapT6(coor,conn), creates a projection matrix that
% allows data from a T3 mesh to be projected to a T6 mesh. coor and conn
% here represent the T6 mesh, where the T3 mesh is the same but only using
% conn(:,1:3).

Ne = size(conn,1);
Nn = size(coor,1);

% build the edge-element matrix
% ===========================


% matrix of size [Nn, Nn]

I = transpose(1:Nn);
J1 = zeros(Nn,1);
J2 = zeros(Nn,1);
S = 0.5 * ones(Nn,1);


% for each element
for ke = 1:Ne

    % elemental connectivity
    con = conn(ke,:);

    % map from vertex node to vertex node
    J1(con(1:3),1) = con(1:3);
    J2(con(1:3),1) = con(1:3);

    % map from edge node to vertex node
    J1(con(4:6),1) = con(1:3);
    J2(con(4:6),1) = con([2, 3, 1]);
end

E1 = sparse(I,J1,S,Nn,Nn);
E2 = sparse(I,J2,S,Nn,Nn);
P = E1 + E2;

% get the nodes that are only vertices
I = unique(conn(:,1:3));

% only keep that part of P
P = P(:,I);

