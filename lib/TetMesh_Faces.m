function [faces, facemap] = TetMesh_Faces(conn)
% faces = TetMesh_Faces(conn),  return the Triangle faces of a tetrahedra
% mesh.
%
% [faces, facemap] = TetMesh_Faces(conn), also return the facemap
% projection matrix


Ne = size(conn,1);
Nf = 4;

facelist = [1, 3, 2; 1 2 4; 2 3 4; 3 1 4];

faces = zeros(Ne * Nf, 3);
facemap = sparse(Ne * Nf, Ne);

% face indices
I = 0 : Nf : (Ne * Nf - 1);
J = 1 : Ne;

for k = 1:Nf
    faces(I + k, :) = conn( :, facelist(k, :) );
    facemap = facemap + sparse( I + k, J, 1, Ne * Nf, Ne);
end