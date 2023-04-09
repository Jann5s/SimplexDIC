function [A, conn] = TetMesh_Fix_insideout(coor, conn)
% [A, conn] = TetMesh_Fix_insideout(coor, conn), check if elements
% are insideout and return a logical vector A wich is 1 for all insideout
% elements. Also return a corrected connectivity matrix conn.

% get the nodal positions per element
V = TetVolume(coor, conn);

A = V < 0;

% flip the triangles
con = [1 3 2 4];
conn(A,:) = conn(A,con);


