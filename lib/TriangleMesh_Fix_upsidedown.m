function [A, conn] = TriangleMesh_Fix_upsidedown(coor, conn)
% [A, conn] = TriangleMesh_Fix_upsidedown(coor, conn), check if elements
% are upsidedown and return a logical vector A wich is 1 for all upsidedown
% elements. Also return a corrected connectivity matrix conn.

% get the nodal positions per element
p1 = coor(conn(:,1),:);
p2 = coor(conn(:,2),:);
p3 = coor(conn(:,3),:);

% compute the surface normal
n = cross(p2-p1,p3-p1);

% find all normals pointing down
A = n(:,3) < 0;

% flip the triangles
con = [1 3 2];
conn(A,:) = conn(A,con);


