function A = TriangleArea(coor,conn)
%  A = A = TriangleArea(coor,conn), compute the signed area of a polygon
%  using the shoelace formula, the area will be positive if the polygon is
%  counter-clockwise. This function assumes the elements to be T3
%  triangles, for T6 elements, only the vertex nodes are used.

X1 = coor(conn(:,1),1);
X2 = coor(conn(:,2),1);
X3 = coor(conn(:,3),1);

Y1 = coor(conn(:,1),2);
Y2 = coor(conn(:,2),2);
Y3 = coor(conn(:,3),2);

A = X1.*Y2 - X2.*Y1 + X2.*Y3 - X3.*Y2 + X3.*Y1 - X1.*Y3;
