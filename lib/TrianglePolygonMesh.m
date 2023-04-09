function [coor, conn] = TrianglePolygonMesh(P,S)
% [coor, conn] = TrianglePolygonMesh(P,S), creates a T3 mesh inside a
% polygon defined by P with the approximate element size of S.
%
% Optionally, P can be a cell-array and S a corresponding size matrix. In
% that case:
% P{1}  : outer contour
% P{2}  : inner contour (e.g. a hole)
% P{3+} : lines along which to force node placement
%
% The bulk element size will be the last value in S, so S can have one
% element more than P to indicate the bulk element size seperately.
%
% This code starts with a simple gridded distribution of nodes that are
% triangulated using delaunay. Then, a simple optimization routine is run
% to more optimally distribute the nodes by minimizing the "electro static
% tension" between them.


if ~iscell(P)
    P = {P};
end

if numel(S) == 1
    S = repmat(S,numel(P)+1,1);
end

% linear interpolation function
lerp = @(x,x1,x2) bsxfun(@times,(1-x),x1) + bsxfun(@times,x,x2);

% Generating a grid of nodes for the bulk
% ======================

% get a bounding box around the polygon
lim = [min(P{1},[],1); max(P{1},[],1)];

% create a grid of nodes inside the domain
Sbulk = S(end);
Nx = ceil((lim(2,1)-lim(1,1))./Sbulk);
Ny = ceil((lim(2,2)-lim(1,2))./Sbulk);
x = linspace(lim(1,1),lim(2,1),Nx+1);
y = linspace(lim(1,2),lim(2,2),Ny+1);
[X, Y] = meshgrid(x(2:Nx),y(2:Ny));

% find points inside
In = polygon_inside(X,Y,P{1}) ~= 0;
% is there a hole
if numel(P) >= 2    
    In = In & polygon_inside(X,Y,P{2}) == 0;
end

% bulk points
bulk = [X(In), Y(In)];

% Generating nodes on the polygon
% ======================

% Create nodes on the polygon edges
edge = [];
w = [];
for l = 1:numel(P)
    % P = polygon_clockwise(P,'fix');
    Np = size(P{l},1);
    
    % for extra edges, don't loop around
    if l > 2
        Np = Np - 1;
    end
    
    % all the edge points for this polygon
    edgep = [];
        
    % loop around vector
    I = [1:Np, 1];
    for k = 1:Np
        p0 = P{l}(I(k),:);
        p1 = P{l}(I(k+1),:);
        dp = p1 - p0;
        L = hypot(dp(1),dp(2));
        N = ceil(L/S(l));
        N = max(N,2);
        t = linspace(0,1,N).';
        t = t(1:N-1);
        edgep = cat(1,edgep,lerp(t,p0,p1));
    end
    
    % for the extra edges, clean up nodes outside the domain
    if l > 2
        % find points inside
        In = polygon_inside(edgep(:,1),edgep(:,2),P{1}) ~= 0;
        % is there a hole
        if numel(P) >= 2
            In = In & polygon_inside(edgep(:,1),edgep(:,2),P{2}) == 0;
        end
        edgep = edgep(In,:);
    end
    
    % add them all the the list of edge nodes
    edge = cat(1,edge,edgep);
    w = cat(1,w,S(l)*ones(size(edgep,1),1));
end


% number of nodes per group
Nedge = size(edge,1);
Nbulk = size(bulk,1);

% indices for each one
Iedge = 1:Nedge;
Ibulk = Nedge + (1:Nbulk);

% combine them together
coor = cat(1,edge,bulk);
w = cat(1,w,Sbulk*ones(Nbulk,1));
w = w ./ max(w);
Nn = Nedge + Nbulk;

% improve the mesh by relaxing it
maxit = 100;
convcrit = 1e-5;
step = 0.1;
conn = delaunay(coor);
conn = clean_zero_area(coor,conn);
conn = convex_elements(P,coor,conn);
Ne = size(conn,1);

ind1 = [2,3,1];
ind2 = [3,1,2];
for it = 1:maxit
    F = zeros(Nn,2);
    for ke = 1:Ne
        con = conn(ke,:);
        con1 = con(ind1);
        con2 = con(ind2);
        F(con,:) = F(con,:) + coor(con1,:)-coor(con,:) + coor(con2,:)-coor(con,:);
    end
    coor(Ibulk,:) = coor(Ibulk,:) + step * [w(Ibulk),w(Ibulk)] .* F(Ibulk,:);
    if mod(it,10) == 0
        conn = delaunay(coor);
        conn = clean_zero_area(coor,conn);
        conn = convex_elements(P,coor,conn);
        Ne = size(conn,1);
    end
    if rms([F(Ibulk,1);F(Ibulk,2)]) < convcrit
        break
    end    
end

% fix any elements that are upside down
conn = fix_upsidedown(coor,conn);
end

function conn = convex_elements(P,coor,conn)
% clean elements that have their centers outside of the polygon
X = (1/3)*(coor(conn(:,1),1) + coor(conn(:,2),1) + coor(conn(:,3),1));
Y = (1/3)*(coor(conn(:,1),2) + coor(conn(:,2),2) + coor(conn(:,3),2));
% find points inside
I = polygon_inside(X,Y,P{1}) ~= 0;
% is there a hole
if numel(P) >= 2    
    I = I & polygon_inside(X,Y,P{2}) == 0;
end
conn = conn(I,:);
end

function conn = clean_zero_area(coor,conn)
X1 = coor(conn(:,1),1);
X2 = coor(conn(:,2),1);
X3 = coor(conn(:,3),1);

Y1 = coor(conn(:,1),2);
Y2 = coor(conn(:,2),2);
Y3 = coor(conn(:,3),2);

A = X1.*Y2 - X2.*Y1 + X2.*Y3 - X3.*Y2 + X3.*Y1 - X1.*Y3;
conn = conn(abs(A)>0,:);
end

function conn = fix_upsidedown(coor,conn)
% get the nodal positions per element
X1 = coor(conn(:,1),1);
X2 = coor(conn(:,2),1);
X3 = coor(conn(:,3),1);

Y1 = coor(conn(:,1),2);
Y2 = coor(conn(:,2),2);
Y3 = coor(conn(:,3),2);

A = X1.*Y2 - X2.*Y1 + X2.*Y3 - X3.*Y2 + X3.*Y1 - X1.*Y3;
A = A < 0;

% flip the triangles
con = [1 3 2];
conn(A,:) = conn(A,con);
end


