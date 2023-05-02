function [S, E] = TetShapefunScattered(coor,conn,pos,varargin)
% [S, E] = TetShapefunScattered(coor,conn,pos), compute T4
% shapefunctions (phi) from the mesh specified by coor and conn on the
% points specified by pos=[x(:),y(:)]. E stores the element id for each
% pixel;
%
% [S, E] = TetShapefunScattered(coor,conn,siz,margin) specify the
% margin inside of which pixels are still considered in the element
%
% Notes:
% - This method uses inversion of the Jacobian to compute the shapefunction
% values for each point in a box surrounding the element. Consequently,
% using the shapefunction values, the points inside the element are
% identified and stored. This method allows non-structured grids for pos,
% but is much slower then a gridded algorithm.
% -


Ne = size(conn,1);
Nc = size(conn,2);
Nn = size(coor,1);
Np = size(pos,1);

% only consider the first three coordinates
pos = pos(:,1:3);
coor = coor(:,1:3);

if Nc ~= 4
    error('incorrect number of nodes per element')
end

% element edge margin
mar = 1e-9;
if (nargin >= 4) && ~isempty(varargin{1})
    mar = varargin{1};
end

% prepare storage for 4 shapefuns per point and an element number
S  = zeros(Np,Nc);
E  = zeros(Np,1);
I  = repmat((1:Np).',1,Nc);
J  = zeros(Np,Nc);


% Method Elements:
% ======================================================================
% the outer loop goes over the elements, this is more efficient if the
% number of elements is low compared to the number of points

for ke = 1:Ne
    
    % elemental connectivity
    con = conn(ke,:);
    
    % get the nodes for this element
    V = coor(con(1:4),:);
    
    % Get the Bounding Box
    % -----------------------
    
    % get a box
    lim(1,:) = min(V,[],1) - mar;
    lim(2,:) = max(V,[],1) + mar;
    
    % find the points inside the box
    Ibox = find(...
        pos(:,1) > lim(1,1) & pos(:,1) < lim(2,1) &...
        pos(:,2) > lim(1,2) & pos(:,2) < lim(2,2) &...
        pos(:,3) > lim(1,3) & pos(:,3) < lim(2,3));
    
    % First find the points in the triangle, then compute SF
    % -----------------------
    
    % find which points are inside the tetrahedron
    Itet = TetInside(V,pos(Ibox,:),mar);
    
    % reduce the list of points
    Ibox = Ibox(Itet);
    
    if isempty(Ibox)
        continue
    end
    
    T = V(1:3,:) - repmat(V(4,:),3,1);
    b = bsxfun(@minus,pos(Ibox,:),V(4,:));
    L = b/T;
    
    % store the shapefunctions
    S(Ibox,1) = L(:,1);
    S(Ibox,2) = L(:,2);
    S(Ibox,3) = L(:,3);
    S(Ibox,4) = 1 - (L(:,1) + L(:,2) + L(:,3));
    E(Ibox,1) = ke;
    J(Ibox,1) = con(1);
    J(Ibox,2) = con(2);
    J(Ibox,3) = con(3);
    J(Ibox,4) = con(4);
end

% sparsify the matrix
Is = S ~= 0;
phi  = sparse(I(Is),J(Is),S(Is) ,Np,Nn);



