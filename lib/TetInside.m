function I = TetInside(V,P,varargin)
% I = TetInside(V,P), compute which points (P) are inside the tetrahedron
% defined by the vertices (V).
% V = [Vx(:), Vy(:), Vz(:)], having 4 rows
% P = [Px(:), Py(:), Pz(:)], having Np rows
%
% I = TetInside(V,P,margin), consider points inside if they are at
% a distance margin outside of the tetrahedron, by default margin = 1e-9

% element edge margin
mar = 1e-9;
if (nargin == 3) && ~isempty(varargin{1})
    mar = varargin{1};
end

% number of vertices and points
Nv = size(V,1);
Np = size(P,1);

assert(Nv == 4,'number of vertices must be 4')

% tetrahedron centroid
% C = mean(V,1);

% positive face definition
faces = [3,2,4 ; 1,3,4 ; 2,1,4 ; 1,2,3] ; 

I = true(Np,1);
% for each face
for k = 1:Nv
    % the face vector indices
    f = faces(k,:);
    
    %     % Median length ( (4/3)*distance to the centroid )
    %     M = (4/3)*sum( (V(k,:)-C).^2 ,2)^0.5;
    
    % get the surface normal
    N = cross(V(f(3),:)-V(f(2),:),V(f(1),:)-V(f(2),:));
    N = N./norm(N);
    
    % the signed distance (dot product)
    D = sum(bsxfun(@times,N,bsxfun(@minus,P,V(f(2),:))),2);
            
    % set if inside (plus margin)
    I = I & D > -mar;
end

