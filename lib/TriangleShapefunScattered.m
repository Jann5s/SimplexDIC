function [phi, varargout] = TriangleShapefunScattered(coor,conn,pos)
% phi = TriangleShapefunScattered(coor,conn,pos), compute T3 or T6
% shapefunctions (phi) from the mesh specified by coor and conn on the
% points specified by pos=[x(:),y(:)].
%
% [phi, E] = TriangleShapefunScattered(coor,conn,pos), also return and
% image with element id's for each point,
%
% [phi, E, phix, phiy] = TriangleShapefunScattered(coor,conn,pos), also
% return phix and phiy, the shapefunction derivatives in x and y
% respectively.
%
% This method uses inversion of the Jacobian to compute the shapefunction
% values for each point in a box surrounding the element. Consequently,
% using the shapefunction values, the points inside the element are
% identified and stored. This method allows non-structured grids for pos,
% but is much slower then a gridded algorithm.


Ne = size(conn,1);
Nn = size(coor,1);
Nc = size(conn,2);
Np = size(pos,1);

quadratic = false;
if Nc == 6
    quadratic = true;
end

% element edge margin
mar = 1e-9;

% prepare storage for 3 or 6 shapefuns per pixel (and indices in phi)
I  = repmat((1:Np).',1,Nc);
J  = zeros(Np,Nc);
S  = zeros(Np,Nc);

if nargout >= 2
    want_E = true;
    E  = zeros(Np,1);
else
    want_E = false;
end
if nargout >= 3
    want_dphi = true;
    Sx = zeros(Np,Nc);
    Sy = zeros(Np,Nc);
else
    want_dphi = false;
end
for ke = 1:Ne

    % elemental connectivity
    con = conn(ke,:);
    
    % get the nodes for this element
    npos = coor(con(1:3),:);
    
    % derivatives
    if want_dphi
        % element area (twice)
        A = det([ones(3,1),npos]);
        
        % T3 derivatives
        l1x = -(npos(3,2)-npos(2,2))./A;
        l1y =  (npos(3,1)-npos(2,1))./A;
        l2x = -(npos(1,2)-npos(3,2))./A;
        l2y =  (npos(1,1)-npos(3,1))./A;
        l3x = -(npos(2,2)-npos(1,2))./A;
        l3y =  (npos(2,1)-npos(1,1))./A;
    end
        
    % get a box
    lim(1,:) = min(npos,[],1) - mar;
    lim(2,:) = max(npos,[],1) + mar;
    
    Ibox = find( pos(:,1) > lim(1,1) & pos(:,1) < lim(2,1) &...
        pos(:,2) > lim(1,2) & pos(:,2) < lim(2,2));
    
    % compute the element local coordinates (barycentric coordinates)
    T = npos(1:2,:) - repmat(npos(3,:),2,1);
    L = bsxfun(@minus,pos(Ibox,:),npos(3,:)) / T;
    L(:,3) = 1 - sum(L,2);
    
    % test which coordinates are in the triangle
    In = all(L >= -mar,2) & all(L <= 1+mar,2);

    % store the shapefunctions
    if want_E
        E(Ibox(In),1) = ke;
    end
    J(Ibox(In),1) = con(1);
    J(Ibox(In),2) = con(2);
    J(Ibox(In),3) = con(3);
    if quadratic
        J(Ibox(In),4) = con(4);
        J(Ibox(In),5) = con(5);
        J(Ibox(In),6) = con(6);

        S(Ibox(In),1) = L(In,1).*(2*L(In,1)-1);
        S(Ibox(In),2) = L(In,2).*(2*L(In,2)-1);
        S(Ibox(In),3) = L(In,3).*(2*L(In,3)-1);
        S(Ibox(In),4) = 4*L(In,1).*L(In,2);
        S(Ibox(In),5) = 4*L(In,2).*L(In,3);
        S(Ibox(In),6) = 4*L(In,3).*L(In,1);
    else
        S(Ibox(In),1) = L(In,1);
        S(Ibox(In),2) = L(In,2);
        S(Ibox(In),3) = L(In,3);        
    end
    
    if want_dphi && quadratic
        Sx(Ibox(In),1) = 4*L(In,1).*l1x - l1x;
        Sx(Ibox(In),2) = 4*L(In,2).*l2x - l2x;
        Sx(Ibox(In),3) = 4*L(In,3).*l3x - l3x;
        Sx(Ibox(In),4) = 4*L(In,1).*l2x + 4*L(In,2)*l1x;
        Sx(Ibox(In),5) = 4*L(In,2).*l3x + 4*L(In,3)*l2x;
        Sx(Ibox(In),6) = 4*L(In,3).*l1x + 4*L(In,1)*l3x;
        
        Sy(Ibox(In),1) = 4*L(In,1).*l1y - l1y;
        Sy(Ibox(In),2) = 4*L(In,2).*l2y - l2y;
        Sy(Ibox(In),3) = 4*L(In,3).*l3y - l3y;
        Sy(Ibox(In),4) = 4*L(In,1).*l2y + 4*L(In,2)*l1y;
        Sy(Ibox(In),5) = 4*L(In,2).*l3y + 4*L(In,3)*l2y;
        Sy(Ibox(In),6) = 4*L(In,3).*l1y + 4*L(In,1)*l3y;
    elseif want_dphi && ~quadratic
        Sx(Ibox(In),1) = l1x;
        Sx(Ibox(In),2) = l2x;
        Sx(Ibox(In),3) = l3x;
        
        Sy(Ibox(In),1) = l1y;
        Sy(Ibox(In),2) = l2y;
        Sy(Ibox(In),3) = l3y;        
    end    
end

% sparsify the matrix
Is = S ~= 0;
phi  = sparse(I(Is),J(Is),S(Is) ,Np,Nn);

if want_E
    varargout{1} = E;
end
if want_dphi
    varargout{2} = sparse(I(Is),J(Is),Sx(Is),Np,Nn);
    varargout{3} = sparse(I(Is),J(Is),Sy(Is),Np,Nn);
end



