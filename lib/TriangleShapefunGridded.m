function [phi, varargout] = TriangleShapefunGridded(coor,conn,siz)
% phi = TriangleShapefunGridded(coor,conn,siz), compute T3 or T6
% shapefunctions (phi) from the mesh specified by coor and conn on an image
% of siz = [rows, columns]; 
%
% Alternatively, siz = [row1,row2,col1,col2]; indicating the first and last
% rows and columns respectively.
%
% [phi, E] = TriangleShapefunGridded(coor,conn,siz), also return and
% image with element id's for each point,
%
% [phi, E, phix, phiy] = TriangleShapefunGridded(coor,conn,siz), also
% return phix and phiy, the shapefunction derivatives in x and y
% respectively.
%
% This method uses an alogrithm that traverses all the pixel lines of an
% element to find the left and right intersections with the element edges.
% The shapefunctions are consequenlty computed by line interpolations from
% the nodes to these intersections and then from the intersections to the
% line of pixels that they span. This type of algorithm is much faster than
% the typical scattered method using inversion of the Jacobian, even more
% so when it is compiled before running it (see the matlab coder).


Ne = size(conn,1);
Nc = size(conn,2);
Nn = size(coor,1);

if numel(siz) == 2
    n = siz(1);
    m = siz(2);
elseif numel(siz) == 4
    n1 = siz(1);
    n2 = siz(2);
    m1 = siz(3);
    m2 = siz(4);
    n = n2-(n1-1);
    m = m2-(m1-1);
    coor(:,1) = coor(:,1)-(m1-1);
    coor(:,2) = coor(:,2)-(n1-1);
else
    error('unexpected input for the image size (3rd argument)')
end
Np = n*m;

quadratic = false;
if Nc == 6
    quadratic = true;
end

% element edge margin
mar = 1e-9;

% elementary triangle (barycentric) coordinates at the nodes
Le = eye(3);

% a list of triangle edges (node pairs)
edges = [1,2;2,3;3,1];

% Store all results as [n,m,3] matrix (3 or 6 shapefun per pixel)
I = repmat(reshape((1:Np).',n,m),1,1,Nc);
J = zeros(n,m,Nc);
S = zeros(n,m,Nc);

if nargout >= 2
    want_E = true;
    E = zeros(n,m,1);
else
    want_E = false;
end
if nargout >= 3
    want_dphi = true;
    Sx = zeros(n,m,Nc);
    Sy = zeros(n,m,Nc);
else
    want_dphi = false;
end

% initialize some temp vars
lj = zeros(1,3);
l1 = zeros(1,3);
l2 = zeros(1,3);
ind = zeros(2,1);

% for each element
for ke = 1:Ne
        
    % elemental connectivity
    con = conn(ke,:);
    
    % get the nodes for this element (only the T3 part)
    pos = coor(con(1:3),:);
    
    if want_dphi
        % element area (twice)
        A = det([ones(3,1),pos]);
        
        % T3 derivatives
        l1x = -(pos(3,2)-pos(2,2))./A;
        l1y =  (pos(3,1)-pos(2,1))./A;
        l2x = -(pos(1,2)-pos(3,2))./A;
        l2y =  (pos(1,1)-pos(3,1))./A;
        l3x = -(pos(2,2)-pos(1,2))./A;
        l3y =  (pos(2,1)-pos(1,1))./A;
    end
    
    % get the limits
    imin = ceil( min(pos(:,2)) - mar);
    imax = floor(max(pos(:,2)) + mar);
    
    % limit to the image height
    imin = max(imin,1);
    imax = min(imax,n);
    
    % for every line of pixels
    for i = imin:imax
        
        % get the left and right point
        for k = 1:3
            % get two vertices (relative to the first)
            o = pos(edges(k,1),:);
            p = pos(edges(k,2),:) - o;
            
            % compute the line coordinates
            t  = (i-o(2))./p(2);
            t(t<0-mar|t>1+mar) = NaN;
            
            % interpolate x on the line coordinates
            lj(k) = o(1) + p(1)*t(:);
            
            % interpolate the element coordinates
            l1(k) = Le(1,edges(k,2))*t + Le(1,edges(k,1))*(1-t);
            l2(k) = Le(2,edges(k,2))*t + Le(2,edges(k,1))*(1-t);
        end
        
        % get the left and right point (one edge is NaN)
        [j1, ind(1)] = min(lj,[],'omitnan');
        [j2, ind(2)] = max(lj,[],'omitnan');
        
        % create a list of x-indices (limited to the image width)
        j = max( ceil(j1-mar),1):min(floor(j2+mar),m);
        if isempty(j)
            continue
        end
        
        % interpolate the shapefunctions
        if any(j - j1 == 0)
            t = 0;
        else
            t = (j-j1)./(j2-j1);
        end
        l1i = l1(1,ind(2))*t + l1(1,ind(1))*(1-t);
        l2i = l2(1,ind(2))*t + l2(1,ind(1))*(1-t);
        l3i = 1 - (l1i+l2i);
                
        % Assemble in the final data structure
        % ------------------
        
        % store the shapefunctions
        if want_E
            E(i,j,1) = ke;
        end
        
        J(i,j,1) = con(1);
        J(i,j,2) = con(2);
        J(i,j,3) = con(3);
        if quadratic
            J(i,j,4) = con(4);
            J(i,j,5) = con(5);
            J(i,j,6) = con(6);

            S(i,j,1) = l1i.*(2*l1i-1);
            S(i,j,2) = l2i.*(2*l2i-1);
            S(i,j,3) = l3i.*(2*l3i-1);
            S(i,j,4) = 4*l1i.*l2i;
            S(i,j,5) = 4*l2i.*l3i;
            S(i,j,6) = 4*l3i.*l1i;
        else
            S(i,j,1) = l1i;
            S(i,j,2) = l2i;
            S(i,j,3) = l3i;
        end
        
        if want_dphi && quadratic
            Sx(i,j,1) = 4*l1i.*l1x - l1x;
            Sx(i,j,2) = 4*l2i.*l2x - l2x;
            Sx(i,j,3) = 4*l3i.*l3x - l3x;
            Sx(i,j,4) = 4*l1i.*l2x + 4*l2i*l1x;
            Sx(i,j,5) = 4*l2i.*l3x + 4*l3i*l2x;
            Sx(i,j,6) = 4*l3i.*l1x + 4*l1i*l3x;

            Sy(i,j,1) = 4*l1i.*l1y - l1y;
            Sy(i,j,2) = 4*l2i.*l2y - l2y;
            Sy(i,j,3) = 4*l3i.*l3y - l3y;
            Sy(i,j,4) = 4*l1i.*l2y + 4*l2i*l1y;
            Sy(i,j,5) = 4*l2i.*l3y + 4*l3i*l2y;
            Sy(i,j,6) = 4*l3i.*l1y + 4*l1i*l3y;
        elseif want_dphi && ~quadratic
            Sx(i,j,1) = l1x;
            Sx(i,j,2) = l2x;
            Sx(i,j,3) = l3x;

            Sy(i,j,1) = l1y;
            Sy(i,j,2) = l2y;
            Sy(i,j,3) = l3y;
        end
    end
end

% Assemble the matrix sparsely
Is = (S ~= 0);
phi = sparse(I(Is),J(Is),S(Is) ,Np,Nn);

if want_E
    varargout{1} = E;
end
if want_dphi
    varargout{2} = sparse(I(Is),J(Is),Sx(Is),Np,Nn);
    varargout{3} = sparse(I(Is),J(Is),Sy(Is),Np,Nn);
end


