function [phi, E] = TetShapefunGridded(coor,conn,siz,varargin)
% [phi, E] = TetShapefunGridded(coor,conn,siz), compute T4
% shapefunctions (S) from the mesh specified by coor and conn on an image
% of siz = [rows, columns, slices]; E stores the element id for each voxel;
%
% [phi, E] = TetShapefunGridded(coor,conn,siz,margin) specify the
% margin inside of which pixels are still considered in the element,
% relative to the element size
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

n = siz(1);
m = siz(2);
d = siz(3);
Np = n*m*d;

if Nc ~= 4
    error('incorrect number of nodes per element')
end

% element edge margin
mar = 1e-9;
if (nargin == 4) && ~isempty(varargin{1})
    mar = varargin{1};
end

% a list of the 6 tetrahedron edges (node pairs)
edges = [1,2;2,3;3,1;1,4;2,4;3,4];

% sorting indices for 4-vertex polygon
vp2 = [1,2,1;2,1,3;3,3,2;4,4,4];
vp3 = vp2([2,3,4,1],:);
vp1 = vp2([4,1,2,3],:);

% Store all results as [n,m,4] matrix (4 shapefun per pixel + E)
S = zeros(n,m,d,Nc);
E = zeros(n,m,d,1);

I = repmat(reshape((1:Np).',n,m,d),1,1,1,Nc);
J = zeros(n,m,d,Nc);

% initialize some temp vars
li = zeros(1,6);
lj = zeros(1,6);
l1 = zeros(1,6);
l2 = zeros(1,6);
l3 = zeros(1,6);
ind = zeros(2,1);

lip = zeros(1,4);
l1p = zeros(1,4);
l2p = zeros(1,4);
l3p = zeros(1,4);

% for each element
for ke = 1:Ne
    
    % elemental connectivity
    con = conn(ke,:);
    
    % get the nodes for this element (only the T3 part)
    pos = coor(con(1:4),1:3);
    
    % enlarge the tetrahedron with the margin
    cen = mean(pos,1);
    Le = zeros(4,4);
    V = zeros(4,3);
    for p = 1:4
        dp = pos(p,:) - cen;
        dn = mar/norm(dp);
        dm = 3*dn/4;
        V(p,:) = pos(p,:) + dn*dp;
        for pp = 1:4
            if p == pp
                Le(p,pp) = 1+dm;
            else
                Le(p,pp) = -dm/3;
            end                
        end
    end
    
    % get the tet bounding box
    lim(1,:) = min(pos,[],1);
    lim(2,:) = max(pos,[],1);
    
    % get the limits
    kmin = ceil (lim(1,3));
    kmax = floor(lim(2,3));
    
    % limit to the image height
    kmin = max(kmin,1);
    kmax = min(kmax,d);
        
    % for every slice in z
    for k = kmin:kmax
        
        % intersect each tetrahedron edge with the z-plane        
        pcnt = 0;
        for e = 1:6
            e1 = edges(e,1);
            e2 = edges(e,2);
            
            if V(e2,3) == V(e1,3)
                continue
            end
            
            % compute the line coordinates
            t = (k-V(e1,3))./(V(e2,3)-V(e1,3));
                        
            if (t >= 0) && (t <= 1)
                pcnt = pcnt + 1;
                % interpolate x and y on the line coordinates
                li(1,pcnt) = Lerp(V(e1,2),V(e2,2),t);
                lj(1,pcnt) = Lerp(V(e1,1),V(e2,1),t);
                
                % interpolate the element coordinates
                l1(1,pcnt) = Lerp(Le(1,e1),Le(1,e2),t);
                l2(1,pcnt) = Lerp(Le(2,e1),Le(2,e2),t);
                l3(1,pcnt) = Lerp(Le(3,e1),Le(3,e2),t);
            end            
        end

        % no intersections
        if (pcnt == 0)
            continue
        end                
        
        % get the left and right point (one edge is NaN)
        [j1, ind1] = min(lj(1:pcnt),[],'omitnan');
        [j2, ind2] = max(lj(1:pcnt),[],'omitnan');
        
        % the relevant columns
        jmin = max( ceil(j1),1);
        jmax = min(floor(j2),m);
        
        if jmax < jmin
            % no columns inside this slice, skip
            continue
        end
               
        % if the polygon contains 4 points, sort to get a quad
        qsort = 1:4;
        if pcnt == 4
            % test if the polygon is sorted in order
            x1 = lj(vp1)-lj(vp3);
            x2 = lj(vp2)-lj(vp3);
            y1 = li(vp1)-li(vp3);
            y2 = li(vp2)-li(vp3);

            % cross product (z-component)
            cp = x1.*y2 - x2.*y1;
            
            % find the version where all cross-products are in the same
            % direction
            Itt = abs(sum(sign(cp)))==4;
            
            if nnz(Itt) == 0
                continue
            end
            qsort = vp2(:,Itt);
        end
        

        % create a loop-around vector for any number of points
        vp = [1:pcnt, 1];        
        
        % for every column of pixels
        for j = jmin:jmax
                        
            % find the end-points of each column
            cnt = 0;
            for p = 1:pcnt
                p1 = qsort(vp(p));
                p2 = qsort(vp(p+1));
                
                if lj(p2) == lj(p1)
                    continue
                end
                
                % compute the column coordinates
                t = (j-lj(p1))./(lj(p2)-lj(p1));
                                
                if (t >= 0) && (t <= 1)                    
                    cnt = cnt + 1;
                    lip(1,cnt) = Lerp(li(p1),li(p2),t);
                    l1p(1,cnt) = Lerp(l1(p1),l1(p2),t);
                    l2p(1,cnt) = Lerp(l2(p1),l2(p2),t);
                    l3p(1,cnt) = Lerp(l3(p1),l3(p2),t);
                end
            end
            
            % no intersections
            if (cnt == 0)
                continue
            end
            
            % get the top and bottom point (one edge is NaN)
            [i1, ind1] = min(lip(1:cnt),[],'omitnan');
            [i2, ind2] = max(lip(1:cnt),[],'omitnan');
            
            % the relevant columns
            imin = max( ceil(i1),1);
            imax = min(floor(i2),n);
            
            if imax < imin
                % no rows inside this slice, skip
                continue
            end
            
            % for a range of y-indices (limited to the image height)
            for i = imin:imax
                
                % interpolate the shapefunctions
                t = (i-i1)./(i2-i1);
                if (i1 == i2)
                    % fix for possible divide by zero
                    t = 0;
                end
                l1i = Lerp(l1p(ind1),l1p(ind2),t);
                l2i = Lerp(l2p(ind1),l2p(ind2),t);
                l3i = Lerp(l3p(ind1),l3p(ind2),t);
                l4i = 1 - (l1i+l2i+l3i);
                                
                % Assemble in the final data structure
                % ------------------
                S(i,j,k,1) = l1i;
                S(i,j,k,2) = l2i;
                S(i,j,k,3) = l3i;
                S(i,j,k,4) = l4i;
                E(i,j,k,1) = ke;

                J(i,j,k,1) = con(1);
                J(i,j,k,2) = con(2);
                J(i,j,k,3) = con(3);
                J(i,j,k,4) = con(4);

            end  % loop over i
        end  % loop over j
    end % loop over k
end % loop over ke

% Assemble the matrix sparsely
Is = (S ~= 0);
phi = sparse(I(Is),J(Is),S(Is) ,Np ,Nn);


function x = Lerp(x1,x2,t)
x = x1*(1-t) + x2*t;
