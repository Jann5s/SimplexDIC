function [coor, conn] = TriangleMesh(P,L,N)
% [coor, conn] = TriangleMesh(P,L,N), generate a structured triangle mesh
% with N(1) rows and N2 columns of equilateral triangles with edge length L
% with P=[px,py] as the bottom left corner.
%
% coor are the nodal coordinates and conn the connectivity, to plot the
% mesh use for example:
% patch('Vertices',coor,'Faces',conn);

% prepare a Q4 grid for the RR and RC cases
nn = N([1, 2])+1;
ne = nn - 1;
Nn = prod(nn);
Ne = prod(ne);

% node index matrix
In = reshape(1:Nn,nn);

% create a Q4 grid with T3s
Ie = reshape(1:Ne,ne);

P = P(:).';

% create the T3 in the rectangular grid, sectioning the Q4s on opposing
% diagonals for the odd rows
conn = zeros(2*Ne,3);
for ki = 1:ne(1)
    for kj = 1:ne(2)
        con = In(ki:ki+1,kj:kj+1);
        ke1 = (Ie(ki,kj)-1)*2 + 1;
        ke2 = (Ie(ki,kj)-1)*2 + 2;
        
        if mod(ki,2) == 1
            % the 'upper' triangle
            conn(ke1,:) = con([1 3 2]);
            % the 'lower' triangle
            conn(ke2,:) = con([4 2 3]);
        else
            % the 'upper' triangle
            conn(ke2,:) = con([1 4 2]);
            % the 'lower' triangle
            conn(ke1,:) = con([1 3 4]);
        end
    end
end

% generate the node grid
x = L*(0:nn(2)-1) + P(1);
y = 0.5*sqrt(3)*L*(0:nn(1)-1) + P(2);
[X, Y] = meshgrid(x,y);

% push all odd rows to the right
I = 2:2:nn(1);
X(I,:) = X(I,:) + 0.5*L;

% store as coordinates
coor = [X(:),Y(:)];