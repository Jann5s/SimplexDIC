function Ee = TriangleStrain(coor,conn,U,varargin)
% E = TriangleStrain(coor,conn,U), compute the strain for each T3 element
% defined by coor and conn using the displacements in U.
%
% E = TriangleStrain(coor,conn,U,type), specify a strain definition,
% options are {'small','GL','DefGrad'}, where the 'small' is the default
%
% Tensor components are always sorted:
% [ E11, E22, E12, E21 ]
%
% This function also accepts T6 elements, but in that case only the T3
% nodes are used to compute the strain, this is not what is normaly done to
% compute strains for T6 elements, use with caution.

% types = {'small','GL','DefGrad'};
straindef = 'small';
if nargin == 4;
    straindef = varargin{1};
end

Nn = size(coor,1);
Ne = size(conn,1);
if size(conn,2) > 3
    warning('These elements have more than 3 nodes, extra nodes are ignored');
end

revertcell = false;
if ~iscell(U)
    U = {U};
    revertcell = true;
end
Ninc = numel(U);

dN = [-1,-(1/sqrt(3)); 1,-(1/sqrt(3)); 0,2*(1/sqrt(3))];


% node coordinates (cancel the 3rd node)
X1 = coor(conn(:,1),:) - coor(conn(:,3),:);
X2 = coor(conn(:,2),:) - coor(conn(:,3),:);

% Jacobian % [J11, J22, J12, J21];
J = zeros(Ne,4);
J(:,1) = dN(1,1)*X1(:,1) + dN(2,1)*X2(:,1);
J(:,2) = dN(1,2)*X1(:,2) + dN(2,2)*X2(:,2);
J(:,3) = dN(1,1)*X1(:,2) + dN(2,1)*X2(:,2);
J(:,4) = dN(1,2)*X1(:,1) + dN(2,2)*X2(:,1);

% transpose and invert
invJt = inv2D(J(:,[1,2,4,3]));

% element derivative
dNe = [invJt(:,1)*dN(1:2,1)',invJt(:,3)*dN(1:2,1)'] + [invJt(:,4)*dN(1:2,2)',invJt(:,2)*dN(1:2,2)'];

% initiate some matrices
E = cell(Ninc,1);

% loop over the displacement increments
for inc = 1:Ninc
    
    % get the displacements (cancel the 3rd node)
    U1 = U{inc}(conn(:,1),:) - U{inc}(conn(:,3),:);
    U2 = U{inc}(conn(:,2),:) - U{inc}(conn(:,3),:);
    
    % DefGrad tensor
    F = dNe(:,[1,3,1,3]).*U1(:,[1,2,2,1]) + dNe(:,[2,4,2,4]).*U2(:,[1,2,2,1]);
    F(:,1:2) = F(:,1:2) + 1;
    
    % Strain
    if strcmpi(straindef,'DefGrad')
        Ee(:,[1,2,4,5]) = F;
    elseif strcmpi(straindef,'small')
        Ee(:,1) = F(:,1) - 1;
        Ee(:,2) = F(:,2) - 1;
        Ee(:,3) = 0.5*(F(:,3) + F(:,4));
        Ee(:,4) = 0.5*(F(:,4) + F(:,3));
    elseif strcmpi(straindef,'GL')
        Ee(:,1) = 0.5*(F(:,1).*F(:,1) + F(:,3).*F(:,3) - 1);
        Ee(:,2) = 0.5*(F(:,4).*F(:,4) + F(:,2).*F(:,2) - 1);
        Ee(:,3) = 0.5*(F(:,1).*F(:,4) + F(:,3).*F(:,2));
        Ee(:,4) = 0.5*(F(:,3).*F(:,1) + F(:,2).*F(:,4));
    end
    
    E{inc} = Ee;
end

if revertcell
    E = E{1};
end


function C = inv2D(A)
% row swap index vectors
swap12 = [4, 3, 2, 1];

% initiate the inverse
N = size(A,1);
C = [ones(N,2),zeros(N,2)];

% find the pivot
[vmax,imax] = max(abs(A(:,[1,4])),[],2);
if any(vmax == 0)
    warning('singular matrices detected')
end

% swap the rows
C(imax == 2,:) = C(imax == 2,swap12);
A(imax == 2,:) = A(imax == 2,swap12);

% modifying row 2 with row 1
C(:,[4,2]) = C(:,[4,2]) - C(:,[1,3]) .* repmat(A(:,4)./A(:,1),1,2);
A(:,[4,2]) = A(:,[4,2]) - A(:,[1,3]) .* repmat(A(:,4)./A(:,1),1,2);

% modifying row 1 with row 2
C(:,[1,3]) = C(:,[1,3]) - C(:,[4,2]) .* repmat(A(:,3)./A(:,2),1,2);
A(:,[1,3]) = A(:,[1,3]) - A(:,[4,2]) .* repmat(A(:,3)./A(:,2),1,2);

% scaling each row
C(:,[1,3]) = C(:,[1,3]) .* repmat(1./A(:,1),1,2);
C(:,[4,2]) = C(:,[4,2]) .* repmat(1./A(:,2),1,2);

function C = eig2D(A)
% solving the characteristic polynomial
b = -A(:,1)-A(:,2);
c =  A(:,1).*A(:,2) - A(:,3).*A(:,4);
d = b.^2 - 4*c;

% fix if near zero values are negative
d(d < eps) = 0;
d = sqrt(d);
eig1 = (-b + d)./2;
eig2 = (-b - d)./2;
C = [max(eig1,eig2),min(eig1,eig2)];

