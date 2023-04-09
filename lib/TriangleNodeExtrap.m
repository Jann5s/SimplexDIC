function B = TriangleNodeExtrap(coor,conn,A)

revertcell = false;
if ~iscell(A)
    A = {A};
    revertcell = true;
end
Ninc = numel(A);

Nn = size(coor,1);
Nc = size(A{1},2);

% Extrapolate to the nodes and store the increment data
for inc = 1:Ninc
    B{inc} = zeros(Nn,Nc);
    NC = zeros(Nn,1);
    for k = 1:3
        B{inc}(conn(:,k),:) = B{inc}(conn(:,k),:) + A{inc};
        NC(conn(:,k)) = NC(conn(:,k)) + 1;
    end
    NC(NC == 0) = 1;
    B{inc} = bsxfun(@rdivide,B{inc},NC);
end

if revertcell
    B = B{1};
end
