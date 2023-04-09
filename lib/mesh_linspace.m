function x = mesh_linspace(a, b, N, c, type)
% x = MESH_LINSPACE(a, b, N), create a node vector from a to b with N nodes
% and the outside elements a factor sqrt(2) larger
%
% x = MESH_LINSPACE(a, b, N, c), instead of sqrt(2) define the factor
% manually
%
% x = MESH_LINSPACE(a, b, L, c, 'length'), instead of defining the number
% of nodes, the 3rd argument now is the element length for the internal
% elements.
%
% This function is a modified version of LINSPACE. The main advantage of
% such a node vector is that when combined with MESH_GEN_STRUCTURED, it
% will create larger elements at the boundaries compensating for their lack
% of robustness because these elements have less neighbors.
%
% See also, mesh_gen_structured

if (nargin <= 3) || isempty(c)
    c = sqrt(2);
end

if nargin <= 4
    type = 'number';
end

if strcmpi(type, 'number')

    if N < 2
        error('must use 2 nodes minimum')
    end

    % the element length of the external nodes
    L = 1 ./ (2 + (N - 3) ./ c);

    x = a + (b - a) * [0, linspace(L, 1 - L, N - 2), 1];

elseif strcmpi(type, 'length')

    % the length of the inside elements
    Li = N;

    % the length of the outside elements
    Lo = c * Li;

    % the length for the inside elements
    Lt = (b - a) - 2 * Lo;

    if Lt < 0
        % only one element fits
        x = [a, b];
        return
    elseif Lt < Li
        % only one extra node fits
        x = [a, 0.5 * (a + b), b];
        return
    end

    % how many fit inside
    Ni = floor(Lt ./ Li);

    % a vector of internal nodes
    xi = (0:Ni) * Li;

    % center the vector
    xi = xi - 0.5 * (xi(Ni + 1) + xi(1));

    % generate the vector
    x = [a, 0.5 * (a + b) + xi, b];

else
    error('unknown type')
end
