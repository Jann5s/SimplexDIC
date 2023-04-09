function cor = SimplexDIC_T4(f,g,init,coor,conn,varargin)
% cor = SimplexDIC_T4(f,g,init,coor,conn), correlate the image g on f using
% a T4 mesh as defined by coor and conn, starting with init as initial
% guess.
% init : [ux(:), uy(:), uz(:)]
% coor : [x(:), y(:), z(:)]
% conn : [n1(:), n2(:), n3(:), n4(:)]
%
% cor = SimplexDIC_T4(f,g,init,coor,conn,opt), where opt is a structure
% with optional fields:
% convcrit: convergence threshold (1e-4)
% maxit:    max number of iterations (20)
% blur:     apply image blur before dic, blur can be a list of blur radii
%           e.g. blur = [20, 10, 1, 0]
% mask:     enable individual pixels using a mask,
% method:   interpolation method (for the last blur step only)
% verbose:  0,1,2
% wantR:    (true) set the false to disable storing of R
% wantU:    (false) set the true to enable storing of Ux, Uy and Uz
% plotflag: (true) set to false to disable plotting

figstr = 'SimplexDIC_T4';

assert(ndims(f) == 3,'expecting a 3D matrix for f')
assert(ndims(g) == 3,'expecting a 3D matrix for g')
assert(all(size(f) == size(g)),'f and g must be the same size')

% set the correct types
if ~isa(f,'double')
    f = double(f);
end
if ~isa(g,'double')
    g = double(g);
end

% the image size
siz = size(f);
n = siz(1);
m = siz(2);
d = siz(3);

% defaults
opt.convcrit = 1e-4;
opt.maxit = 20;
opt.blur = 0;
opt.verbose = 2;
opt.LM = 1e-3; % levenberg marquardt
opt.tikhonov = 1e-9; % tikhonov
opt.mask = true(n,m);
opt.method = 'spline';
opt.plotflag = true;
opt.wantR = true;
opt.wantU = false;
opt.Hax = [];

% overwrite options
if (nargin == 6) && ~isempty(varargin{1}) && isstruct(varargin{1})
    fields = fieldnames(varargin{1});
    for k = 1:numel(fields)
        field = (fields{k});
        if ~isempty(field)
            opt.(field) = varargin{1}.(field);
        end
    end
end

% create simple voxel coordinates
x = 1:siz(2);
y = 1:siz(1);
z = 1:siz(3);
[X, Y, Z] = meshgrid(x,y,z);

% blur levels
blur = opt.blur;
Nb = numel(blur);

% number of nodes
Nn = size(coor,1);

% number of voxels
Np = n * m * d;

if opt.verbose >= 1
    fprintf('%s ------------------------\n',figstr);
end

% shape functions
if opt.verbose == 2
    fprintf('Generating the shapefunctions...\n')
end
[phi, E] = TetShapefunGridded(coor,conn,siz);

% combine the mesh and the mask
mask = (E ~= 0) & opt.mask;

% options
method = 'linear';

% DOF indices
Ix = 1:3:3*Nn;
Iy = 2:3:3*Nn;
Iz = 3:3:3*Nn;

% initial guess
a = zeros(3*Nn,1);
a(Ix) = init(:,1);
a(Iy) = init(:,2);
a(Iz) = init(:,3);

% tikhonov reference
aref = a;

if opt.verbose == 1
    fprintf('%4s, %3s, %10s, %10s\n','blur','it','res','da');
end
for bk = 1:Nb
    if opt.verbose == 2
        fprintf('Starting blur %d (%g px)...\n',bk,opt.blur(bk))
        fprintf('   %3s, %10s, %10s\n','it','res','da');
    end
    
    % blur the images
    if blur(bk) > 0
        fb = imgaussfilt3(f,blur(bk));
        gb = imgaussfilt3(g,blur(bk));
    else
        fb = f;
        gb = g;
    end
    
    % gradient
    [fx, fy, fz] = gradient(fb);
    
    % Least Squares matrix (Gauss)
    Lx = bsxfun(@times,phi(mask,:),fx(mask));
    Ly = bsxfun(@times,phi(mask,:),fy(mask));
    Lz = bsxfun(@times,phi(mask,:),fz(mask));
    
    % Hessian
    M = spalloc(3*Nn, 3*Nn, 3*3*3*3*Nn);
    M(Ix,Ix) = transpose(Lx)*Lx;
    M(Ix,Iy) = transpose(Lx)*Ly;
    M(Ix,Iz) = transpose(Lx)*Lz;
    M(Iy,Ix) = transpose(Ly)*Lx;
    M(Iy,Iy) = transpose(Ly)*Ly;
    M(Iy,Iz) = transpose(Ly)*Lz;
    M(Iz,Ix) = transpose(Lz)*Lx;
    M(Iz,Iy) = transpose(Lz)*Ly;
    M(Iz,Iz) = transpose(Lz)*Lz;
    
    % get the first eigenvector of M
    lambda = eigs(M,1);
    
    % regularization matrices
    Mlm  = lambda * opt.LM       * speye(3*Nn);
    Mtik = lambda * opt.tikhonov * speye(3*Nn);
    
    % add regularization
    M = M + Mlm + Mtik;
    
    if bk == Nb
        method = opt.method;
    end
    
    % prepare the interpolator
    G = griddedInterpolant(Y,X,Z,gb,method,'none');
    
    for it = 1:opt.maxit
        Ux = reshape(phi*a(Ix),n,m,d);
        Uy = reshape(phi*a(Iy),n,m,d);
        Uz = reshape(phi*a(Iz),n,m,d);
        
        gt = G(Y+Uy,X+Ux,Z+Uz);
        res = fb - gt;
        
        b = zeros(3*Nn,1);
        b(Ix,1) = transpose(Lx) * res(mask);
        b(Iy,1) = transpose(Ly) * res(mask);
        b(Iz,1) = transpose(Lz) * res(mask);
        
        % add tikhonov regularization
        b = b + lambda * opt.tikhonov * (a - aref);
        
        % solve the system
        da = M \ b;
        
        % update the dof
        a = a + da;
        
        res(~mask) = NaN;
        r = rms(res(mask));
        
        if opt.verbose == 2
            fprintf('   %3d, %10.3e, %10.3e\n',it,r,rms(da));
        end
        
        if rms(da) < opt.convcrit
            break
        end
    end
    
    if opt.verbose == 1
        fprintf('%4d, %3d, %10.3e, %10.3e\n',blur(bk),it,rms(res(Iroi)),rms(da));
    end
end

% final results
Ux = reshape(phi*a(Ix),n,m,d);
Uy = reshape(phi*a(Iy),n,m,d);
Uz = reshape(phi*a(Iz),n,m,d);
gt = G(Y+Uy,X+Ux,Z+Uz);
res = f - gt;
res(~mask) = NaN;
r = rms(res(mask));

cor.r = r;
cor.a = [a(Ix),a(Iy),a(Iz)];
if opt.wantR
    cor.res = single(res);
end
if opt.wantU
    cor.Ux = single(Ux);
    cor.Uy = single(Uy);
    cor.Uz = single(Uz);
end



