function cor = SimplexDIC_T3(f,g,init,coor,conn,varargin)
% cor = SimplexDIC_T3(f,g,init,coor,conn), correlate the image g on f using
% a T3 mesh as defined by coor and conn, starting with init as initial
% guess. 
% init : [ux(:), uy(:)]
% coor : [x(:), y(:)]
% conn : [n1(:), n2(:), n3(:)]
%
% cor = SimplexDIC_T3(f,g,init,coor,conn,opt), where opt is a structure
% with optional fields: 
% convcrit: convergence threshold (1e-4)
% maxit:    max number of iterations (20)
% blur:     apply image blur before dic, blur can be a list of blur radii 
%           e.g. blur = [20, 10, 1, 0]
% mask:     enable individual pixels using a mask,
% method:   interpolation method (for the last blur step only)
% verbose:  0,1,2
% wantR:    (true) set the false to disable storing of R
% wantU:    (false) set the true to enable storing of Ux and Uy
% wantE:    (false) set the true to enable storing of Exx, Eyy and Exy
% plotflag: (true) set to false to disable plotting
% CLim:     the color range used for plotting

if ~ismatrix(f) || ~ismatrix(g)
    error('expecting f and g is 2D matrices')
end

% set the correct types
if ~isa(f,'double')
    f = double(f);
end
if ~isa(g,'double')
    g = double(g);
end
if ~all(size(f) == size(g)) 
    error('f and g must be the same size')
end

% the image size
siz = size(f);
n = siz(1);
m = siz(2);

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
opt.wantE = false;
opt.CLim = 0.1*quantile(f(:),0.99) * [-1, 1];
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

% create simple coordinates
x = 1:siz(2);
y = 1:siz(1);
[X, Y] = meshgrid(x,y);

% blur levels
blur = opt.blur;
Nb = numel(blur);

% number of nodes
Nn = size(coor,1);

% number of pixels
Np = n * m;

% quadratic elements?
if size(conn,2) == 3
    figstr = 'SimplexDIC_T3';
    Iface = 1:3;
elseif size(conn,2) == 6
    figstr = 'SimplexDIC_T6';
    Iface = [1, 4, 2, 5, 3, 6];
else
    error('connectivity (conn) should have 3 or 6 columns')
end

if opt.verbose >= 1
    fprintf('%s ------------------------\n',figstr);
end

% shape functions
if opt.verbose == 2
    fprintf('Generating the shapefunctions...\n')
end
if opt.wantE
    [phi, E, phix, phiy] = TriangleShapefunGridded(coor,conn,siz);
else
    [phi, E] = TriangleShapefunGridded(coor,conn,siz);
end

% combine the mesh and the mask
mask = (E ~= 0) & opt.mask;

% options
method = 'linear';

% DOF indices
Ix = 1:2:2*Nn;
Iy = 2:2:2*Nn;

% initial guess
a = zeros(2*Nn,1);
a(Ix) = init(:,1);
a(Iy) = init(:,2);

% tikhonov reference
aref = a;

if opt.plotflag
    if isempty(opt.Hax) || ~ishandle(opt.Hax)
        figure('Name',figstr);
        colormap(colorbrewer(51,'RdBu'));
        ha = axes('Position',[0.05,0.05,0.9,0.9],'DataAspectRatio',[1 1 1],'XColor','none','YColor','none','Xlim',[1, m],'Ylim',[1, n],'NextPlot','Add','YDir','reverse','Color','k','Clipping','off');
        hi = imagesc(f-g);
        hp = patch('Vertices',coor,'Faces',conn(:,Iface),'EdgeColor','y','FaceColor','none','EdgeAlpha',0.2);
        set(hi,'AlphaData',1 - 0.3*(~mask));
        set(ha,'CLim',opt.CLim);
        colorbar
    else
        ha = opt.Hax;
        set(get(ha,'Parent'),'Name',figstr);
        hi = findobj(ha,'Type','image');
        hp = findobj(ha,'Type','image','Tag','Mesh');
    end
end

if opt.verbose == 1
    fprintf('%4s, %3s, %10s, %10s\n','blur','it','res','da');
end
for bk = 1:Nb
    if opt.verbose == 2
        fprintf('Starting blur %d (%g px)...\n',bk,opt.blur(bk))
        fprintf('   %3s, %10s, %10s\n','it','res','da');
    end
    
    fb = image_blur(f,blur(bk));
    gb = image_blur(g,blur(bk));
    
    % gradient
    [fx, fy] = gradient(fb);
    
    % Least Squares matrix (Gauss)
    Lx = bsxfun(@times,phi(mask,:),fx(mask));
    Ly = bsxfun(@times,phi(mask,:),fy(mask));
    
    % Hessian
    M = zeros(2*Nn,2*Nn);
    M(Ix,Ix) = transpose(Lx)*Lx;
    M(Ix,Iy) = transpose(Lx)*Ly;
    M(Iy,Ix) = transpose(Ly)*Lx;
    M(Iy,Iy) = transpose(Ly)*Ly;
    
    % get the first eigenvector of M
    lambda = eigs(M,1);
    
    % regularization matrices
    Mlm  = lambda * opt.LM       * speye(2*Nn,2*Nn);
    Mtik = lambda * opt.tikhonov * speye(2*Nn,2*Nn);
    
    % add regularization
    M = M + Mlm + Mtik;

    if bk == Nb
        method = opt.method;
    end
    
    % prepare the interpolator
    G = griddedInterpolant(Y,X,gb,method,'none');    
    
    for it = 1:opt.maxit
        Ux = reshape(phi*a(Ix),n,m);
        Uy = reshape(phi*a(Iy),n,m);
        
        gt = G(Y+Uy,X+Ux);
        res = fb - gt;
        
        b = zeros(2*Nn,1);
        b(Ix,1) = transpose(Lx) * res(mask);
        b(Iy,1) = transpose(Ly) * res(mask);
        
        % add tikhonov regularizatoin
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
        if opt.plotflag && ishandle(hi)
            set(hi,'CData',res);
            set(hp,'Vertices',coor + [a(Ix), a(Iy)]);
            drawnow;
        end
        
        if rms(da) < opt.convcrit
            break
        end
    end
    
    if opt.verbose == 1
        fprintf('%4d, %3d, %10.3e, %10.3e\n',blur(bk),it,rms(res(Iroi)),rms(da));
    end
    if opt.plotflag && ishandle(hi)
        set(hi,'CData',res);
        drawnow;
    end
    
end

% final results
Ux = reshape(phi*a(Ix),n,m);
Uy = reshape(phi*a(Iy),n,m);
gt = G(Y+Uy,X+Ux);
res = f - gt;
res(~mask) = NaN;
r = rms(res(mask));

if opt.plotflag
    if ishandle(hi(1))
        set(hi(1),'CData',res);
    end
    drawnow;
end

cor.r = r;
cor.a = [a(Ix),a(Iy)];
if opt.wantR
    cor.res = single(res);
end
if opt.wantU
    cor.Ux = single(Ux);
    cor.Uy = single(Uy);
end
if opt.wantE
    cor.Exx = single(reshape(phix*a(Ix),n,m));
    cor.Eyy = single(reshape(phiy*a(Iy),n,m));
    cor.Exy = single(reshape(phix*a(Iy),n,m));
end



