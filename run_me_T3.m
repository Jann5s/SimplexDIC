clear; close all;

% add some the SimplexDIC library functions to the matlab path
addpath('lib')

% read the images
imagepath = 'images/virtimage2D';
f = imread(fullfile(imagepath,'virtimage_00.png'));
g = imread(fullfile(imagepath,'virtimage_03.png'));

% generate a color image (RGB) for later plotting
frgb = repmat(f,1,1,3);

% image size
[n, m] = size(f);

% scaling to [0, 1]
fmax = single(intmax(class(f)));
f = single(f) ./ fmax;
g = single(g) ./ fmax;

% adding noise
f = f + 1e-3 * randn(size(f));
g = g + 1e-3 * randn(size(g));

% create a mesh
% =======================================================

% generate a Hexagonal pattern mesh using TriangleMesh

% origin (top-left)
P = 0.1 * [m, n];

% element edge length
L = 50;

% element height
H = (0.5*sqrt(3)*L);

% number of elements
N(1) = floor( (0.8*n) / H );
N(2) = floor( (0.8*m) / L );

% generate the mesh
[coor3, conn3] = TriangleMesh(P,L,N);

% compute the element size
A = TriangleArea(coor3,conn3);

% extrapolate to the nodes
An = TriangleNodeExtrap(coor3,conn3,A);

% Plot the mesh
% -----------------------

figure;
ha(1) = subplot(1,2,1);
patch('Vertices',coor3,'Faces',conn3,'FaceVertexCData',A,'Edgecolor','k','FaceColor','flat','Marker','.');
set(gca,'YDir','reverse');
colorbar
title('Element Size')

ha(2) = subplot(1,2,2);
patch('Vertices',coor3,'Faces',conn3,'FaceVertexCData',An,'Edgecolor','k','FaceColor','interp','Marker','.');
set(gca,'YDir','reverse');
colorbar
title('Element Size per node')

drawnow


% Run T3 DIC:
% ==================================

% number of nodes
Nn = size(coor3,1);

% initial guess
init = zeros(Nn, 2);

% options
opt3.blur = [3, 1, 0];
opt3.CLim = 0.01 * [-1, 1];
opt3.wantU = true;
opt3.wantE = true;

% DIC
cor3 = SimplexDIC_T3(f, g, init, coor3, conn3, opt3);

% Plotting
% ---------------------------
ULim = 0.5*[-1,1];

figure;
subplot(2,2,1)
imagesc(cor3.Ux);
title('Ux')
axis image
colorbar
caxis(ULim)

subplot(2,2,2)
imagesc(cor3.Uy);
title('Uy')
axis image
colorbar
caxis(ULim)

subplot(2,2,3)
imagesc(frgb);
patch('Vertices',coor3,'Faces',conn3,'FaceVertexCData',cor3.a(:,1),'EdgeColor','Interp','FaceColor','none');
title('Ux')
axis image
colorbar
caxis(ULim)

subplot(2,2,4)
imagesc(frgb);
patch('Vertices',coor3,'Faces',conn3,'FaceVertexCData',cor3.a(:,2),'EdgeColor','Interp','FaceColor','none');
title('Uy')
axis image
colorbar
caxis(ULim)



