clear; close all;

% read the images
imagepath = 'images/virtimage3D';
f = raw_read(fullfile(imagepath,'virtimage3D_00_l_uint8_160x160x160.raw'));
g = raw_read(fullfile(imagepath,'virtimage3D_03_l_uint8_160x160x160.raw'));

% convert to floating point
f = double(f);
g = double(g);

% image size
[n, m, d] = size(f);

% adding noise
f = f + 1e-3 * randn(size(f));
g = g + 1e-3 * randn(size(g));

% create a mesh
% =======================================================

% Create a grid of points
mar = 0.1;
Nnodes = [4, 5, 6];
xn = m * linspace(mar, 1 - mar, Nnodes(1));
yn = n * linspace(mar, 1 - mar, Nnodes(2));
zn = d * linspace(mar, 1 - mar, Nnodes(3));

% grid the coordinates
[Xn, Yn, Zn] = meshgrid(xn, yn, zn);

% store as a list of coordinates
coor = [Xn(:), Yn(:), Zn(:)];

% Connect the nodes (delaunay doesn't like very regular grids unless you add Qz)
conn = delaunayn(coor,{'Qt','Qbb','Qc','Qx','Qz'});

% fix insideout elements
[I, conn] = TetMesh_Fix_insideout(coor, conn);

% compute the element size
V = TetVolume(coor,conn);

% Remove zero volume elements
I = V > 1e-9;
conn = conn(I,:);
V = V(I,:);

% extrapolate to the nodes (JN: sum(V) ~= sum(Vn), probalby a bug)
Vn = TetNodeExtrap(coor,conn,V);

% Plot the mesh
% -----------------------

% create T3 elements from the T4 faces for plotting
[faces, facemap] = TetMesh_Faces(conn);

figure;
ha(1) = subplot(1,2,1);
patch('Vertices',coor,'Faces',faces,'FaceVertexCData',facemap * V,'Edgecolor','k','FaceColor','flat','Marker','.', 'FaceAlpha', 0.5);
set(gca,'YDir','reverse');
colorbar
title('Element Size')
view(3)

ha(2) = subplot(1,2,2);
patch('Vertices',coor,'Faces',faces,'FaceVertexCData',Vn,'Edgecolor','k','FaceColor','interp','Marker','.', 'FaceAlpha', 0.5);
set(gca,'YDir','reverse');
colorbar
title('Element Size per node')
view(3)

drawnow

% Run T4 DIC:
% ==================================

% number of nodes
Nn = size(coor,1);

% initial guess
init = zeros(Nn,3);

% options
opt.blur = [3, 1, 0];
opt.wantU = true;

% DIC
cor = SimplexDIC_T4(f, g, init, coor, conn, opt);

% Plotting
ULim = 1.5*[-1,1];

figure;
subplot(1,3,1)
patch('Vertices',coor,'Faces',faces,'FaceVertexCData',cor.a(:,1),'EdgeColor','Interp','FaceColor','none');
title('Ux')
axis image
colorbar
caxis(ULim)
view(3)

subplot(1,3,2)
patch('Vertices',coor,'Faces',faces,'FaceVertexCData',cor.a(:,2),'EdgeColor','Interp','FaceColor','none');
title('Uy')
axis image
colorbar
caxis(ULim)
view(3)

subplot(1,3,3)
patch('Vertices',coor,'Faces',faces,'FaceVertexCData',cor.a(:,3),'EdgeColor','Interp','FaceColor','none');
title('Uz')
axis image
colorbar
caxis(ULim)
view(3)

