function A = image_coarsegrain(A,sps)
% B = IMAGE_COARSEGRAIN(A), return a coarsegrained image where each group
% of 2x2 pixels is averaged to 1 superpixel (or 2x2x2 for 3D images). The
% algorithm always starts at the top most pixel (i.e. A(1,1)). Remember
% that this changes the pixel coordinates, it is advised to use the related
% function mesh_coarsegrain to correct the mesh accordingly.
%
% B = IMAGE_COARSEGRAIN(A,sps), create superpixels of size (sps x sps) or 
% (sps x sps x sps) for 3d images.
%
% B = IMAGE_COARSEGRAIN(A,[sps1,sps2]), create superpixels of size (sps1 x sps2)
% for 2d images.
% 
% B = IMAGE_COARSEGRAIN(A,[sps1,sps2,sps3]), create superpixels of size (sps1 x
% sps2 x sps3) for 3d images.
% 
% see also, mesh_coarsegrain
%
% Tags, image,

% 2015-08, New version which allows coarsegraining differently in each
% direction

% This algorithim superpixels the image one direction at the time. This way
% the image is already shrinking during the coarse graining making each
% extra direction cheaper.

if nargin < 2
  sps = 2;
end

% convert to doubles to preserve accuracy
cls = class(A);
A = double(A);

if numel(sps) == 1
    sps = repmat(sps,1,ndims(A));
end

spsOrig = sps;
if any(sps > size(A))
  sz1 = size(A);
  sps(sps > size(A)) = sz1(sps > size(A));
end

% perform the superpixeling in order of superpixel size
[~, ind] = sort(sps, 'descend');
for i = ind
  % "extra" regions are discarded if the image is not a multiple of the
  % superpixel size
  I = 1 : sps(i) : floor(size(A,i)/sps(i))*sps(i);

  % superpixel the image in the first direction
  sz2 = size(A);
  sz2(i) = length(I);
  B = zeros(sz2);
  for j = 1:sps(i)
    switch i
      case 1; B = B + A(I,:,:);
      case 2; B = B + A(:,I,:);
      case 3; B = B + A(:,:,I);
    end
    % go the the next subpixel in this row
    I = I + 1;
  end
  % overwrite A
  A = B / spsOrig(i);
end

% cast back to the original type
A = cast(A,cls);