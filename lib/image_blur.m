function A = image_blur(A,R)
% B = image_blur(A, R) blurs the image A with a gaussian kernel of radius R
% 
% If R takes the form :
%  [R1, R2, R3]  the kernel is anisotropic (different radius in each dimension)
%  [R1, R2]      the kernel is anisotropic (different radius in each dimension)
%  [R]           the kernel is isotropic, equivalent to [R, R] or [R, R, R]
% As such, R(i) is applied over the i-th dimension of the image
% 
% Notes:
% - images are mirrored to deal with the edges.
% - internally, doubles are used, but B is returned as the same type as A.
%   When A is an integer, this may lead to intensity clipping.
%
% see also, image_grayscale
%
% Tags, image,
%

% Changelog: 
% 2016-10, AMQ
% 2015-03, JN

% convert to doubles to preserve accuracy
cls = class(A);
A = double(A);

if nargin < 2
	error('must specify a radius or radii');
end

if numel(R) == 1
    R = R*ones(ndims(A));
end

% create the kernel
% ------------------------
kernel = 'gaussian';

for i = 1:ndims(A)
    
    if strcmpi(kernel,'gaussian')
        % box size
        s = 2*ceil(4*R(i)+0.5)-1;
        x = (1:s)-mean(1:s);
        x = x(:);
        K = (1/(sqrt(2*pi)*R(i)))*exp(-(x.^2)./(2*R(i)^2));
    end
    % mirroring the image
    % ------------------------
    [am,an,ap] = size(A);
    km         = length(K);
    
    % mirror distance
    Nm = floor(km/2);
    
    % skip a direction if it has too few pixels
    if (Nm >= am) || (R(i) == 0)
        if ismatrix(A)
            % transpose the image and blur again
            A = A';
        elseif ndims(A) == 3
            % transpose the image and blur again
            A = permute(A,[2 3 1]);
        end
        continue
    end
    
    Amir = zeros(am+2*Nm,an,ap);
    
    % image rows and columns in mirrored image
    Im = Nm+1:Nm+am;
        
    % image rows and columns in the original image
    Jm = 1:am;
    
    % mirror parts in mirrored image
    Iml = Nm:-1:1; % left columns
    Imr = 2*Nm+am:-1:Nm+1+am; % right columns
    
    % mirror parts in original image
    Jml = 2:Nm+1;
    Jmr = (am:am+Nm-1)-Nm;

    % Create an image with mirrored edges
    Amir([Iml,Im,Imr],:,:) = A([Jml,Jm,Jmr],:,:);
    
    % Apply the blur using convolution
    % ------------------------

    if ismatrix(A)
        A = conv2(Amir,K,'valid');

        % transpose the image and blur again
        A = A';
    elseif ndims(A) == 3
        A = convn(Amir,K,'valid');

        % transpose the image and blur again
        A = permute(A,[2 3 1]);
    end
    
    
end

% Recasting
% ------------------------

% cast back to the original type
A = cast(A,cls);


