function cor = SimplexDIC_T6(f,g,init,coor,conn,varargin)
% cor = SimplexDIC_T6(f,g,init,coor,conn), correlate the image g on f using
% a T6 mesh as defined by coor and conn, starting with init as initial
% guess.
% init : [ux(:), uy(:)]
% coor : [x(:), y(:)]
% conn : [n1(:), n2(:), n3(:), n4(:), n5(:), n6(:)]
%
% Note, this code works with a special form of T6 elements where the
% quadratic nodes are assumed to be in the center of straight element
% edges.
%
% cor = SimplexDIC_T6(f,g,init,coor,conn,opt), where opt is a structure
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

cor = SimplexDIC_T3(f,g,init,coor,conn,varargin{:});


