function savepng(filename,varargin)
% savepng(filename), save the current figure as a file
%
% savepng(filename,H), save the figure with handle H as a file
%
% savepng(filename,[],dpi), save at dpi resolution
%
% Tags, save, plot
%

H = gcf;
r = 200;

if nargin == 2
    if ~isempty(varargin{1}) && ishandle(varargin{1})
        H = varargin{1};
    end
elseif nargin == 3
    if ~isempty(varargin{1}) && ishandle(varargin{1})
        H = varargin{1};
    else
        r = varargin{1};
    end
    
    if ~isempty(varargin{2}) && ishandle(varargin{2})
        H = varargin{2};
    else
        r = varargin{2};
    end
end
C = get(H,'Color');
set(H,'Color',[1 1 1]);


rstr = sprintf('-r%d',r);

% fix the extention, to be always .png (small case)
filename = [regexprep(filename,'.png$','','ignorecase') '.png'];

% get the original figure position (and size)
savepos = get(H,'Position');
% set the paper position to 1 inch per 100 pixels
set(H,'PaperUnits','inches','PaperPosition',savepos.*[0 0 1e-2 1e-2])
set(H,'PaperSize',savepos(3:4).*[1e-2 1e-2])

% save pdf
print(H,filename,'-dpng',rstr)

set(H,'Color',C);

