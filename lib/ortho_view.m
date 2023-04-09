function varargout = ortho_view(A, varargin)
% ORTHO_VIEW(A), plot three cross-sections of the 3D matrix A. The
% presented images are interactive. Clicking anywhere in one of the three
% views will change which slice is shown in the other two depending on
% where is clicked.
%
% H = ORTHO_VIEW(A), return the figure handle of the generated figure
%
% Arrows: left/right move the current view by one, up/down more the current
%         view by 10.
% H     : hide/show the settings panel
%
% Options:
% H = ORTHO_VIEW(A, 'name', value) provide options using inputParser
%
% Available options are:
% ind       : [y, x, z] specify which slices to show initially, if empty,
%             the central slices are chosen
% handle    : provide an existing ortho_view figure and update it using
%             the new image A
% clim      : [cmin, cmax], enforce color limits, sets climmode to Manual
% climmode  : automatic color limits, a string with one of the names listed
%             in the <CLim mode> setting in the settings panel
% cmap      : color map, a string with one of the names listed in the
%             <CMap> setting in the settings panel
% flatten   : flatten the 3D image to three 2D images by applying a flatten
%             function, listed in the <Flatten mode>
% HidePanel : true/false, hide the settings panel
% Labels    : {'y','x','z'}, provide cellstring array of size 3 for the
%             labels in the order of the matrix indices
% FigHeight : [fh] the height of the figure in pixels, if empty 70% of the
%             monitor height is used
%
% Keyboard shortcuts:
% Space : autoplay all the slices of the current active view (the one with
%         the mouse over it) starting and ending at the current slice.
% Left/Right : move by one slice
% Up/Down : move by 10 slices
% H : hide the options panel
% M, L, F : cycle the colorMap, colorLim or Flatten modes, use shift to
% reverse cycle
%
% About the colorbar:
% The colorbar is always limited by the extremes of the color limits and
% the data values. Values that are outside of the color limits saturate to
% the extreme colors of the colormap. The color limits are also indicated
% with vertical dashed lines. The histogram of the data is plotted with a
% solid line. Moreover, the quantiles [0.001, 0.01, 0.1, 0.9, 0.99, 0.999]
% are shown as dots and the quantiles [0, 0.5, 1] (i.e. min, median, max)
% are shown as diamonds in the colorbar. The only exception to this is when
% the image contains only one color, in that case the colorbar will be
% centered around this value with 0.01 lower and upper limits.
%
% See also, stack_view, slice_plot

%#ok<*AGROW>
% Suppress warning for the entire file: "The variable 'xxx' appears to
% change size on every loop iteration. Consider preallocating for speed."
% AMQ: It would be nice not to do this...

assert(nargin > 0, 'ortho_view requires at least one argument: a 3D image')

% define a linear interpolation function
lerp = @(x,x1,x2) bsxfun(@times,(1-x),x1) + bsxfun(@times,x,x2);

Asiz = size(A);

if numel(Asiz) ~= 3
    error('expecting a 3D matrix')
end

% Compute the histogram
if isinteger(A)
    [Ahist, Ahistedges] = histcounts(A,'Normalization','count','BinMethod','integers');
else
    [Ahist, Ahistedges] = histcounts(A,'Normalization','count','BinMethod','auto');
end

% normalize the histogram
Ahist = Ahist ./ sum(Ahist);

% get the hist values
Ahistval = 0.5 * (Ahistedges(1:end-1) + Ahistedges(2:end));

Quantiles = [0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999];
if numel(Ahist) == 1
    % for when A just has one value
    Aqtl = Ahistval * ones(size(Quantiles));

else

    % Compute the quantiles using the histogram
    Q = cumsum(Ahist);
    Aqtl = zeros(size(Quantiles));
    for k = 1:numel(Quantiles)
        I2 = find(Q >= Quantiles(k),1);

        % force I2 to be within Q
        if isempty(I2)
            I2 = numel(Q);
        end

        % force I2 to be at least 2 and max numel(Q)
        I2 = max(I2, 2);

        % get the interpolation value
        x2 = Q(I2);
        x1 = Q(I2-1);
        x = (Quantiles(k) - x1) / (x2 - x1);

        % interpolate
        Aqtl(k) = lerp(x,Ahistedges(I2-1),Ahistedges(I2));
    end
end

% add the image limits to the quantils
Quantiles = [0, Quantiles, 1];
Aqtl = [double(min(A(:))), Aqtl, double(max(A(:)))];

% scale the histogram to [0, 1];
Ahist = Ahist ./ max(Ahist);

% Define the colormodes
climmodes{1, 1} = 'Auto 100%';
climmodes{2, 1} = 'Auto 99.9%';
climmodes{3, 1} = 'Auto 99.0%';
climmodes{4, 1} = 'Auto 90.0%';
climmodes{5, 1} = 'Centered 100%';
climmodes{6, 1} = 'Centered 99.9%';
climmodes{7, 1} = 'Centered 99.0%';
climmodes{8, 1} = 'Centered 90.0%';
climmodes{9, 1} = 'Manual';

% Define the colormaps
cmaps{ 1, 1} = 'Grays';
cmaps{ 2, 1} = 'Spectral';
cmaps{ 3, 1} = 'RdBu';
cmaps{ 4, 1} = 'RdGy';
cmaps{ 5, 1} = 'RdYlBu';
cmaps{ 6, 1} = 'YlGnBu';
cmaps{ 7, 1} = 'PuBuGn';
cmaps{ 8, 1} = 'YlOrBr';
cmaps{ 9, 1} = 'YlOrRd';
cmaps{10, 1} = 'Paired';
Ncolors = 51;

% a list of contrasting line colors
cmap_linecol{1, 1} = 'b';
cmap_linecol{2, 1} = 'k';
cmap_linecol{3, 1} = 'k';
cmap_linecol{4, 1} = 'k';
cmap_linecol{5, 1} = 'k';
cmap_linecol{6, 1} = 'k';
cmap_linecol{7, 1} = 'b';
cmap_linecol{8, 1} = 'k';
cmap_linecol{9, 1} = 'k';
cmap_linecol{10, 1} = 'k';

flattenfun{1, 1} = 'None';
flattenfun{2, 1} = 'Max';
flattenfun{3, 1} = 'Min';
flattenfun{4, 1} = 'Median';
flattenfun{5, 1} = 'Mean';
flattenfun{6, 1} = 'RMS';
flattenfun{7, 1} = 'STD';

Aflat{1} = zeros(Asiz([1,2]),'like',A);
Aflat{2} = zeros(Asiz([2,3]),'like',A);
Aflat{3} = zeros(Asiz([3,1]),'like',A);


% ====================
prs = inputParser;
addParameter(prs, 'ind', [], @isnumeric);
addParameter(prs, 'handle', [], @ishandle);
addParameter(prs, 'clim', [], @isnumeric);
addParameter(prs, 'climmode', climmodes{1}, @(x) any(ismember(x, climmodes)));
addParameter(prs, 'cmap', cmaps{1}, @(x) any(ismember(x, cmaps)));
addParameter(prs, 'flatten', flattenfun{1}, @(x) any(ismember(x, flattenfun)));
addParameter(prs, 'HidePanel', false, @islogical);
addParameter(prs, 'Labels', {'y', 'x', 'z'}, @iscellstr);
addParameter(prs, 'FigHeight', [], @isnumeric);
parse(prs, varargin{:});
opt = prs.Results;

Hf = [];

if ~isempty(opt.handle) && (numel(opt.handle) == 1) && ishandle(opt.handle)
    Hf = opt.handle;

    if ~strcmpi(Hf.Tag, 'ortho_view')
        error('The provided figure must be a ortho_view figure')
    end

    Hf = opt.handle;
end

ind = round(0.5 * (1 + Asiz));

if ~isempty(opt.ind) && (numel(opt.ind) == 3)
    ind = opt.ind;
end

OptionsVisible = 'on';

if opt.HidePanel
    OptionsVisible = 'off';
end

% find the clim mode
climmode = find(ismember(climmodes, opt.climmode));

if isempty(climmode)
    error('unknown clim mode %s', opt.climmode);
end

if ~isempty(opt.clim)
    climmode = 9;
end

% find the colormap
cmap = find(ismember(cmaps, opt.cmap));

if isempty(climmode)
    error('unknown cmap %s', opt.cmap);
end

% which indices are represented per axes
axind(1, :) = [2, 3]; % u,v
axind(2, :) = [2, 1];
axind(3, :) = [3, 1];

% which index is not represented per axes
axother = [1, 3, 2];

aspect(1) = Asiz(2) ./ Asiz(3);
aspect(2) = Asiz(2) ./ Asiz(1);
aspect(3) = Asiz(3) ./ Asiz(1);

xlim(1, :) = [0.5, Asiz(2) + 0.5]; %x
xlim(2, :) = [0.5, Asiz(2) + 0.5]; %x
xlim(3, :) = [0.5, Asiz(3) + 0.5]; %z

ylim(1, :) = [0.5, Asiz(3) + 0.5]; %z
ylim(2, :) = [0.5, Asiz(1) + 0.5]; %y
ylim(3, :) = [0.5, Asiz(1) + 0.5]; %y

Isiz(1, :) = [Asiz(3), Asiz(2)];
Isiz(2, :) = [Asiz(1), Asiz(2)];
Isiz(3, :) = [Asiz(1), Asiz(3)];

pointer{1} = 'arrow';
pointer{2} = 'cross';

markercolors{1, 1} = [0.5, 0.0, 0.0];
markercolors{2, 1} = [0.0, 0.5, 0.0];
markercolors{3, 1} = [0.0, 0.0, 0.5];

currentaxes = 0;
lastaxes = 1;

clim = [-1, 1];
MarkerSize = 10;


indlabel = [opt.Labels{1}, opt.Labels{2}, opt.Labels{3}, ': [%d, %d, %d] = %.2f'];


% abort flag
abortplay = false;
playing = false;

% create a new figure
if isempty(Hf)

    % figure
    % ----------------------
    % get the monitor size and the active monitor
    Hr = groot;
    mon = Hr.MonitorPositions;
    Nmon = size(mon, 1);

    if Nmon > 1
        p0 = Hr.PointerLocation;
        Imon = find(p0(1) >= mon(:, 1) & p0(1) < mon(:, 1) + mon(:, 3) & p0(2) >= mon(:, 2) & p0(2) < mon(:, 2) + mon(:, 4));

        if ~isempty(Imon)
            mon = mon(Imon, :);
        else
            mon = mon(1, :);
        end

    end

    if isempty(opt.FigHeight)
        fheight = 0.7 * mon(4);
    else
        fheight = opt.FigHeight;
    end

    % fig positions
    % 1 3
    % 2
    % 4 = colorbar
    % 5 = panel

    ftop = 0.06 * fheight;
    fbot = 0.06 * fheight;
    fsep = 0.01 * fheight;
    fleft = 0.06 * fheight;
    fright = 0.06 * fheight;

    fvert = fheight - (ftop + fbot + fsep);
    a1height = (Isiz(1,1) ./ (Isiz(1,1) + Isiz(2,1))) * fvert;
    a2height = (Isiz(2,1) ./ (Isiz(1,1) + Isiz(2,1))) * fvert;
    a3height = a2height;
    a4height = 0.08 * fheight;
    a5height = a1height - (a4height + 2*fsep + 4*fsep);

    a1bot = fbot;
    a2bot = a1bot + a1height + fsep;
    a3bot = a2bot;
    a4bot = a3bot - (a4height + fsep);
    a5bot = a1bot;

    a1width = a1height * aspect(1);
    a2width = a1width;
    a3width = a3height * aspect(3);
    a4width = a3width;
    a5width = a3width;

    fwidth = fleft + fright + a2width + fsep + a3width;

    a1left = fleft;
    a2left = a1left;
    a3left = a2left + a2width + fsep;
    a4left = a3left;
    a5left = a3left;

    gcapos(1, :) = [a1left, a1bot, a1width, a1height];
    gcapos(2, :) = [a2left, a2bot, a2width, a2height];
    gcapos(3, :) = [a3left, a3bot, a3width, a3height];
    gcapos(4, :) = [a4left, a4bot, a4width, a4height];
    gcapos(5, :) = [a5left, a5bot, a5width, a5height];

    gcfpos(4) = fheight; % Fig height
    gcfpos(3) = fwidth; % Fig width
    gcfpos(1) = mon(1) + (mon(3) - gcfpos(3)) / 2;
    gcfpos(2) = mon(2) + (mon(4) - gcfpos(4)) / 2;

    Hf = figure;
    Hf.Position = gcfpos;
    Hf.Units = 'pixels';
    Hf.Name = sprintf('%s [%d,%d,%d]', inputname(1), Asiz);
    Hf.MenuBar = 'none';
    Hf.ToolBar = 'figure';
    Hf.Tag = 'ortho_view';

    % add some options
    Hpanel = uipanel('Title', 'Options', 'Units', 'pixel', 'Position', gcapos(5, :), 'Tag', 'v3v_panel');
    set(Hpanel, 'Units', 'Normalized');
    set(Hpanel, 'Visible', OptionsVisible);

    x0 = linspace(0.01, 0.99, 6);
    y0 = linspace(0.01, 0.99, 9);
    % dx = 0.98 * (x0(2) - x0(1));
    dy = 0.98 * (y0(2) - y0(1));

    % CLim selection button
    row = 8;
    uicontrol('String', 'CLim mode (L)', ...
        'Style', 'text', ...
        'HorizontalAlignment', 'left', ...
        'Units', 'normalized', ...
        'Position', [x0(1), y0(row), 0.98 * (x0(3) - x0(1)), dy], ...
        'Parent', Hpanel);
    Hb.clim = uicontrol('String', climmodes, ...
        'ToolTipString', 'Select the color limit mode', ...
        'Style', 'popupmenu', ...
        'Value', climmode, ...
        'Units', 'normalized', ...
        'Position', [x0(3), y0(row), (x0(6) - x0(3)), dy], ...
        'Tag', 'v3v_clim', ...
        'call', @setclim, ...
        'Parent', Hpanel);
    row = row - 1;
    uicontrol('String', 'CMap (M)', ...
        'Style', 'text', ...
        'HorizontalAlignment', 'left', ...
        'Units', 'normalized', ...
        'Position', [x0(1), y0(row), 0.98 * (x0(3) - x0(1)), dy], ...
        'Parent', Hpanel);
    Hb.cmap = uicontrol('String', cmaps, ...
        'ToolTipString', 'Select a color map', ...
        'Style', 'popupmenu', ...
        'Value', cmap, ...
        'Units', 'normalized', ...
        'Position', [x0(3), y0(row), (x0(6) - x0(3)), dy], ...
        'Tag', 'v3v_cmap', ...
        'call', @setcmap, ...
        'Parent', Hpanel);
    row = row - 1;
    uicontrol('String', 'Flatten mode (F)', ...
        'Style', 'text', ...
        'HorizontalAlignment', 'left', ...
        'Units', 'normalized', ...
        'Position', [x0(1), y0(row), 0.98 * (x0(3) - x0(1)), dy], ...
        'Parent', Hpanel);
    Hb.flatten = uicontrol('String', flattenfun, ...
        'ToolTipString', 'This will apply the selected function on each of the image dimensions, creating three 2D images', ...
        'Style', 'popupmenu', ...
        'Value', 1, ...
        'Units', 'normalized', ...
        'Position', [x0(3), y0(row), (x0(6) - x0(3)), dy], ...
        'Tag', 'v3v_cmap', ...
        'call', @setflatten, ...
        'Parent', Hpanel);
    uicontrol('String', {
        'Click in any of the axes to update the position of the others'
        'The intersection locations are marked on the edge of the axes'
        ''
        'Left/Right: move by one slice'
        'Up/Down: move by 10 slices'
        'space: autoplay'
        '(the above 3 apply to the last active axis)'
        'h: hide this panel'
        }, ...
        'Style', 'text', ...
        'HorizontalAlignment', 'left', ...
        'Units', 'normalized', ...
        'Position', [x0(1), y0(1), (x0(6) - x0(1)), (y0(6) - y0(1))], ...
        'Parent', Hpanel);

    getclim();

    if ~verLessThan('matlab', '9.5')
        addToolbarExplorationButtons(Hf)
    end

    Hf.WindowButtonMotionFcn = @freemove;
    Hf.WindowButtonDownFcn = @mousebutdown;
    Hf.WindowButtonUpFcn = @mousebutup;
    Hf.WindowScrollWheelFcn = @mousescroll;
    Hf.KeyPressFcn = @keypress;

    for kk = 1:3
        Ha(kk) = axes('Units', 'pixels', 'Position', gcapos(kk, :), 'NextPlot', 'Add', 'Box', 'On', 'XTick', [], 'Ytick', [], 'Tag', sprintf('v3v_ax_%d', kk));
        set(Ha(kk), 'XLim', xlim(kk, :), 'YLim', ylim(kk, :));

        if ~verLessThan('matlab', '9.5')
            set(Ha(kk).Toolbar, 'Visible', 'off');
        end

        Hi(kk) = imagesc(zeros(Isiz(kk, :)), 'Parent', Ha(kk), 'Tag', sprintf('v3v_im_%d', kk));
        Hp(kk, 1) = plot(nan(2, 1), nan(2, 1), 'd', 'color', markercolors{axind(kk, 2)}, 'markerfacecolor', 'w', 'MarkerSize', MarkerSize, 'LineWidth', 2, 'Tag', sprintf('v3v_p1_%d', kk));
        Hp(kk, 2) = plot(nan(2, 1), nan(2, 1), 'd', 'color', markercolors{axind(kk, 1)}, 'markerfacecolor', 'w', 'MarkerSize', MarkerSize, 'LineWidth', 2, 'Tag', sprintf('v3v_p2_%d', kk));
    end
    set(Ha, 'YDir', 'reverse', 'DataAspectRatio', [1 1 1]);
    set(Ha,'CLim',clim);

    % Colorbar
    Cbarlim = [min(Aqtl(1),clim(1)), max(Aqtl(9),clim(2))];
    Cbar = linspace(Cbarlim(1),Cbarlim(2),256);
    kk = 4;
    Hca = axes('Units', 'pixels', 'Position', gcapos(kk, :), 'NextPlot', 'Add', 'Box', 'On', 'Ytick', [], 'Tag', 'v3v_ax_Hca');
    Hci = imagesc(Cbarlim, [0,1], Cbar, 'Tag', 'v3v_ax_Hci');
    Hch = plot(Ahistval, Ahist, 'k-', 'LineWidth', 2, 'Tag', 'v3v_ax_Hch');
    Hcl = line(clim .* [1 ; 1], repmat([0; 1],1,2),'Color','k','LineStyle','--', 'LineWidth', 1, 'Tag', 'v3v_ax_Hcl');
    Hcq = plot(Aqtl([2,3,4,6,7,8]), repmat(0.5,1,6), 'k.', Aqtl([1,5,9]), repmat(0.5,1,3), 'kd', 'MarkerSize', MarkerSize, 'LineWidth', 2, 'Tag', 'v3v_ax_Hcq');
    set(Hca,'CLim',clim);
    set(Hca,'XLim',Cbarlim);
    set(Hca,'YLim',[0, 1]);


    for kk = 1:3
        xlabel(Ha(kk), opt.Labels{axind(kk, 1)});
        ylabel(Ha(kk), opt.Labels{axind(kk, 2)});
    end

    set(Ha(2), 'XAxisLocation', 'top');
    set(Ha(3), 'XAxisLocation', 'top', 'YAxisLocation', 'right');

    %     set(Ha,'DataAspectRatio',[1,1,1])
    set(Ha, 'Units', 'Normalized');
    set(Hca, 'Units', 'Normalized');

    Ht = uicontrol('String', '', ...
        'Style', 'text', ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'Position', [5, 2, 200, 25], ...
        'Parent', Hf);

    setcmap();
    redraw();

else

    Hf.Name = sprintf('%s [%d,%d,%d]', inputname(1), Asiz);

    for kk = 1:3
        Ha(kk) = findobj(Hf, 'Type', 'Axes', 'Tag', sprintf('v3v_ax_%d', kk));
        Hi(kk) = findobj(Hf, 'Type', 'Image', 'Tag', sprintf('v3v_im_%d', kk));
        Hp(kk, 1) = findobj(Hf, 'Type', 'Line', 'Tag', sprintf('v3v_p1_%d', kk));
        Hp(kk, 2) = findobj(Hf, 'Type', 'Line', 'Tag', sprintf('v3v_p2_%d', kk));
    end

    Hca = findobj(Hf, 'Type', 'Axes', 'Tag', 'v3v_ax_Hca');
    Hci = findobj(Hf, 'Type', 'Image', 'Tag', 'v3v_ax_Hci');
    Hch = findobj(Hf, 'Type', 'Line', 'Tag', 'v3v_ax_Hch');
    Hcl = findobj(Hf, 'Type', 'Line', 'Tag', 'v3v_ax_Hcl');
    Hcq = findobj(Hf, 'Type', 'Line', 'Tag', 'v3v_ax_Hcq');

    Ht = findobj(Hf, 'Type', 'Text', 'Tag', 'v3v_t');
    Hb.clim = findobj(Hf, 'Type', 'UIControl', 'Tag', 'v3v_clim');
    Hb.cmap = findobj(Hf, 'Type', 'UIControl', 'Tag', 'v3v_cmap');
    Hpanel = findobj(Hf, 'Type', 'Panel', 'Tag', 'v3v_panel');
    clim = Hci.XData;

    setcmap();
    redraw();

end

if nargout == 1
    varargout{1} = Hf;
end

% ----------------------------------------------
    function redraw(varargin)
        ind = min(max(round(ind), 1), Asiz);

        if Hb.flatten.Value > 1
            set(Hp, 'Visible', 'off');
            set(Ht, 'String', sprintf('Flatten mode: %s', flattenfun{Hb.flatten.Value}));

            set(Hi(1), 'Cdata', permute(Aflat{1}(1, :, :), [3, 2, 1]));
            set(Hi(2), 'Cdata', permute(Aflat{3}(:, :, 1), [1, 2, 3]));
            set(Hi(3), 'Cdata', permute(Aflat{2}(:, 1, :), [1, 3, 2]));
            
        else
            set(Hp, 'Visible' ,'on');

            set(Ht, 'String', sprintf(indlabel, ind, A(ind(1), ind(2), ind(3))));

            set(Hi(1), 'Cdata', permute(A(ind(1), :, :), [3, 2, 1]));
            set(Hi(2), 'Cdata', permute(A(:, :, ind(3)), [1, 2, 3]));
            set(Hi(3), 'Cdata', permute(A(:, ind(2), :), [1, 3, 2]));

            set(Hp(1, 1), 'Xdata', xlim(1, :), 'YData', ind(3) * [1, 1]);
            set(Hp(1, 2), 'Ydata', ylim(1, :), 'XData', ind(2) * [1, 1]);

            set(Hp(2, 1), 'Xdata', xlim(2, :), 'YData', ind(1) * [1, 1]);
            set(Hp(2, 2), 'Ydata', ylim(2, :), 'XData', ind(2) * [1, 1]);

            set(Hp(3, 1), 'Xdata', xlim(3, :), 'YData', ind(1) * [1, 1]);
            set(Hp(3, 2), 'Ydata', ylim(3, :), 'XData', ind(3) * [1, 1]);
        end
        
        drawnow
    end

% ----------------------------------------------
    function p = getpoint(varargin)

        if currentaxes ~= 0
            lastaxes = currentaxes;
        end

        currentaxes = 0;

        for k = 1:3
            p = get(Ha(k), 'CurrentPoint');
            p = p(1, 1:2);
            xl = get(Ha(k), 'Xlim');
            yl = get(Ha(k), 'Ylim');

            % put the interactive border inside the frame
            xl = xl + 0.01 * diff(xl) * [1, - 1];
            yl = yl + 0.01 * diff(yl) * [1, - 2];

            if (p(1) > xl(1)) && (p(1) < xl(2)) && (p(2) > yl(1)) && (p(2) < yl(2))
                currentaxes = k;
                break
            end

        end

    end

% ----------------------------------------------
    function freemove(varargin)
        getpoint();

        if currentaxes == 0
            set(Hf, 'Pointer', pointer{1});
            return
        else
            set(Hf, 'Pointer', pointer{2});
        end

    end

% ----------------------------------------------
    function dragmove(varargin)
        
        p = getpoint();

        if currentaxes == 0
            set(Hf, 'Pointer', pointer{1});
            return
        else
            set(Hf, 'Pointer', pointer{2});
        end

        ind(axind(currentaxes, :)) = round(p);
        redraw();
    end

% ----------------------------------------------
    function mousebutdown(varargin)
        Hf.WindowButtonMotionFcn = @dragmove;
    end

% ----------------------------------------------
    function mousebutup(varargin)
        Hf.WindowButtonMotionFcn = @freemove;

        p = getpoint();

        if currentaxes == 0
            set(Hf, 'Pointer', pointer{1});
            return
        else
            set(Hf, 'Pointer', pointer{2});
        end

        % dispable compression mode
        Hb.flatten.Value = 1;
        
        ind(axind(currentaxes, :)) = round(p);
        redraw();
    end

% ----------------------------------------------
    function mousescroll(~, varargin)
        event = varargin{1};
        a = event.VerticalScrollCount;

        % p = getpoint();

        if currentaxes == 0
            set(Hf, 'Pointer', pointer{1});
            return
        else
            set(Hf, 'Pointer', pointer{2});
        end

        ind(axother(currentaxes)) = round(ind(axother(currentaxes)) + a);
        redraw();

    end

% ----------------------------------------------
    function keypress(varargin)
        event = varargin{2};

        if strcmpi(event.Key, 'h')

            if strcmpi(Hpanel.Visible, 'on')
                Hpanel.Visible = 'off';
            else
                Hpanel.Visible = 'on';
            end
            return
        end

        if strcmpi(event.Key, 'f')
            val = Hb.flatten.Value;
            if isempty(event.Modifier)
                val = val + 1;
            else
                val = val - 1;
            end
            val = mod(val-1,numel(flattenfun))+1;
            Hb.flatten.Value = val;
            setflatten();
            return
        end

        if strcmpi(event.Key, 'm')
            val = Hb.cmap.Value;
            if isempty(event.Modifier)
                val = val + 1;
            else
                val = val - 1;
            end
            val = mod(val-1,numel(cmaps))+1;
            Hb.cmap.Value = val;
            setcmap();
            return
        end

        if strcmpi(event.Key, 'l')
            val = Hb.clim.Value;
            if isempty(event.Modifier)
                val = val + 1;
            else
                val = val - 1;
            end
            val = mod(val-1,numel(climmodes))+1;
            Hb.clim.Value = val;
            setclim();
            return
        end
        

        if strcmpi(event.Key, 'space')
            % space -> autoplay
            % if already playing, -> abort play
            % if not yet playing, -> start play
            if playing
                abortplay = ~abortplay;
            else
                autoplay();
            end

            return
        end

        a = 0;

        if strcmpi(event.Key, 'leftarrow')
            a = -1;
        elseif strcmpi(event.Key, 'rightarrow')
            a = 1;
        elseif strcmpi(event.Key, 'uparrow')
            a = 10;
        elseif strcmpi(event.Key, 'downarrow')
            a = -10;
        end

        ind(axother(lastaxes)) = round(ind(axother(lastaxes)) + a);
        redraw();
    end

% ----------------------------------------------
    function autoplay(varargin)
        playing = true;
        i = axother(lastaxes);
        j = ind(i);
        N = Asiz(i);
        I = mod((1:N) + j - 1, N) + 1;

        for j = I

            if abortplay
                abortplay = false;
                playing = false;
                return
            end

            ind(i) = j;
            redraw();
        end

        playing = false;
    end

% ----------------------------------------------
    function getclim(varargin)

        if ~isempty(opt.clim)
            clim = opt.clim;
            return
        end

        val = Hb.clim.Value;

        switch val
            case 1 % 'Auto 100%';
                clim = Aqtl([1, 9]);
            case 2 % 'Auto 99.9%';
                clim = Aqtl([2, 8]);
            case 3 % 'Auto 99.0%';
                clim = Aqtl([3, 7]);
            case 4 % 'Auto 90.0%';
                clim = Aqtl([4, 6]);
            case 5 % 'Centered 100%';
                clim = max(abs(Aqtl([1, 9]))) * [-1, 1];
            case 6 % 'Centered 99.9%';
                clim = max(abs(Aqtl([2, 8]))) * [-1, 1];
            case 7 % 'Centered 99.0%';
                clim = max(abs(Aqtl([3, 7]))) * [-1, 1];
            case 8 % 'Centered 90.0%';
                clim = max(abs(Aqtl([4, 6]))) * [-1, 1];
            case 9
                % manual, do nothing
                return
            otherwise
                error('Unknown method.')
        end

        if clim(1) == clim(2)
            clim = clim(1) + 0.01 * [-1, 1];
        end

    end

% ----------------------------------------------
    function setclim(varargin)

        if Hb.clim.Value == 9
            prompt = {'CLim(1):', 'CLim(2):'};
            name = 'Input the color limits';
            defaultanswer = {sprintf('%g', clim(1)), sprintf('%g', clim(2))};
            numlines = 1;
            answer = inputdlg(prompt, name, numlines, defaultanswer);

            if isempty(answer)
                return
            end

            clim = [eval(answer{1}), eval(answer{2})];
            clim = sort(clim);
        end

        opt.clim = [];

        getclim();
        set(Ha, 'CLim', clim);
        set(Hca, 'CLim', clim);

        set(Hcl(1),'XData',clim(1) * [1, 1]);
        set(Hcl(2),'XData',clim(2) * [1, 1]);

        Cbarlim = [min(Aqtl(1),clim(1)), max(Aqtl(9),clim(2))];
        Cbar = linspace(Cbarlim(1),Cbarlim(2),256);
        set(Hca,'XLim',Cbarlim);
        set(Hci,'XData',Cbarlim,'CData',Cbar);
    end

% ----------------------------------------------
    function setflatten(varargin)
        val = Hb.flatten.Value;

        if val == 1
            redraw();
            return
        elseif val == 2
            fun = @(A,i) max(A,[],i,'omitnan');
        elseif val == 3
            fun = @(A,i) min(A,[],i,'omitnan');
        elseif val == 4
            fun = @(A,i) max(A,[],i,'omitnan');
        elseif val == 5
            fun = @(A,i) median(A,i,'omitnan');
        elseif val == 6
            fun = @(A,i) mean(A,i,'omitnan');
        elseif val == 7
            fun = @(A,i) rms(A,i,'omitnan');
        elseif val == 8
            fun = @(A,i) std(A,i,'omitnan');
        end

        for k = 1:3
            Aflat{k} = fun( A , k);
        end
        redraw();
    end


% ----------------------------------------------
    function setcmap(varargin)
        val = Hb.cmap.Value;
        colormap(colorbrewer(cmaps{val}, Ncolors));

        set(Hch,'Color',cmap_linecol{val});
        set(Hcl,'Color',cmap_linecol{val});
        set(Hcq,'Color',cmap_linecol{val});
    end

end
