function plotRegisteredEnsembleStandalone(X, AtomName, varargin)
% Plot an ensemble of 3D structures with semi-transparent atoms by atom type.
%
% Inputs
% ------
% X
%   Either:
%     [N x 3 x Nsample]
%     [N x 3 x Ntime x Ntraj]
%
% AtomName
%   [N x 1] atom names, e.g. {'C','H','O',...}
%
% Name-value options
% ------------------
% 'plotIdx'     : atom indices to plot (default all)
% 'timeIndex'   : which time to use if X is 4D (default 1)
% 'alphaVal'    : marker transparency (default 0.15)
% 'markerScale' : marker size (default 110)
% 'titleStr'    : figure title
% 'ax'          : existing axes handle to plot into (default = create new figure)
%
% Examples
% --------
% plotRegisteredEnsembleStandalone(xyzt_final(:,:,1:10), AtomNames)
% plotRegisteredEnsembleStandalone(xyzt_final, AtomNames, 'timeIndex', 1)
% plotRegisteredEnsembleStandalone(xyzt_final, AtomNames, 'timeIndex', 1, 'plotIdx', 1:38)
%
% Reuse same figure in a loop:
%   fig = figure('Color','w');
%   ax = axes('Parent',fig);
%   for nt = 1:size(xyzt_final,3)
%       cla(ax)
%       plotRegisteredEnsembleStandalone(xyzt_final(1:38,:,nt,1:10), AtomNames(1:38), ...
%           'ax', ax, 'titleStr', sprintf('t = %d', nt));
%       drawnow
%   end

    AtomName = string(AtomName(:));
    N = numel(AtomName);

    p = inputParser;
    p.addParameter('plotIdx', 1:N, @(x) isnumeric(x) && isvector(x));
    p.addParameter('timeIndex', 1, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('alphaVal', 0.15, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
    p.addParameter('markerScale', 110, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('titleStr', 'Registered structures', @(x) ischar(x) || isstring(x));
    p.addParameter('ax', [], @(x) isempty(x) || isgraphics(x,'axes'));
    p.parse(varargin{:});
    S = p.Results;

    S.plotIdx = S.plotIdx(:).';

    % Convert input to [N x 3 x Nsample]
    if ndims(X) == 3
        if size(X,1) ~= N || size(X,2) ~= 3
            error('For 3D input, X must be [N x 3 x Nsample].');
        end
        Xplot = X;
    elseif ndims(X) == 4
        if size(X,1) ~= N || size(X,2) ~= 3
            error('For 4D input, X must be [N x 3 x Ntime x Ntraj].');
        end
        Xplot = squeeze(X(:,:,S.timeIndex,:));   % [N x 3 x Ntraj]
    else
        error('X must be [N x 3 x Nsample] or [N x 3 x Ntime x Ntraj].');
    end

    % Use existing axes if provided, otherwise create figure/axes
    if isempty(S.ax)
        figure('Color','w');
        ax = axes;
    else
        ax = S.ax;
    end

    typeList = unique(AtomName(S.plotIdx), 'stable');
    colorMap = defaultAtomColorsStandalone(typeList);

    hold(ax,'on');
    axis(ax,'equal');
    grid(ax,'on');
    view(ax,3);
    xlabel(ax,'x');
    ylabel(ax,'y');
    zlabel(ax,'z');
    title(ax,S.titleStr);

    hLeg = gobjects(numel(typeList),1);

    for itype = 1:numel(typeList)
        thisType = typeList(itype);
        idx = S.plotIdx(AtomName(S.plotIdx) == thisType);
        clr = colorMap(char(thisType));

        for isamp = 1:size(Xplot,3)
            Xi = Xplot(idx,:,isamp);
            scatter3(ax, Xi(:,1), Xi(:,2), Xi(:,3), ...
                S.markerScale, repmat(clr, numel(idx), 1), ...
                'filled', ...
                'MarkerFaceAlpha', S.alphaVal, ...
                'MarkerEdgeAlpha', max(0.05, 0.5*S.alphaVal));
        end

        hLeg(itype) = scatter3(ax, nan, nan, nan, S.markerScale, clr, 'filled');
    end

    legend(ax, hLeg, cellstr(typeList), 'Location', 'bestoutside');
    camlight(ax,'headlight');
    lighting(ax,'gouraud');
end

%% ------------------------------------------------------------------------
function cmap = defaultAtomColorsStandalone(typeList)

    cmap = containers.Map;
    base = containers.Map;

    base('H')  = [0.90 0.90 0.90];
    base('C')  = [0.20 0.20 0.20];
    base('N')  = [0.15 0.25 0.90];
    base('O')  = [0.90 0.10 0.10];
    base('P')  = [1.00 0.50 0.00];
    base('S')  = [0.95 0.80 0.10];
    base('Cl') = [0.10 0.70 0.10];
    base('Br') = [0.60 0.20 0.10];
    base('I')  = [0.50 0.00 0.70];
    base('Pt') = [0.60 0.60 0.70];

    for k = 1:numel(typeList)
        key = char(typeList(k));
        if isKey(base, key)
            cmap(key) = base(key);
        else
            rng(sum(double(key)), 'twister');
            cmap(key) = 0.2 + 0.6*rand(1,3);
        end
    end
end
