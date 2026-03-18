function [Data_reg, info] = registerStructures3D(Data, AtomName, varargin)
% REGISTERSTRUCTURES3D
% Rigid registration of molecular xyz structures using weighted Kabsch alignment.
%
% Supports:
%   Data size [N x 3 x Nsample]
%   Data size [N x 3 x Ntime x Ntraj]
%
% Default behavior:
%   - all atoms are used for alignment
%   - weights are atomic masses from AtomicMass
%   - no symmetry permutations are assumed
%
% The transform is determined from:
%   - the only frame, if Data is [N x 3 x Nsample]
%   - the first time point, if Data is [N x 3 x Ntime x Ntraj]
%
% Then the same transform is applied to all time points of that sample/trajectory.

    %% ---------- atom names ----------
    AtomName = string(AtomName(:));
    N = numel(AtomName);

    %% ---------- data shape ----------
    if ndims(Data) == 3
        % [N x 3 x Nsample]
        if size(Data,1) ~= N || size(Data,2) ~= 3
            error('For 3D input, Data must be [N x 3 x Nsample].');
        end
        is3Dinput = true;
        Nsample = size(Data,3);
        Ntime = 1;
        Ntraj = Nsample;
    elseif ndims(Data) == 4
        % [N x 3 x Ntime x Ntraj]
        if size(Data,1) ~= N || size(Data,2) ~= 3
            error('For 4D input, Data must be [N x 3 x Ntime x Ntraj].');
        end
        is3Dinput = false;
        Ntime = size(Data,3);
        Ntraj = size(Data,4);
    else
        error('Data must be [N x 3 x Nsample] or [N x 3 x Ntime x Ntraj].');
    end

    %% ---------- parse inputs ----------
    p = inputParser;
    p.addParameter('alignIdx', 1:N, @(x) isnumeric(x) && isvector(x));
    p.addParameter('plotIdx', 1:N, @(x) isnumeric(x) && isvector(x));
    p.addParameter('referenceIndex', 1, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('referenceCoords', [], @(x) isnumeric(x) && isequal(size(x), [N,3]));
    p.addParameter('weights', 'mass', @(x) isnumeric(x) || ischar(x) || isstring(x));
    p.addParameter('symmetryPermutations', {}, @(x) iscell(x));
    p.addParameter('showPlot', true, @(x) islogical(x) || isnumeric(x));
    p.addParameter('alphaVal', 0.15, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
    p.addParameter('markerScale', 110, @(x) isnumeric(x) && isscalar(x));
    p.parse(varargin{:});
    S = p.Results;

    S.alignIdx = S.alignIdx(:).';
    S.plotIdx  = S.plotIdx(:).';

    %% ---------- get masses from AtomicMass ----------
    mass = zeros(N,1);
    for i = 1:N
        mass(i) = AtomicMass(AtomName(i));
    end

    %% ---------- weights ----------
    if isnumeric(S.weights)
        wAll = S.weights(:);
        if numel(wAll) ~= N
            error('Numeric weights must have length N.');
        end
    else
        switch lower(string(S.weights))
            case "mass"
                wAll = mass;
            case "unit"
                wAll = ones(N,1);
            otherwise
                error('weights must be ''mass'', ''unit'', or a numeric vector of length N.');
        end
    end

    wAlign = wAll(S.alignIdx);

    %% ---------- reference structure ----------
    if isempty(S.referenceCoords)
        if is3Dinput
            Xref_full = Data(:,:,S.referenceIndex);
        else
            Xref_full = Data(:,:,1,S.referenceIndex);
        end
    else
        Xref_full = S.referenceCoords;
    end

    Xref_align = Xref_full(S.alignIdx,:);
    cref = weightedCenter(Xref_align, wAlign);
    Xref_align_c = Xref_align - cref;

    %% ---------- symmetry ----------
    nA = numel(S.alignIdx);
    if isempty(S.symmetryPermutations)
        S.symmetryPermutations = {1:nA};
    end

    for k = 1:numel(S.symmetryPermutations)
        pk = S.symmetryPermutations{k}(:).';
        if numel(pk) ~= nA || ~isequal(sort(pk), 1:nA)
            error('Each symmetry permutation must be a permutation of 1:numel(alignIdx).');
        end
    end

    %% ---------- allocate ----------
    Data_reg = nan(size(Data), 'like', Data);

    info = struct;
    info.mass = mass;
    info.weights = wAll;
    info.alignIdx = S.alignIdx;
    info.referenceIndex = S.referenceIndex;
    info.referenceCenter = cref;
    info.rotation = nan(3,3,Ntraj);
    info.sourceCenter = nan(Ntraj,3);
    info.translation = nan(Ntraj,3);
    info.rmsd = nan(Ntraj,1);
    info.bestPermutation = cell(Ntraj,1);

    %% ---------- register each structure / trajectory ----------
    for itraj = 1:Ntraj

        if is3Dinput
            X0_full = Data(:,:,itraj);
        else
            X0_full = Data(:,:,1,itraj);
        end

        X0_align = X0_full(S.alignIdx,:);
        cX = weightedCenter(X0_align, wAlign);
        X0_align_c = X0_align - cX;

        bestR = eye(3);
        bestPerm = 1:nA;
        bestRMSD = inf;

        for ip = 1:numel(S.symmetryPermutations)
            perm = S.symmetryPermutations{ip};
            Xtest = X0_align_c(perm,:);
            [R, rmsdVal] = kabschWeighted(Xtest, Xref_align_c, wAlign(perm));

            if rmsdVal < bestRMSD
                bestRMSD = rmsdVal;
                bestR = R;
                bestPerm = perm;
            end
        end

        if is3Dinput
            Data_reg(:,:,itraj) = (Data(:,:,itraj) - cX) * bestR + cref;
        else
            for it = 1:Ntime
                Data_reg(:,:,it,itraj) = (Data(:,:,it,itraj) - cX) * bestR + cref;
            end
        end

        info.rotation(:,:,itraj) = bestR;
        info.sourceCenter(itraj,:) = cX;
        info.translation(itraj,:) = cref - cX*bestR;
        info.rmsd(itraj) = bestRMSD;
        info.bestPermutation{itraj} = bestPerm;
    end

    %% ---------- plot ----------
    if S.showPlot
        if is3Dinput
            Xplot = Data_reg;                 % [N x 3 x Nsample]
            titleStr = 'Registered structures';
        else
            Xplot = squeeze(Data_reg(:,:,1,:)); % [N x 3 x Ntraj]
            titleStr = 'Registered structures at t = 1';
        end

        plotRegisteredEnsemble(Xplot, AtomName, ...
            'plotIdx', S.plotIdx, ...
            'alphaVal', S.alphaVal, ...
            'markerScale', S.markerScale, ...
            'titleStr', titleStr);
    end
end

%% ========================================================================
function c = weightedCenter(X, w)
    w = w(:);
    c = sum(X .* w, 1) ./ sum(w);
end

%% ========================================================================
function [R, rmsdVal] = kabschWeighted(X, Y, w)
% X,Y are centered [N x 3] arrays
% finds R minimizing ||X*R - Y|| with weights w

    w = w(:);
    H = X' * (Y .* w);

    [U,~,V] = svd(H);
    R = U * diag([1 1 sign(det(U*V'))]) * V';

    D = X*R - Y;
    rmsdVal = sqrt(sum(w .* sum(D.^2,2)) / sum(w));
end

%% ========================================================================
function plotRegisteredEnsemble(X, AtomName, varargin)
% X is [N x 3 x Nsample]

    p = inputParser;
    p.addParameter('plotIdx', 1:size(X,1), @(x) isnumeric(x) && isvector(x));
    p.addParameter('alphaVal', 0.15, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('markerScale', 110, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('titleStr', 'Registered structures', @(x) ischar(x) || isstring(x));
    p.parse(varargin{:});
    S = p.Results;

    AtomName = string(AtomName(:));
    S.plotIdx = S.plotIdx(:).';

    [~,~,Nsample] = size(X);

    typeList = unique(AtomName(S.plotIdx), 'stable');
    colorMap = defaultAtomColors(typeList);

    figure('Color','w');
    ax = axes('NextPlot','add');
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

        for isamp = 1:Nsample
            Xi = X(idx,:,isamp);
            scatter3(ax, Xi(:,1), Xi(:,2), Xi(:,3), ...
                S.markerScale, repmat(clr,numel(idx),1), ...
                'filled', ...
                'MarkerFaceAlpha', S.alphaVal, ...
                'MarkerEdgeAlpha', max(0.05, 0.5*S.alphaVal));
        end

        hLeg(itype) = scatter3(ax, nan, nan, nan, S.markerScale, clr, 'filled');
    end

    legend(hLeg, cellstr(typeList), 'Location','bestoutside');
    camlight headlight;
    lighting gouraud;
end

%% ========================================================================
function cmap = defaultAtomColors(typeList)

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