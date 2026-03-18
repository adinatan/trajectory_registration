clear all
close all
%%
load test_traj_dataset.mat
 
% Solute subset used to define the rigid transform
idSol = 1:38;
AtomName = AtomNames(idSol);

% Register using solute only, from the first time point of each trajectory
[Data_reg, info] = registerStructures3D(xyzt(idSol,:,1,:), AtomName, ...
    'weights', 'mass', 'showPlot', true);

%  %  or use a subset of atoms 
%  alignIdx = [ 1 2 13  16 29  ];   
% [Data_reg, info] = registerStructures3D(Data, AtomName, ...
%     'alignIdx', alignIdx, ...
%     'weights', 'mass', ...
%     'showPlot', true);
 

% Use the registered solute reference to define one common final frame
Xref = Data_reg(:,:,1,1);          % registered reference solute structure
m = info.mass(:);                  % masses for the solute atoms

% Mass-weighted center of mass of the registered reference
c0 = sum(Xref .* m, 1) ./ sum(m);

% Center reference at origin
Xc = Xref - c0;

% Mass-weighted covariance -> principal axes
C = (Xc' * (Xc .* m)) / sum(m);
[V,D] = eig(C);
[~,ord] = sort(diag(D), 'descend');
V = V(:,ord);

% Largest principal axis -> z
zOld = V(:,1); 
zOld = zOld / norm(zOld);

% Second principal axis -> x, made orthogonal to z
xOld = V(:,2) - zOld * (zOld.' * V(:,2));
xOld = xOld / norm(xOld);

% Third axis -> y, then re-orthogonalize x
yOld = cross(zOld, xOld);
yOld = yOld / norm(yOld);
xOld = cross(yOld, zOld);
xOld = xOld / norm(xOld);

% Fix the sign of z so the reference points consistently toward +z
if mean(Xc * zOld) < 0
    zOld = -zOld;
    yOld = -yOld;   % keep a right-handed frame
end

% Row-vector convention: Xnew = (X - c0) * R0
R0 = [xOld, yOld, zOld];

% Apply the per-trajectory registration + the common final frame
[Natom,~,Ntime,Ntraj] = size(xyzt);
xyzt_final = nan(size(xyzt), 'like', xyzt);

for n4 = 1:Ntraj
    R = info.rotation(:,:,n4);     % trajectory-specific rotation
    t = info.translation(n4,:);    % trajectory-specific translation

    % Stack all time points of one trajectory into one 2D array
    X = permute(xyzt(:,:,:,n4), [1 3 2]);   % [Natom x Ntime x 3]
    X = reshape(X, [], 3);                  % [Natom*Ntime x 3]

    % First register to the common solute reference, then center/orient globally
    X = (X * R + t - c0) * R0;

    % Reshape back to [Natom x 3 x Ntime]
    X = reshape(X, Natom, Ntime, 3);
    xyzt_final(:,:,:,n4) = permute(X, [1 3 2]);
end

return % comment this out if you want to make that movie
%% make a small movie of how it all evolves
 
fig = figure('Color','w');
ax = axes('Parent',fig);
sb = 3.6; % box size in Ang
v = VideoWriter('ensemble_movie.mp4','MPEG-4');
v.FrameRate=size(xyzt_final,3)/8;
open(v)
%idSol
for nt = 1:size(xyzt_final,3)
    cla(ax)

    plotRegisteredEnsembleStandalone(xyzt_final(idSol,:,nt,:), AtomNames(idSol), ...
        'ax', ax, 'titleStr', sprintf('t = %d fs', (nt-1)*100));

    xlim([-sb sb])
    ylim([-sb sb])
    zlim([-sb sb])
    axis square
    view(10,15)
    drawnow
    writeVideo(v, getframe(fig));
end

close(v)