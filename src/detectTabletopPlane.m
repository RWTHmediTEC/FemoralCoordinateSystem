function LMIdx = detectTabletopPlane(femur, side, HJC, LMIdx, varargin)

% Definition of the tabletop plane:
% - Resection of the femoral neck including the head
% - Positioning of the posterior side of the femur on a table
% - The three contact points define the tabletop plane

% Inputs
p = inputParser;
addRequired(p,'femur',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addRequired(p,'side',@(x) any(validatestring(x,{'R','L'})));
addOptional(p,'visualization',true,@islogical);
parse(p,femur,side,varargin{:});

femur = p.Results.femur;
visu = logical(p.Results.visualization);

%% Inital transformation (iTFM)
iTFM = Tabletop(femur, side, HJC, LMIdx, 'visu', visu);

% Transform the mesh by the inital TFM
iFemur = transformPoint3d(femur, iTFM);

NeckAxis = createLine3d(iFemur.vertices(LMIdx.NeckAxis(1),:),iFemur.vertices(LMIdx.NeckAxis(2),:));
NeckOrthogonal = createLine3d(iFemur.vertices(LMIdx.NeckOrthogonal(1),:),iFemur.vertices(LMIdx.NeckOrthogonal(2),:));
[~, NeckAxis(1:3), ~] = distanceLines3d(NeckAxis, NeckOrthogonal);
ICN = iFemur.vertices(LMIdx.IntercondylarNotch,:);

if visu
    % Patch properties
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [223, 206, 161]/255;
    patchProps.FaceAlpha = 0.5;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    [~,axH] = visualizeMeshes(iFemur, patchProps);
    drawAxis3d(axH, 35, 1.5)
    
    % Axes
    drawVector3d(NeckAxis(1:3),NeckAxis(4:6))
    
    anatomicalViewButtons('RAS')
end

%% Resect condyles
% Get the length of the femur
iLength = abs(max(iFemur.vertices(:,3)))+abs(min(iFemur.vertices(:,3)));
% Cut off the distal part
DISTAL_FACTOR = 1/6;
distalPlane=[0 0 DISTAL_FACTOR*iLength+min(iFemur.vertices(:,3)), 1 0 0, 0 1 0];
distalPart = cutMeshByPlane(iFemur, distalPlane, 'part','below');
% Cut the distal part into the medial and lateral condyle
sagittalPlane=createPlane(ICN, [1 0 0]);
[LCMesh, ~, MCMesh]  = cutMeshByPlane(distalPart, sagittalPlane);
if visu
    patchProps.FaceAlpha = 1;
    patchProps.FaceColor = 'r';
    patch(LCMesh, patchProps);
    patch(MCMesh, patchProps);
end

%% Resect proximal femur
PROXIMAL_FACTOR = 2/6;
proximalPlane=[0 0 -PROXIMAL_FACTOR*iLength, 1 0 0, 0 1 0];
proximalPart = cutMeshByPlane(iFemur, proximalPlane, 'part','above');
% Resect the neck and the head
if NeckAxis(6)>0; NeckAxis(4:6)=-NeckAxis(4:6); end
neckPlane = createPlane(NeckAxis(1:3), NeckAxis(4:6));
proximalPart = cutMeshByPlane(proximalPart, neckPlane, 'part','above');
if visu
    patch(proximalPart, patchProps);
end

%% Refinement of tabletop plane
[tempROT, MPC, LPC, PTC] = refineTableTopPlane(MCMesh, LCMesh, proximalPart);
refROT = tempROT;
while ~isequal(eye(4), tempROT)
    MCMesh = transformPoint3d(MCMesh, tempROT);
    LCMesh = transformPoint3d(LCMesh, tempROT);
    proximalPart = transformPoint3d(proximalPart, tempROT);
    [tempROT, MPC, LPC, PTC] = refineTableTopPlane(MCMesh, LCMesh, proximalPart);
    refROT = tempROT*refROT;
end

% Combine all transformations
TFM = refROT*iTFM;

% The femur in the AFCS
femurCS = transformPoint3d(femur, TFM);

% Sanity checks
tabletopPlane = createPlane(LPC, MPC, PTC);
assert(sum(isBelowPlane(MCMesh.vertices, tabletopPlane)) <= 1)
assert(any(isBelowPlane(LCMesh.vertices, tabletopPlane)) <= 1)
assert(any(isBelowPlane(proximalPart.vertices, tabletopPlane)) <= 1)

% Get the vertex indices of the tabletop points
MPC_Idx = find(ismembertol(femurCS.vertices, MPC, 'ByRows',true'));
LPC_Idx = find(ismembertol(femurCS.vertices, LPC, 'ByRows',true'));
PTC_Idx = find(ismembertol(femurCS.vertices, PTC, 'ByRows',true'));

switch side
    case 'R'
        LMIdx.MedialPosteriorCondyle = MPC_Idx;
        LMIdx.LateralPosteriorCondyle = LPC_Idx;
    case 'L'
        LMIdx.MedialPosteriorCondyle = LPC_Idx;
        LMIdx.LateralPosteriorCondyle = MPC_Idx;
end
LMIdx.PosteriorTrochantericCrest = PTC_Idx;

%% visualization
Tabletop(femur, side, HJC, LMIdx, 'visu', visu);


end

function [ROT, MPC, LPC, PTC] = refineTableTopPlane(MCMesh, LCMesh, proximalPart)

% Get the index of the most posterior point of the condyle
[~, MPC_Idx] = min(MCMesh.vertices(:,2));
[~, LPC_Idx] = min(LCMesh.vertices(:,2));
% Get the most posterior point of the condyle
MPC = MCMesh.vertices(MPC_Idx,:);
LPC = LCMesh.vertices(LPC_Idx,:);
PosteriorCondyleAxis=createLine3d(MPC, LPC);

% Get the index of the most posterior point of the trochanteric crest
[~, PTC_Idx] = min(proximalPart.vertices(:,2));
PTC = proximalPart.vertices(PTC_Idx,:);

% Create the tabletop plane
tabletopPlane = createPlane(MPC, PTC, LPC);
tabletopNormal = planeNormal(tabletopPlane);
% Normal has to point in anterior direction
if tabletopNormal(2) < 0; tabletopNormal=-tabletopNormal; end
X = normalizeVector3d(PosteriorCondyleAxis(4:6));
Y = normalizeVector3d(tabletopNormal);
Z = normalizeVector3d(crossProduct3d(X, Y));

ROT = [[X; Y; Z; 0 0 0], [0 0 0 1]'];

end