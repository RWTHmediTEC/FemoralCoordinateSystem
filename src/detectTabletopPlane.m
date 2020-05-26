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
visu = p.Results.visualization;

%% Landmarks
MPC = femur.vertices(LMIdx.MedialPosteriorCondyle,:);
LPC = femur.vertices(LMIdx.LateralPosteriorCondyle,:);
ICN = femur.vertices(LMIdx.IntercondylarNotch,:);

%% Construction of an inital system
MechanicalAxis = createLine3d(ICN, HJC);

%% Inital transformation (iTFM)
% Connection of the most posterior points of the condyles
PosteriorCondyleAxis = createLine3d(MPC, LPC);

Z = normalizeVector3d(MechanicalAxis(4:6));
Y = normalizeVector3d(crossProduct3d(MechanicalAxis(4:6), PosteriorCondyleAxis(4:6)));
X = normalizeVector3d(crossProduct3d(Y, Z));

iTFM = [[X;Y;Z],[0 0 0]'; [0 0 0 1]]*createTranslation3d(-HJC);
% If it is a left femur, mirror to right side
if strcmp(side, 'L')
    mirrorTFM = eye(4); mirrorTFM(1,1) = -1;
    iTFM = mirrorTFM*createRotationOz(pi)*iTFM;
end

% Transform the mesh by the inital TFM
iFemur = transformPoint3d(femur, iTFM);
if strcmp(side, 'L')
    iFemur.faces = fliplr(iFemur.faces);
end

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

% Sanity check
tabletopPlane = createPlane(LPC, MPC, PTC);
assert(sum(isBelowPlane(MCMesh.vertices, tabletopPlane)) <= 1)
assert(any(isBelowPlane(LCMesh.vertices, tabletopPlane)) <= 1)
assert(any(isBelowPlane(proximalPart.vertices, tabletopPlane)) <= 1)

% Get the vertex indices of the tabletop points
MPC_Idx = find(ismembertol(femurCS.vertices, MPC, 'ByRows',true'));
LPC_Idx = find(ismembertol(femurCS.vertices, LPC, 'ByRows',true'));
PTC_Idx = find(ismembertol(femurCS.vertices, PTC, 'ByRows',true'));

LMIdx.MedialPosteriorCondyle = MPC_Idx;
LMIdx.LateralPosteriorCondyle = LPC_Idx;
LMIdx.PosteriorTrochantericCrest = PTC_Idx;

%% visualization
if visu
    % Patch properties
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [223, 206, 161]/255;
    patchProps.FaceAlpha = 0.75;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    [~,axH] = visualizeMeshes(femurCS, patchProps);
    drawAxis3d(axH,35,1.5)
    
    % Landmarks
    drawPoint3d([MPC; LPC; PTC],'MarkerFaceColor','k','MarkerEdgeColor','k')
    % Axes
    drawLine3d(transformLine3d(NeckAxis, TFM))
    % Tabletop patch
    patchProps.LineStyle='none';
    patchProps.LineWidth = 2;
    patchProps.Marker='o';
    patchProps.MarkerFaceColor='k';
    patchProps.MarkerEdgeColor='k';
    patchProps.FaceColor='k';
    patchProps.FaceAlpha=0.75;
    patchProps.EdgeColor='none';
    
    tablePatch.vertices=[...
        femurCS.vertices(MPC_Idx,:);...
        femurCS.vertices(LPC_Idx,:);...
        femurCS.vertices(PTC_Idx,:)];
    tablePatch.faces=1:3;
    
    patch(tablePatch, patchProps)
    
    textPosX=1/3*(MPC(1)+LPC(1)+PTC(1));
    textPosY=MPC(2);
    textPosZ=1/3*(MPC(3)+LPC(3)+PTC(3));
    
    text(textPosX, textPosY, textPosZ, 'Tabletop plane','Rotation',90)
    
    text(tablePatch.vertices(:,1),tablePatch.vertices(:,2),tablePatch.vertices(:,3),...
        {'MPC';'LPC';'PTC'})
    
    anatomicalViewButtons('RAS')
end

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