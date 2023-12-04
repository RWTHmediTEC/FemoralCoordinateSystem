function LM_Idx = detectTabletopPlane(femur, side, FHC, neckAxis, LM_Idx, varargin)
%DETECTTABLETOPPLANE identifies the table top plane of the femur
% 
%   Definition of the tabletop plane:
%	  - Resection of the femoral neck including the head (in case 
%     retroversion is present)
%     - Position the posterior side of the femur on a table
%     - The three contact points define the tabletop plane
%
%   INPUT
%     femur: femur mesh
%     side: 'L' or 'R': left or right femur
%     FHC: femoral head center
%     neckAxis: neckAxis(1:3) is the center of the neck isthmus and 
%     neckAxis(4:6) is the normal of the neck istmus plane (see FNA.m)
%     LM_Idx: Struct containing the vertex indices of medial and lateral
%     epicondyle mapped from the template.
%
%   OUTPUT
%     LM_Idx = A struct containing the vertex indices of the femoral mesh 
%     defining the three contact points with the fields:
%     - MedialPosteriorCondyle
%     - LateralPosteriorCondyle
%     - PosteriorTrochantericCrest
%
% AUTHOR: Maximilian C. M. Fischer
% COPYRIGHT (C) 2020-2023 Maximilian C. M. Fischer
% LICENSE: EUPL v1.2
%

% Inputs
p = inputParser;
addRequired(p,'femur',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addRequired(p,'side',@(x) any(validatestring(x,{'R','L'})));
addOptional(p,'visualization',true,@islogical);
parse(p,femur,side,varargin{:});

femur = p.Results.femur;
visu = logical(p.Results.visualization);

%% Inital transformation (iTFM)
iTFM = Tabletop(femur, side, FHC, LM_Idx, 'visu', false);

% Transform the mesh by the inital TFM
iFemur = transformPoint3d(femur, iTFM);
neckAxis = transformLine3d(neckAxis, iTFM);
MPC = iFemur.vertices(LM_Idx.MedialPosteriorCondyle,:);
LPC = iFemur.vertices(LM_Idx.LateralPosteriorCondyle,:);

%% Resect condyles
% Get the length of the femur
iLength = abs(max(iFemur.vertices(:,3)))+abs(min(iFemur.vertices(:,3)));
% Cut off the distal part
DISTAL_FACTOR = 1/6;
distalPlane = [0 0 DISTAL_FACTOR*iLength+min(iFemur.vertices(:,3)), 1 0 0, 0 1 0];
distalPart = cutMeshByPlane(iFemur, distalPlane, 'part','below');
% Cut the distal part into the medial and lateral condyle
sagittalPlane = createPlane(midPoint3d(MPC,LPC), [1 0 0]);
[LCMesh, ~, MCMesh]  = cutMeshByPlane(distalPart, sagittalPlane);

%% Resect proximal femur
PROXIMAL_FACTOR = 2/6;
proximalPlane = [0 0 -PROXIMAL_FACTOR*iLength, 1 0 0, 0 1 0];
[proximalPart, ~, shaft] = cutMeshByPlane(iFemur, proximalPlane);
shaft = cutMeshByPlane(shaft, distalPlane, 'part','above');
% Resect the neck and the head
if neckAxis(6)>0; neckAxis(4:6)=-neckAxis(4:6); end
neckPlane = createPlane(neckAxis(1:3), neckAxis(4:6));
[proximalPart,~,head] = cutMeshByPlane(proximalPart, neckPlane);

if visu
    % Patch properties
    patchProps.EdgeColor = 'none';
    patchProps.FaceAlpha = 0.5;
    patchProps.FaceLighting = 'gouraud';
    [~,axH,figH] = visualizeMeshes([head,shaft], patchProps); 
    % drawAxis3d(axH, 35, 1.5)
    figH.Name = 'Detection of the table top plane';
    figH.NumberTitle = 'Off';
    
    patchProps.FaceAlpha = 1;
    patchProps.FaceColor = 'r';
    patch(axH, LCMesh, patchProps);
    patch(axH, MCMesh, patchProps);
    % drawPoint3d(midPoint3d(MPC,LPC),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k')
    patch(axH, proximalPart, patchProps);
    % Neck axis
    drawArrow3d(axH, neckAxis(1:3),-normalizeVector3d(neckAxis(4:6))*40,'g') % reversed for better visu
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
assert(sum(isBelowPlane(LCMesh.vertices, tabletopPlane)) <= 1)
assert(sum(isBelowPlane(proximalPart.vertices, tabletopPlane)) <= 1)

% Get the vertex indices of the tabletop points
switch side
    case 'R'
        MPC_Idx = find(ismembertol(femurCS.vertices, MPC, 'ByRows',1));
        LPC_Idx = find(ismembertol(femurCS.vertices, LPC, 'ByRows',1));
    case 'L'
        LPC_Idx = find(ismembertol(femurCS.vertices, MPC, 'ByRows',1));
        MPC_Idx = find(ismembertol(femurCS.vertices, LPC, 'ByRows',1));
end
PTC_Idx = find(ismembertol(femurCS.vertices, PTC, 'ByRows',true'));

LM_Idx.MedialPosteriorCondyle = MPC_Idx;
LM_Idx.LateralPosteriorCondyle = LPC_Idx;
LM_Idx.PosteriorTrochantericCrest = PTC_Idx;

%% Visualization
if visu
    % Landmarks
    MPC = iFemur.vertices(MPC_Idx,:);
    LPC = iFemur.vertices(LPC_Idx,:);
    PTC = iFemur.vertices(PTC_Idx,:);
    text(axH, MPC(:,1),MPC(:,2)-5,MPC(:,3)-5,{'MPC'},'FontSize',14,'HorizontalAlignment','Right')
    text(axH, LPC(:,1),LPC(:,2)-5,LPC(:,3)-5,{'LPC'},'FontSize',14,'HorizontalAlignment','Right')
    text(axH, PTC(:,1),PTC(:,2)-5,PTC(:,3)+5,{'PTC'},'FontSize',14,'HorizontalAlignment','Left')
    
    % Tabletop patch
    patchProps.LineStyle='-';
    patchProps.LineWidth = 1;
    patchProps.Marker='o';
    patchProps.MarkerFaceColor='k';
    patchProps.MarkerEdgeColor='k';
    patchProps.FaceColor='k';
    patchProps.FaceAlpha=0.75;
    patchProps.EdgeColor='k';
    
    tablePatch.vertices = [MPC;LPC;PTC];
    tablePatch.faces = 1:3;
    
    patch(axH, tablePatch, patchProps)
    
    textPosX = 1/3*(MPC(1)+LPC(1)+PTC(1));
    textPosY = MPC(2);
    textPosZ = 1/3*(MPC(3)+LPC(3)+PTC(3));
    
    text(axH, textPosX, textPosY-2, textPosZ, 'TTP','Rotation',0,'FontSize',16,'FontWeight','bold')
    
    anatomicalViewButtons(axH, 'RAS')
    
    % % For publication
    % axis(axH, 'off')
    % eqEllipsoid = equivalentEllipsoid(iFemur.vertices);
    % ttpNormal = planeNormal(createPlane(LPC, MPC, PTC));
    % camTar = mean(iFemur.vertices);
    % set(axH, 'CameraTarget', camTar);
    % set(axH, 'CameraPosition', camTar+ttpNormal*1000);
    % set(axH, 'CameraUpVector', normalizeVector3d(crossProduct3d(eqEllipsoid(7:9),ttpNormal)));
    % set(axH, 'CameraViewAngle',15)
    % set(figH, 'GraphicsSmoothing','off')
    % export_fig('Figure5', '-tif', '-r300')
    
    Tabletop(femur, side, FHC, LM_Idx, 'visu', visu);
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