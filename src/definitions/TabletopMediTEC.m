function [TFM, LMIdx] = TabletopMediTEC(femur, side, HJC, LMIdx, varargin)

% Y axis: Normal of the tabletop plane.
%         Definition of the tabletop plane:
%           - Resection of the femoral neck including the head
%           - Positioning of the posterior side of the femur on a table
%           - The three contact points define the tabletop plane
% Z axis: Projection of the mechanical axis on the tabletop plane. The
%         mechanical axis is defined by intercondylar notch and the hip
%         joint center
% X axis: Orthogonal to Y and Z axis

% inputs
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
NeckAxis = createLine3d(femur.vertices(LMIdx.NeckAxis(1),:),femur.vertices(LMIdx.NeckAxis(2),:));
NeckOrthogonal = createLine3d(femur.vertices(LMIdx.NeckOrthogonal(1),:),femur.vertices(LMIdx.NeckOrthogonal(2),:));
[~, NeckAxis(1:3), ~] = distanceLines3d(NeckAxis, NeckOrthogonal);

%% Construction of an inital system
MechanicalAxis=createLine3d(ICN, HJC);

%% Inital transformation (iTFM)
% Connection of the most posterior points of the condyles
PosteriorCondyleAxis = createLine3d(MPC, LPC);

Z = normalizeVector3d(MechanicalAxis(4:6));
Y = normalizeVector3d(crossProduct3d(MechanicalAxis(4:6), PosteriorCondyleAxis(4:6)));
X = normalizeVector3d(crossProduct3d(Y, Z));
iTFM = inv([[inv([X; Y; Z]), HJC']; [0 0 0 1]]);

% If it is a left femur, rotate 180° around the Z axis
if strcmp(side, 'L')
    iTFM=createRotationOz(pi)*iTFM; %#ok<MINV>
end

%% Refinement 1
% Transform the mesh by the inital TFM
iMesh=transformPoint3d(femur, iTFM);
% Get the length of the femur
iLength = abs(max(iMesh.vertices(:,3)))+abs(min(iMesh.vertices(:,3)));
% Cut off the distal part
DISTAL_FACTOR = 1/6;
distalPlane=[0 0 DISTAL_FACTOR*iLength+min(iMesh.vertices(:,3)), 1 0 0, 0 1 0];
distalPart = cutMeshByPlane(iMesh, distalPlane, 'part','below');
% Cut the distal part into the medial and lateral condyle
sagittalPlane=createPlane(transformPoint3d(ICN, iTFM), [1 0 0]);
[LCMesh, ~, MCMesh]  = cutMeshByPlane(distalPart, sagittalPlane);

% Start refinement: find the most posterior points of the condyles and
% rotate them into the new posterior condyle line
[tempROT, MPC, LPC] = refinePosteriorCondyleAxis(MCMesh, LCMesh);
ref1ROT = tempROT;
while ~isequal(eye(4), tempROT)
    MCMesh.vertices=transformPoint3d(MCMesh.vertices, tempROT);
    LCMesh.vertices=transformPoint3d(LCMesh.vertices, tempROT);
    [tempROT, MPC, LPC] = refinePosteriorCondyleAxis(MCMesh, LCMesh);
    ref1ROT = tempROT*ref1ROT;
end

%% Refinement 2
% Cut off the proximal part
PROXIMAL_FACTOR = 2/6;
proximalPlane=[0 0 -PROXIMAL_FACTOR*iLength, 1 0 0, 0 1 0];
proximalPart = cutMeshByPlane(iMesh, proximalPlane, 'part','above');
% Resect the neck and the head
iNeckAxis = transformLine3d(NeckAxis, iTFM);
if iNeckAxis(6)>0; iNeckAxis(4:6)=-iNeckAxis(4:6); end
neckPlane = createPlane(iNeckAxis(1:3), iNeckAxis(4:6));
proximalPart = cutMeshByPlane(proximalPart, neckPlane, 'part','above');
% Transform the proximal part by the 1st refinement ROT
proximalPart.vertices = transformPoint3d(proximalPart.vertices, ref1ROT);
% Get the most posterior point of the trochanteric crest
[~, PTC_Idx] = min(proximalPart.vertices(:,2));
PTC=proximalPart.vertices(PTC_Idx,:);
% Create the tabletop plane
tabletopPlane = createPlane(MPC, LPC, PTC);
tabletopNormal = planeNormal(tabletopPlane);
% Normal has to point in anterior direction
if tabletopNormal(2) < 0; tabletopNormal=-tabletopNormal; end
% Create the 2nd refinement ROT
X = [1 0 0];
Y = normalizeVector3d(tabletopNormal);
Z = normalizeVector3d(crossProduct3d(X, Y));
ref2ROT = [[X; Y; Z; 0 0 0], [0 0 0 1]'];

% Combine all transformations
TFM=ref2ROT*ref1ROT*iTFM;

% The femur in the AFCS
femurCS = transformPoint3d(femur, TFM);
% Sanity check:
% Projection of the mechanical axis on the tabletop plane == Z axis?
MAproj=projLineOnPlane(transformLine3d(MechanicalAxis,TFM),[0 0 0 1 0 0 0 0 1]);
assert(all(ismembertol([0 0 1],normalizeVector3d(MAproj(4:6)))))

% Position of the landmarks in the bone CS
MPC = transformPoint3d(MPC,ref2ROT);
LPC = transformPoint3d(LPC,ref2ROT);
PTC = transformPoint3d(PTC,ref2ROT);
% Get the vertex indices of the tabletop points
MPC_Idx = find(ismembertol(femurCS.vertices, MPC,'ByRows',true'));
LPC_Idx = find(ismembertol(femurCS.vertices, LPC,'ByRows',true'));
PTC_Idx = find(ismembertol(femurCS.vertices, PTC,'ByRows',true'));

if strcmp(side, 'L') % Switch MPC & LPC
    MPC_Idx(2)=LPC_Idx; LPC_Idx=MPC_Idx(1); MPC_Idx(1)=[];
end

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
    [~,axH]=visualizeMeshes(femurCS, patchProps);
    drawAxis3d(axH,35,1.5)
    
    % Landmarks
    drawPoint3d(femurCS.vertices(PTC_Idx,:),'MarkerFaceColor','k','MarkerEdgeColor','k')
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
    
    text(textPosX,textPosY,textPosZ,'Tabletop plane','Rotation',90)
    
    text(tablePatch.vertices(:,1),tablePatch.vertices(:,2),tablePatch.vertices(:,3),...
        {'MPC';'LPC';'PTC'})
    
    anatomicalViewButtons('RAS')
end

end

function [ROT, MPC, LPC] = refinePosteriorCondyleAxis(MCMesh, LCMesh)

% Get the index of the most posterior point of the condyle
[~, MPC_Idx] = min(MCMesh.vertices(:,2));
[~, LPC_Idx] = min(LCMesh.vertices(:,2));
% Get the most posterior point of the condyle
MPC=MCMesh.vertices(MPC_Idx,:);
LPC=LCMesh.vertices(LPC_Idx,:);
% Construct the rotation into the most posterior points
PosteriorCondyleAxis=createLine3d(MPC, LPC);
Z = [0 0 1];
Y = normalizeVector3d(crossProduct3d(Z, PosteriorCondyleAxis(4:6)));
X = normalizeVector3d(crossProduct3d(Y, Z));
ROT = [[X; Y; Z; 0 0 0], [0 0 0 1]'];

end