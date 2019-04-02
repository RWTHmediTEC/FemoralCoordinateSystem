function [TFM, LMIdx] = Tabletop(femur, side, HJC, LMIdx, varargin)

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
ShaftAxis = createLine3d(femur.vertices(LMIdx.ShaftAxis(1),:),femur.vertices(LMIdx.ShaftAxis(2),:));

%% Construction of P1
% P1
[~, P1, ~] = distanceLines3d(NeckAxis, ShaftAxis);
FemoralMidLine=createLine3d(ICN, P1);

%% inital transformation
% Connection of the most posterior points of the condyles
PosteriorCondyleAxis = createLine3d(MPC, LPC);

Z = normalizeVector3d(FemoralMidLine(4:6));
Y = normalizeVector3d(crossProduct3d(FemoralMidLine(4:6), PosteriorCondyleAxis(4:6)));
X = normalizeVector3d(crossProduct3d(Y, Z));
iTFM = inv([[inv([X; Y; Z]), HJC']; [0 0 0 1]]);

if strcmp(side, 'L')
    iTFM=createRotationOz(pi)*iTFM; %#ok<MINV>
end

%% refinement 1
% transform the mesh by the inital TFM
iMesh=transformPoint3d(femur, iTFM);
% get the length of the femur
iLength = abs(max(iMesh.vertices(:,3)))+abs(min(iMesh.vertices(:,3)));
% cut off the distal part
DISTAL_FACTOR = 1/6;
distalPlane=[0 0 DISTAL_FACTOR*iLength+min(iMesh.vertices(:,3)), 1 0 0, 0 1 0];
distalPart = cutMeshByPlane(iMesh, distalPlane, 'part','below');
% cut the distal part into the medial and lateral condyle
sagittalPlane=createPlane(transformPoint3d(ICN, iTFM), [1 0 0]);
[LCMesh, ~, MCMesh]  = cutMeshByPlane(distalPart, sagittalPlane);

% start refinement: find the most posterior points of the condyles and
% rotate them into the new posterior condyle line
[tempROT, MPC, LPC] = refinePosteriorCondyleAxis(MCMesh, LCMesh);
ref1ROT = tempROT;
while ~isequal(eye(4), tempROT)
    MCMesh.vertices=transformPoint3d(MCMesh.vertices, tempROT);
    LCMesh.vertices=transformPoint3d(LCMesh.vertices, tempROT);
    [tempROT, MPC, LPC] = refinePosteriorCondyleAxis(MCMesh, LCMesh);
    ref1ROT = tempROT*ref1ROT;
end

%% refinement 2
% cut off the proximal part
PROXIMAL_FACTOR = 2/6;
proximalPlane=[0 0 -PROXIMAL_FACTOR*iLength, 1 0 0, 0 1 0];
proximalPart = cutMeshByPlane(iMesh, proximalPlane, 'part','above');
% cut off the head
iNeckAxis = transformLine3d(NeckAxis, iTFM);
if iNeckAxis(6)>0; iNeckAxis(4:6)=-iNeckAxis(4:6); end
neckPlane = createPlane(iNeckAxis(1:3), iNeckAxis(4:6));
proximalPart = cutMeshByPlane(proximalPart, neckPlane, 'part','above');
% transform the proximal part by the 1st refinement ROT
proximalPart.vertices = transformPoint3d(proximalPart.vertices, ref1ROT);
% Get the most posterior point of the trochanteric crest
[~, PTC_Idx] = min(proximalPart.vertices(:,2));
PTC=proximalPart.vertices(PTC_Idx,:);
% Create the tabletop plane
tabletopPlane = createPlane(MPC, LPC, PTC);
tabletopNormal = planeNormal(tabletopPlane);
if tabletopNormal(2) < 0; tabletopNormal=-tabletopNormal; end
% Create the 2nd refinement ROT
X = [1 0 0];
Y = normalizeVector3d(tabletopNormal);
Z = normalizeVector3d(crossProduct3d(X, Y));
ref2ROT = [[X; Y; Z; 0 0 0], [0 0 0 1]'];

% combination
TFM=ref2ROT*ref1ROT*iTFM;

% The femur in the AFCS
femurCS = transformPoint3d(femur, TFM);

% Position of the landmarks in the bone CS
MPC = transformPoint3d(MPC,ref2ROT);
LPC = transformPoint3d(LPC,ref2ROT);
PTC = transformPoint3d(PTC,ref2ROT);
% Get the index of the tabletop points
MPC_Idx = find(ismembertol(femurCS.vertices, MPC,'ByRows',true'));
LPC_Idx = find(ismembertol(femurCS.vertices, LPC,'ByRows',true'));
PTC_Idx = find(ismembertol(femurCS.vertices, PTC,'ByRows',true'));

if side == 'L' % Switch MPC & LPC
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
    visualizeMeshes(femurCS, patchProps)
    
    % Coordinate system
    Q.C = [1 0 0; 0 1 0; 0 0 1];
    QDScaling = distancePoints3d(MPC, LPC);
    Q.P = repmat([0, 0, 0], 3, 1);
    Q.D = QDScaling*[1 0 0; 0 1 0; 0 0 1];
    [~] = quiver3D(Q.P, Q.D, Q.C);
    
    % Landmarks
    drawPoint3d(femurCS.vertices(PTC_Idx,:),'MarkerFaceColor','k','MarkerEdgeColor','k')
    % Axes
    edgeProps.LineStyle='-';
    edgeProps.LineWidth = 2;
    edgeProps.Marker='o';
    edgeProps.MarkerFaceColor='k';
    edgeProps.MarkerEdgeColor='k';
    edgeProps.FaceColor='none';
    edgeProps.EdgeColor='k';
    
    tablePatch.vertices=[...
        femurCS.vertices(MPC_Idx,:);...
        femurCS.vertices(LPC_Idx,:);...
        femurCS.vertices(PTC_Idx,:)];
    tablePatch.faces=[1 2 3];
    
    patch(tablePatch, edgeProps)
    
    drawLine3d(transformLine3d(NeckAxis, TFM))
    
    text(tablePatch.vertices(:,1),tablePatch.vertices(:,2),tablePatch.vertices(:,3),...
        {'MPC';'LPC';'PTC'})
    
    medicalViewButtons('RAS')
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