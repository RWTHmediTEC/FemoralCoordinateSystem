function [TFM, MPC_Idx, LPC_Idx] = WuBergmannComb(femur, side, HJC, MPC, LPC, ICN, varargin)

% A Combination of Wu2002 and Bergmann2016
%   Mechanical axis: Connection of intercondylar notch and hip joint center
%   Medial-Lateral direction based on the posterior condyle axis

% inputs
p = inputParser;
addRequired(p,'femur',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addRequired(p,'side',@(x) any(validatestring(x,{'R','L'})));
addOptional(p,'visualization',true,@islogical);
parse(p,femur,side,varargin{:});

femur = p.Results.femur;
visu = p.Results.visualization;

%% inital transformation
% mechanical axis is the connection of intercondylar notch and hip joint center
MechanicalAxis = createLine3d(ICN, HJC);
% connection of the most posterior points of the condyles
PosteriorCondyleAxis = createLine3d(MPC, LPC);

Y = normalizeVector3d(MechanicalAxis(4:6));
X = normalizeVector3d(crossProduct3d(MechanicalAxis(4:6), PosteriorCondyleAxis(4:6)));
Z = normalizeVector3d(crossProduct3d(X, Y));
iTFM = inv([[inv([X; Y; Z]), HJC']; [0 0 0 1]]);

if strcmp(side, 'L')
    iTFM=createRotationOy(pi)*iTFM;
end

%% refinement
% transform the mesh by the inertial TFM
iMesh=transformPoint3d(femur, iTFM);
% get the length of the femur
iLength = abs(max(iMesh.vertices(:,2)))+abs(min(iMesh.vertices(:,2)));
% cut off the distal part
DISTAL_FACTOR = 1/6;
distalPlane=[0 DISTAL_FACTOR*iLength+min(iMesh.vertices(:,2)) 0, 0 0 1, 1 0 0];
distalPart = cutMeshByPlane(iMesh, distalPlane, 'part','below');
% cut the distal part into the medial and lateral condyle
sagittalPlane=[0 0 0 1 0 0 0 1 0];
[LCMesh, ~, MCMesh]  = cutMeshByPlane(distalPart, sagittalPlane);

% start refinement: find the most posterior points of the condyles and
% rotate them into the new posterior condyle line
tempRot = refinePosteriorCondyleAxis(MCMesh, LCMesh);
refRot = tempRot;
while ~isequal(eye(4), tempRot)
    MCMesh.vertices=transformPoint3d(MCMesh.vertices, tempRot);
    LCMesh.vertices=transformPoint3d(LCMesh.vertices, tempRot);
    [tempRot, MPC, LPC]  = refinePosteriorCondyleAxis(MCMesh, LCMesh);
    refRot = tempRot*refRot;
end

% combination
TFM=refRot*iTFM;

% The femur in the AFCS
femurCS = transformPoint3d(femur, TFM);

% Get the index of the most posterior point of the condyle
MPC_Idx = find(ismembertol(femurCS.vertices, MPC,'ByRows',true'));
LPC_Idx = find(ismembertol(femurCS.vertices, LPC,'ByRows',true'));

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
    drawPoint3d(transformPoint3d(ICN, TFM),'MarkerFaceColor','k','MarkerEdgeColor','k')
    % Axes
    drawLine3d(transformLine3d(MechanicalAxis, TFM),'k')
    edgeProps.LineStyle='-';
    edgeProps.Color='k';
    edgeProps.Marker='o';
    edgeProps.MarkerEdgeColor='k';
    edgeProps.MarkerFaceColor='k';
    drawEdge3d(MPC,LPC,edgeProps)
    
    viewButtonsASR
end
end

function [ROT, MPC, LPC] = refinePosteriorCondyleAxis(MCMesh, LCMesh)

% Get the index of the most posterior point of the condyle
[~, MPC_Idx] = min(MCMesh.vertices(:,1));
[~, LPC_Idx] = min(LCMesh.vertices(:,1));
% Get the most posterior point of the condyle
MPC=MCMesh.vertices(MPC_Idx,:);
LPC=LCMesh.vertices(LPC_Idx,:);
% Construct the rotation into the most posterior points
PosteriorCondyleAxis=createLine3d(MPC, LPC);
Y = [0 1 0];
X = normalizeVector3d(crossProduct3d(Y, PosteriorCondyleAxis(4:6)));
Z = normalizeVector3d(crossProduct3d(X, Y));
ROT = [[X; Y; Z; 0 0 0], [0 0 0 1]'];

end