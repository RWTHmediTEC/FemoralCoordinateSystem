function [TFM, MPC_Idx, LPC_Idx] = Bergmann2016(femur, side, ...
    HJC, MPC, LPC, ICN, NeckAxis_Idx, ShaftAxis, varargin)

% 2016 - Bergmann et al. - Standardized Loads Acting in Hip Implants:
%   "The origin of this coordinate system is located at the centre of the 
%   femoral head. The +z axis points upward and is defined by the line 
%   connecting the two points where the curved femoral mid-line intersected 
%   with the neck axis (P1) and where it passes the intercondylar notch 
%   (P2). The +x axis points laterally and is oriented parallel to the 
%   proximal contour of the condyles. The +y axis points in the anterior 
%   direction."

% inputs
p = inputParser;
addRequired(p,'femur',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addRequired(p,'side',@(x) any(validatestring(x,{'R','L'})));
addOptional(p,'visualization',true,@islogical);
parse(p,femur,side,varargin{:});

femur = p.Results.femur;
visu = p.Results.visualization;

%% Construction of P1
NeckAxis = createLine3d(femur.vertices(NeckAxis_Idx(1),:),femur.vertices(NeckAxis_Idx(2),:));
% P1
[~, P1, ~] = distanceLines3d(NeckAxis, ShaftAxis);
FemoralMidLine=createLine3d(ICN, P1);

%% inital transformation
% Connection of the most posterior points of the condyles
PosteriorCondyleAxis = createLine3d(MPC, LPC);

Z = normalizeVector3d(FemoralMidLine(4:6));
Y = normalizeVector3d(vectorCross3d(FemoralMidLine(4:6), PosteriorCondyleAxis(4:6)));
X = normalizeVector3d(vectorCross3d(Y, Z));
iTFM = inv([[inv([X; Y; Z]), HJC']; [0 0 0 1]]);

if strcmp(side, 'L')
    iTFM=createRotationOz(pi)*iTFM; %#ok<MINV>
end

%% refinement
% transform the mesh by the inertial TFM
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
    drawPoint3d(transformPoint3d(P1, TFM),'MarkerFaceColor','k','MarkerEdgeColor','k')
    % Axes
    edgeProps.LineStyle='-';
    edgeProps.Color='k';
    
    drawEdge3d(clipLine3d(transformLine3d(FemoralMidLine, TFM),...
        boundingBox3d(femurCS.vertices)),edgeProps)
    drawEdge3d(...
        femurCS.vertices(NeckAxis_Idx(1),:),...
        femurCS.vertices(NeckAxis_Idx(2),:), edgeProps);
    edgeProps.Marker='o';
    edgeProps.MarkerEdgeColor='k';
    edgeProps.MarkerFaceColor='k';
    drawEdge3d(...
        femurCS.vertices(MPC_Idx,:),...
        femurCS.vertices(LPC_Idx,:), edgeProps)
    
    viewButtonsRAS
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
Y = normalizeVector3d(vectorCross3d(Z, PosteriorCondyleAxis(4:6)));
X = normalizeVector3d(vectorCross3d(Y, Z));
ROT = [[X; Y; Z; 0 0 0], [0 0 0 1]'];

end