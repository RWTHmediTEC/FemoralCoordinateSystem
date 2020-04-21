function [TFM, LMIdx] = Bergmann2016(femur, side, HJC, LMIdx, varargin)

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

% Landmarks
MPC = femur.vertices(LMIdx.MedialPosteriorCondyle,:);
LPC = femur.vertices(LMIdx.LateralPosteriorCondyle,:);
ICN = femur.vertices(LMIdx.IntercondylarNotch,:);
NeckAxis = createLine3d(femur.vertices(LMIdx.NeckAxis(1),:),femur.vertices(LMIdx.NeckAxis(2),:));
ShaftAxis = createLine3d(femur.vertices(LMIdx.ShaftAxis(1),:),femur.vertices(LMIdx.ShaftAxis(2),:));

%% Construction of P1
[~, P1, ~] = distanceLines3d(NeckAxis, ShaftAxis);
StraightFemurAxis=createLine3d(ICN, P1);

%% inital transformation
% Connection of the most posterior points of the condyles
PosteriorCondyleAxis = createLine3d(MPC, LPC);

Z = normalizeVector3d(StraightFemurAxis(4:6));
Y = normalizeVector3d(crossProduct3d(StraightFemurAxis(4:6), PosteriorCondyleAxis(4:6)));
X = normalizeVector3d(crossProduct3d(Y, Z));
iTFM = inv([[inv([X; Y; Z]), HJC']; [0 0 0 1]]);

if strcmp(side, 'L')
    iTFM=createRotationOz(pi)*iTFM; %#ok<MINV>
end

%% refinement
% transform the mesh by the initial TFM
iMesh=transformPoint3d(femur, iTFM);
% get the length of the femur
iLength = abs(max(iMesh.vertices(:,3)))+abs(min(iMesh.vertices(:,3)));
% cut off the distal part
DISTAL_FACTOR = 1/6;
distalPlane=[0 0 DISTAL_FACTOR*iLength+min(iMesh.vertices(:,3)), 1 0 0, 0 1 0];
distalPart = cutMeshByPlane(iMesh, distalPlane, 'part','below');
% cut the distal part into the medial and lateral condyle
sagittalPlane=createPlane(transformPoint3d(ICN, iTFM), [1 0 0]);
[LCMesh, ~, MCMesh] = cutMeshByPlane(distalPart, sagittalPlane);

% start refinement: find the most posterior points of the condyles and
% rotate them into the new posterior condyle line
[tempRot, MPC, LPC] = refinePosteriorCondyleAxis(MCMesh, LCMesh);
refRot = tempRot;
while ~isequal(eye(4), tempRot)
    MCMesh.vertices=transformPoint3d(MCMesh.vertices, tempRot);
    LCMesh.vertices=transformPoint3d(LCMesh.vertices, tempRot);
    [tempRot, MPC, LPC] = refinePosteriorCondyleAxis(MCMesh, LCMesh);
    refRot = tempRot*refRot;
end

% combination
TFM=refRot*iTFM;

% The femur in the AFCS
femurCS = transformPoint3d(femur, TFM);

% Get the index of the most posterior point of the condyle
switch side
    case 'R'
        LMIdx.MedialPosteriorCondyle = find(ismembertol(femurCS.vertices, MPC,'ByRows',true'));
        LMIdx.LateralPosteriorCondyle = find(ismembertol(femurCS.vertices, LPC,'ByRows',true'));
    case 'L'
        LMIdx.MedialPosteriorCondyle = find(ismembertol(femurCS.vertices, LPC,'ByRows',true'));
        LMIdx.LateralPosteriorCondyle = find(ismembertol(femurCS.vertices, MPC,'ByRows',true'));
end

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
    drawPoint3d(transformPoint3d(P1, TFM),'MarkerFaceColor','k','MarkerEdgeColor','k')
    % Axes
    edgeProps.LineStyle='-';
    edgeProps.Color='k';
    
    drawEdge3d(axH,clipLine3d(transformLine3d(StraightFemurAxis, TFM),...
        boundingBox3d(femurCS.vertices)),edgeProps)
    drawEdge3d(axH,...
        femurCS.vertices(LMIdx.NeckAxis(1),:),...
        femurCS.vertices(LMIdx.NeckAxis(2),:), edgeProps);
    edgeProps.Marker='o';
    edgeProps.MarkerEdgeColor='k';
    edgeProps.MarkerFaceColor='k';
    
    PC=[femurCS.vertices(LMIdx.MedialPosteriorCondyle,:);...
        femurCS.vertices(LMIdx.LateralPosteriorCondyle,:)];
    drawEdge3d(axH,PC(1,:),PC(2,:), edgeProps)
    
    text(axH,PC(:,1),PC(:,2),PC(:,3),{'MPC';'LPC'})
    
    anatomicalViewButtons(axH,'RAS')
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