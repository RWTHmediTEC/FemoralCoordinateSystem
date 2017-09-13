function [TFM, MPC_Idx, LPC_Idx] = Bergmann2016(femur, side, ...
    HJC, MPC, LPC, ICN, NeckAxis_Idx, Shaft, varargin)

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
% Shaft: fit ellipsoid to the shaft
ShaftEllipsoid = inertiaEllipsoid(Shaft);
ShaftVector = transformVector3d([1 0 0], eulerAnglesToRotation3d(ShaftEllipsoid(7:9)));
ShaftAxis = [ShaftEllipsoid(1:3) ShaftVector];
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
    iTFM=createRotationOz(pi)*iTFM;
end

%% refinement
% transform the mesh by the inertial TFM
iMesh.faces=femur.faces;
iMesh.vertices=transformPoint3d(femur.vertices, iTFM);
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
femurCS.vertices = transformPoint3d(femur.vertices, TFM);
femurCS.faces = femur.faces;

% Get the index of the most posterior point of the condyle
MPC_Idx = find(ismembertol(femurCS.vertices, MPC,'ByRows',true'));
LPC_Idx = find(ismembertol(femurCS.vertices, LPC,'ByRows',true'));

%% visualization
if visu
    % New figure
    monitorsPosition = get(0,'MonitorPositions');
    FigHandle = figure('Units','pixels','renderer','opengl', ...
        'Color', 'w','ToolBar','figure',...
        'WindowScrollWheelFcn',@zoomWithWheel,...
        'WindowButtonDownFcn',@rotateWithLeftMouse);
    if     size(monitorsPosition,1) == 1
        set(FigHandle,'OuterPosition',monitorsPosition(1,:));
    elseif size(monitorsPosition,1) == 2
        set(FigHandle,'OuterPosition',monitorsPosition(2,:));
    end
    hold on
    title({'The femur in the femoral coordinate system (Bergmann2016)';...
        'Left mouse - Rotate | Mouse wheel - Zoom'})
    cameratoolbar('SetCoordSys','none')
    axis equal; axis on; 
    xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]');
    lightH(1) = light; light('Position', -1*(get(lightH(1),'Position')));
    view(160, 15)
    
    % Coordinate system
    Q.C = [1 0 0; 0 1 0; 0 0 1];
    QDScaling = distancePoints3d(MPC, LPC);
    Q.P = repmat([0, 0, 0], 3, 1);
    Q.D = QDScaling*[1 0 0; 0 1 0; 0 0 1];
    [~] = quiver3D(Q.P, Q.D, Q.C);
    
    % Patch properties
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [223, 206, 161]/255;
    patchProps.FaceAlpha = 0.75;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    patch(femurCS, patchProps)
    
    % Landmarks
    drawPoint3d(transformPoint3d(P1, TFM),'MarkerFaceColor','k','MarkerEdgeColor','k')
    % Axes
    edgeProps.LineStyle='-';
    edgeProps.Color='k';
    edgeProps.Marker='o';
    edgeProps.MarkerFaceColor='k';
    edgeProps.MarkerEdgeColor='k';
    
    drawLine3d(transformLine3d(FemoralMidLine, TFM),'k')
    drawEdge3d(MPC,LPC, edgeProps)
    drawEdge3d(...
        femurCS.vertices(NeckAxis_Idx(1),:),...
        femurCS.vertices(NeckAxis_Idx(2),:), edgeProps);
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


function zoomWithWheel(~,evnt)
if evnt.VerticalScrollCount > 0
    CVA_old = get(gca,'CameraViewAngle');
    CVA_new = CVA_old + 1;
    draw
elseif evnt.VerticalScrollCount < 0
    CVA_old = get(gca,'CameraViewAngle');
    CVA_new = CVA_old - 1;
    draw
end
    function draw
        set(gca,'CameraViewAngle',CVA_new)
        drawnow
    end
end

function rotateWithLeftMouse(src,~)
if strcmp(get(src,'SelectionType'),'normal')
    cameratoolbar('SetMode','orbit')
end
end