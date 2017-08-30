function TFM = femoralCS_Wu2002(femur, side, HJC, MEC, LEC, varargin)

% inputs
p = inputParser;
addRequired(p,'femur',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addRequired(p,'side',@(x) any(validatestring(x,{'R','L'})));
addOptional(p,'visualization',true,@islogical);
parse(p,femur,side,varargin{:});

femur = p.Results.femur;
visu = p.Results.visualization;

% algorithm
MEC_LEC_midPoint=midPoint3d(MEC, LEC);

MechanicalAxis = createLine3d(MEC_LEC_midPoint, HJC);

CondyleAxis = createLine3d(MEC, LEC);

Y = normalizeVector3d(MechanicalAxis(4:6));

X = normalizeVector3d(vectorCross3d(MechanicalAxis(4:6), CondyleAxis(4:6)));

Z = normalizeVector3d(vectorCross3d(X, Y));

TFM = inv([[inv([X; Y; Z]), HJC']; [0 0 0 1]]);

if strcmp(side, 'L')
    TFM=createRotationOy(pi)*TFM;
end

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
    title({'The femur in the femoral coordinate system (Wu2002)';...
        'Left mouse - Rotate | Mouse wheel - Zoom'})
    cameratoolbar('SetCoordSys','none')
    axis equal; axis on; 
    xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]');
    lightH(1) = light; light('Position', -1*(get(lightH(1),'Position')));
    view(90,0)
    
    % Coordinate system
    Q.C = [1 0 0; 0 1 0; 0 0 1];
    QDScaling = distancePoints3d(MEC, LEC);
    Q.P = repmat([0, 0, 0], 3, 1);
    Q.D = QDScaling*[1 0 0; 0 1 0; 0 0 1];
    [~] = quiver3D(Q.P, Q.D, Q.C);
    
    % Patch properties
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [223, 206, 161]/255;
    patchProps.FaceAlpha = 0.75;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    
    % The pelvis in the APCS
    femurCS.vertices = transformPoint3d(femur.vertices, TFM);
    femurCS.faces = femur.faces;
    patch(femurCS, patchProps)
    
    % Landmarks
    drawPoint3d(transformPoint3d([MEC;LEC;MEC_LEC_midPoint], TFM),...
        'MarkerFaceColor','k','MarkerEdgeColor','k')
    
    drawLine3d(transformLine3d(MechanicalAxis, TFM),'k')
    drawEdge3d([transformPoint3d(MEC,TFM),transformPoint3d(LEC,TFM)])
    
    % Viewpoint
    set(gca, 'CameraUpVector',[0 1 0])
    CamNormal=get(gca, 'CameraPosition')-get(gca, 'CameraTarget');
    CamDist=vectorNorm3d(CamNormal);
    set(gca, 'CameraPosition', get(gca, 'CameraTarget')+...
        [0.95, 0.05, 0.3]*CamDist);
end
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