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
    % Patch properties
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [223, 206, 161]/255;
    patchProps.FaceAlpha = 0.75;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    % The femur in the AFCS
    femurCS = transformPoint3d(femur, TFM);
    visualizeMeshes(femurCS, patchProps)
    
    % Coordinate system
    Q.C = [1 0 0; 0 1 0; 0 0 1];
    QDScaling = distancePoints3d(MEC, LEC);
    Q.P = repmat([0, 0, 0], 3, 1);
    Q.D = QDScaling*[1 0 0; 0 1 0; 0 0 1];
    [~] = quiver3D(Q.P, Q.D, Q.C);
    
    % Landmarks
    drawPoint3d(transformPoint3d(MEC_LEC_midPoint, TFM),...
        'MarkerFaceColor','k','MarkerEdgeColor','k')
    drawLine3d(transformLine3d(MechanicalAxis, TFM),'k')
    edgeProps.LineStyle='-';
    edgeProps.Color='k';
    edgeProps.Marker='o';
    edgeProps.MarkerEdgeColor='k';
    edgeProps.MarkerFaceColor='k';
    drawEdge3d([transformPoint3d(MEC,TFM),transformPoint3d(LEC,TFM)],edgeProps)
    
    viewButtonsASR
end
end