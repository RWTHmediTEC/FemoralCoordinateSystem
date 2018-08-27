function TFM = Wu2002(femur, side, HJC, LMIdx, varargin)

% 2002 - Wu et al. - ISB recommendation on definitions of joint coordinate 
% systems of various joints for the reporting of human joint motion Part 1
% Femoral coordinate system—xyz:
% o: The origin coincident with the right (or left) hip center of rotation, 
%    coincident with that of the pelvic coordinate system (o) in the 
%    neutral configuration.
% y: The line joining the midpoint between the medial and lateral femoral 
%    epicondyles (FEs) and the origin, and pointing cranially.
% z: The line perpendicular to the y-axis, lying in the plane defined by 
%    the origin and the two FEs, pointing to the right.
% x: The line perpendicular to both y- and z-axis, pointing anteriorly 
% (Cappozzo et al., 1995)

% inputs
p = inputParser;
addRequired(p,'femur',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addRequired(p,'side',@(x) any(validatestring(x,{'R','L'})));
addOptional(p,'visualization',true,@islogical);
parse(p,femur,side,varargin{:});

femur = p.Results.femur;
visu = p.Results.visualization;

%% Landmarks
MEC = femur.vertices(LMIdx.MedialEpicondyle,:);
LEC = femur.vertices(LMIdx.LateralEpicondyle,:);

%% Coordinate system
% Midpoint between the epicondyles
MEC_LEC_midPoint=midPoint3d(MEC, LEC);
% Mechanical axis is the connection of EC midpoint and hip joint center
MechanicalAxis = createLine3d(MEC_LEC_midPoint, HJC);
% Connection of the epicondyles
EpicondyleAxis = createLine3d(MEC, LEC);

Y = normalizeVector3d(MechanicalAxis(4:6));
X = normalizeVector3d(crossProduct3d(MechanicalAxis(4:6), EpicondyleAxis(4:6)));
Z = normalizeVector3d(crossProduct3d(X, Y));

TFM = inv([[inv([X; Y; Z]), HJC']; [0 0 0 1]]);

if strcmp(side, 'L'); TFM=createRotationOy(pi)*TFM; end %#ok<MINV>

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
    EC_TFM=[transformPoint3d(MEC,TFM);transformPoint3d(LEC,TFM)];
    drawEdge3d(EC_TFM(1,:),EC_TFM(2,:),edgeProps)
    text(EC_TFM(:,1),EC_TFM(:,2),EC_TFM(:,3),{'MEC';'LEC'})
    
    medicalViewButtons('ASR')
end

end