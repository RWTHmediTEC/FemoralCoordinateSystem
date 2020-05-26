function TFM = WuBergmannComb(femur, side, HJC, LMIdx, varargin)

% A Combination of Wu2002 and Bergmann2016
%   Mechanical axis: Connection of intercondylar notch and hip joint center
%   Mediolateral direction based on the posterior condyle axis

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

%% inital transformation
% mechanical axis is the connection of intercondylar notch and hip joint center
MechanicalAxis = createLine3d(ICN, HJC);
% connection of the most posterior points of the condyles
PosteriorCondyleAxis = createLine3d(MPC, LPC);

Y = normalizeVector3d(MechanicalAxis(4:6));
X = normalizeVector3d(crossProduct3d(MechanicalAxis(4:6), PosteriorCondyleAxis(4:6)));
Z = normalizeVector3d(crossProduct3d(X, Y));

TFM = [[X; Y; Z],[0 0 0]'; [0 0 0 1]]*createTranslation3d(-HJC);
% If it is a left femur, rotate 180° around the Y axis
if strcmp(side, 'L')
    TFM=createRotationOy(pi)*TFM;
end

% The femur in the AFCS
femurCS = transformPoint3d(femur, TFM);

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
    drawPoint3d(transformPoint3d(ICN, TFM),'MarkerFaceColor','k','MarkerEdgeColor','k')
    % Axes
    drawLine3d(transformLine3d(MechanicalAxis, TFM),'k')
    edgeProps.LineStyle='-';
    edgeProps.Color='k';
    edgeProps.Marker='o';
    edgeProps.MarkerEdgeColor='k';
    edgeProps.MarkerFaceColor='k';
    
    PC=[femurCS.vertices(LMIdx.MedialPosteriorCondyle,:);...
        femurCS.vertices(LMIdx.LateralPosteriorCondyle,:)];
    drawEdge3d(PC(1,:),PC(2,:), edgeProps)
    
    text(PC(:,1),PC(:,2),PC(:,3),{'MPC';'LPC'})
    
    anatomicalViewButtons('ASR')
end

end