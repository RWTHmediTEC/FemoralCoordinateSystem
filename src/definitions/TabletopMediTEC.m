function TFM  = TabletopMediTEC(femur, side, HJC, LMIdx, varargin)

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
PTC = femur.vertices(LMIdx.PosteriorTrochantericCrest,:);

%% Axes
MechanicalAxis=createLine3d(ICN, HJC);
switch side
    case 'R'
        tabletopPlane = createPlane(MPC, PTC, LPC);
    case 'L'
        tabletopPlane = createPlane(MPC, LPC, PTC);
    otherwise
        error('Invalid side identifier!')
end
tabletopNormal = planeNormal(tabletopPlane);
projMechanicalAxis = projLineOnPlane(MechanicalAxis,tabletopPlane);
Z = normalizeVector3d(projMechanicalAxis(4:6));
Y = normalizeVector3d(tabletopNormal);
X = normalizeVector3d(crossProduct3d(Y, Z));

TFM = [[X; Y; Z],[0 0 0]'; [0 0 0 1]]*createTranslation3d(-HJC);

%% visualization
if visu
    % The femur in the AFCS
    femurCS = transformPoint3d(femur, TFM);
    
    % Position of the landmarks in the bone CS
    MPC = transformPoint3d(MPC, TFM);
    LPC = transformPoint3d(LPC, TFM);
    PTC = transformPoint3d(PTC, TFM);
    
    % Patch properties
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [223, 206, 161]/255;
    patchProps.FaceAlpha = 0.75;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    [~,axH]=visualizeMeshes(femurCS, patchProps);
    drawAxis3d(axH,35,1.5)
    
    % Landmarks
    drawPoint3d(transformPoint3d([ICN; HJC], TFM), ...
        'MarkerFaceColor','k','MarkerEdgeColor','k')
    
    % Axes
    drawLine3d(transformLine3d(MechanicalAxis, TFM),'k')
    
    % Tabletop patch
    patchProps.LineStyle='-';
    patchProps.Marker='o';
    patchProps.MarkerFaceColor='k';
    patchProps.MarkerEdgeColor='k';
    patchProps.FaceColor='k';
    patchProps.FaceAlpha=0.5;
    patchProps.EdgeColor='k';
    
    tablePatch.vertices=[MPC; LPC; PTC];
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