function TFM  = MediTEC(femur, side, HJC, LMIdx, varargin)

% Z axis: The mechanical axis defined by intercondylar notch and the hip
%         joint center
% X axis: Orthogonal to the normal of the tabletop plane and Z axis
%         Definition of the tabletop plane:
%           - Resection of the femoral neck including the head
%           - Positioning of the posterior side of the femur on a table
%           - The three contact points define the tabletop plane
% Y axis: Orthogonal to the Z and X axis

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
MechanicalAxis = createLine3d(ICN, HJC);
switch side
    case 'R'
        tabletopPlane = createPlane(MPC, PTC, LPC);
    case 'L'
        tabletopPlane = createPlane(MPC, LPC, PTC);
    otherwise
        error('Invalid side identifier!')
end
tabletopNormal = planeNormal(tabletopPlane);

Z = normalizeVector3d(MechanicalAxis(4:6));
X = normalizeVector3d(crossProduct3d(tabletopNormal, Z));
Y = normalizeVector3d(crossProduct3d(Z, X));

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
    patchProps.FaceAlpha = 0.75;
    [~, axH] = visualizeMeshes(femurCS, patchProps);
    drawAxis3d(axH, 35,1.5)
    
    % Landmarks
    LM_Idx = struct2cell(LMIdx);
    LM_Idx = cell2mat(LM_Idx(structfun(@(x) length(x) == 1, LMIdx)));
    drawPoint3d(axH, femurCS.vertices(LM_Idx, :),...
        'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
    drawPoint3d(axH, transformPoint3d([ICN; HJC], TFM), ...
        'MarkerFaceColor','k','MarkerEdgeColor','k')
    
    % Axes
    drawLine3d(axH, transformLine3d(MechanicalAxis, TFM),'k')
    
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
    
    patch(axH, tablePatch, patchProps)
    
    textPosX=1/3*(MPC(1)+LPC(1)+PTC(1));
    textPosY=MPC(2);
    textPosZ=1/3*(MPC(3)+LPC(3)+PTC(3));
    text(axH, textPosX,textPosY,textPosZ,'Tabletop plane','Rotation',90)
    
    text(axH, tablePatch.vertices(:,1),tablePatch.vertices(:,2),tablePatch.vertices(:,3),...
        {'MPC';'LPC';'PTC'})
    
    anatomicalViewButtons(axH, 'RAS')
end

end