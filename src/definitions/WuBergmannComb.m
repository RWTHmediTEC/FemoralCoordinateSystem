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

%% visualization
if visu
    % The femur in the AFCS
    femurCS = transformPoint3d(femur, TFM);
    
    % Patch properties
    patchProps.FaceAlpha = 0.75;
    [~, axH] = visualizeMeshes(femurCS, patchProps);
    drawAxis3d(axH, 35,1.5)
    
    % Landmarks
    LM_Idx = struct2cell(LMIdx);
    LM_Idx = cell2mat(LM_Idx(structfun(@(x) length(x) == 1, LMIdx)));
    drawPoint3d(axH, femurCS.vertices(LM_Idx, :),...
        'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
    drawPoint3d(axH, transformPoint3d(ICN, TFM),'MarkerFaceColor','k','MarkerEdgeColor','k')
    % Axes
    drawLine3d(axH, transformLine3d(MechanicalAxis, TFM),'k')
    edgeProps.LineStyle='-';
    edgeProps.Color='k';
    edgeProps.Marker='o';
    edgeProps.MarkerEdgeColor='k';
    edgeProps.MarkerFaceColor='k';
    
    PC=[femurCS.vertices(LMIdx.MedialPosteriorCondyle,:);...
        femurCS.vertices(LMIdx.LateralPosteriorCondyle,:)];
    drawEdge3d(axH, PC(1,:),PC(2,:), edgeProps)
    
    text(axH, PC(:,1),PC(:,2),PC(:,3),{'MPC';'LPC'})
    
    anatomicalViewButtons(axH, 'ASR')
end

end