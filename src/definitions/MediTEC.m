function TFM  = MediTEC(femur, side, HJC, LMIdx, varargin)

% Z axis: The mechanical axis defined by intercondylar notch and the hip
%         joint center
% Y axis: Orthogonal to the Z axis and posterior condylar axis
% X axis: Orthogonal to the Y and Z axis

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

%% Axes
mechanicalAxis = createLine3d(ICN, HJC);
switch side
    case 'R'
        posteriorCondylarAxis = createLine3d(MPC, LPC);
    case 'L'
        posteriorCondylarAxis = createLine3d(LPC, MPC);
    otherwise
        error('Invalid side identifier!')
end

Z = normalizeVector3d(mechanicalAxis(4:6));
Y = normalizeVector3d(crossProduct3d(Z, posteriorCondylarAxis(4:6)));
X = normalizeVector3d(crossProduct3d(Y, Z));

TFM = [[X; Y; Z],[0 0 0]'; [0 0 0 1]]*createTranslation3d(-HJC);

%% visualization
if visu
    % The femur in the AFCS
    femurCS = transformPoint3d(femur, TFM);
    
    % Position of the landmarks in the bone CS
    ICN = transformPoint3d(ICN, TFM);
    HJC = transformPoint3d(HJC, TFM);
    MPC = transformPoint3d(MPC, TFM);
    LPC = transformPoint3d(LPC, TFM);

    
    % Patch properties
    patchProps.FaceAlpha = 0.75;
    [~, axH] = visualizeMeshes(femurCS, patchProps);
    drawAxis3d(axH, 35,1.5)
    
    % Landmarks
    LM_Idx = struct2cell(LMIdx);
    LM_Idx = cell2mat(LM_Idx(structfun(@(x) length(x) == 1, LMIdx)));
    drawPoint3d(axH, femurCS.vertices(LM_Idx, :),...
        'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
    drawPoint3d(axH, [ICN; HJC; MPC; LPC], ...
        'MarkerFaceColor','k','MarkerEdgeColor','k')
    
    % Axes
    drawLine3d(axH, transformLine3d(mechanicalAxis, TFM),'k')
    drawLine3d(axH, transformLine3d(posteriorCondylarAxis, TFM),'k')
    
    drawLabels3d(axH, [MPC; LPC; ICN],{'MPC';'LPC';'ICN'})
    
    anatomicalViewButtons(axH, 'RAS')
end

end