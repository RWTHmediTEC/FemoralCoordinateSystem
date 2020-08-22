function TFM  = MediTEC(femur, side, HJC, LMIdx, varargin)

% Z axis: The mechanical axis defined by intercondylar notch and the hip
%         joint center
% Y axis: Orthogonal to the Z axis and posterior condylar axis
% X axis: Orthogonal to the Y and Z axis

% inputs
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addRequired(p,'femur',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addRequired(p,'side',@(x) any(validatestring(x,{'R','L'})));
addOptional(p,'visualization',true,logParValidFunc);
parse(p,femur,side,varargin{:});

femur = p.Results.femur;
visu = p.Results.visualization;

%% Landmarks
if isfield(LMIdx,'MedialPosteriorCondyle')
    MPC = femur.vertices(LMIdx.MedialPosteriorCondyle,:);
elseif isfield(LMIdx,'MPC')
    MPC = femur.vertices(LMIdx.MPC,:);
else
    error('Medial posterior condyle (MPC) is missing!')
end
if isfield(LMIdx,'LateralPosteriorCondyle')
    LPC = femur.vertices(LMIdx.LateralPosteriorCondyle,:);
elseif isfield(LMIdx,'LPC')
    LPC = femur.vertices(LMIdx.LPC,:);
else
    error('Lateral posterior condyle (LPC) is missing!')
end
if isfield(LMIdx,'IntercondylarNotch')
    ICN = femur.vertices(LMIdx.IntercondylarNotch,:);
elseif isfield(LMIdx,'ICN')
    ICN = femur.vertices(LMIdx.ICN,:);
else
    error('Intercondylar notch (ICN) is missing!')
end

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
    drawPoint3d(axH, [ICN; HJC; MPC; LPC], ...
        'MarkerFaceColor','k','MarkerEdgeColor','k')
    drawLabels3d(axH, [MPC; LPC; ICN],{'MPC';'LPC';'ICN'})
    % Other landmarks
    LMIdx_cell = struct2cell(LMIdx);
    LM = cell2mat(LMIdx_cell(structfun(@(x) length(x) == 1, LMIdx)));
    drawPoint3d(axH, femurCS.vertices(LM, :),...
        'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
    
    % Axes
    drawLine3d(axH, transformLine3d(mechanicalAxis, TFM),'Color','k','LineWidth',2)
    drawLine3d(axH, transformLine3d(posteriorCondylarAxis, TFM),'Color','k','LineWidth',2)
    % % Other axes
    % Axes = cell2mat(LMIdx_cell(structfun(@(x) length(x) == 2, LMIdx)));
    % for a=1:size(Axes,1)
    %     drawLine3d(axH, createLine3d(femurCS.vertices(Axes(a,1),:),femurCS.vertices(Axes(a,2),:)),'k')
    % end
    
    anatomicalViewButtons(axH, 'RAS')
end

end