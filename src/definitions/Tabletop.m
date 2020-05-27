function TFM = Tabletop(femur, side, HJC, LMIdx, varargin)

% Femoral bone CS based on the tabletop plane:
%   - X axis: Axis through the medial and lateral posterior condyle 
%   - Y axis: Normal of the tabletop plane
%   - Z axis: Orthogonal to Y and Z axis

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
PTC = femur.vertices(LMIdx.PosteriorTrochantericCrest,:);

%% Transformation into bone CS
switch side
    case 'R'
        PosteriorCondyleAxis = createLine3d(MPC, LPC);
        tabletopPlane = createPlane(MPC, PTC, LPC);
    case 'L'
        PosteriorCondyleAxis = createLine3d(LPC, MPC);
        tabletopPlane = createPlane(MPC, LPC, PTC);
    otherwise
        error('Invalid side identifier!')
end
tabletopNormal = planeNormal(tabletopPlane);
X = normalizeVector3d(PosteriorCondyleAxis(4:6));
Y = normalizeVector3d(tabletopNormal);
Z = normalizeVector3d(crossProduct3d(X, Y));

TFM = [[X; Y; Z],[0 0 0]'; [0 0 0 1]]*createTranslation3d(-HJC);

%% visualization
if visu
    % The femur in the AFCS
    femurCS = transformPoint3d(femur, TFM);
    
    % Position of the landmarks in the bone CS
    MPC = transformPoint3d(MPC,TFM);
    LPC = transformPoint3d(LPC,TFM);
    PTC = transformPoint3d(PTC,TFM);
    
    % Patch properties
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [223, 206, 161]/255;
    patchProps.FaceAlpha = 0.75;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    [~, axH] = visualizeMeshes(femurCS, patchProps);
    drawAxis3d(axH, 35,1.5)
    
    % Landmarks
    LM_Idx = struct2cell(LMIdx);
    LM_Idx = cell2mat(LM_Idx(structfun(@(x) length(x) == 1, LMIdx)));
    drawPoint3d(axH, femurCS.vertices(LM_Idx, :),...
        'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
    
    % Tabletop patch
    patchProps.LineStyle='-';
    patchProps.Marker='o';
    patchProps.MarkerFaceColor='k';
    patchProps.MarkerEdgeColor='k';
    patchProps.FaceColor='k';
    patchProps.FaceAlpha=0.5;
    patchProps.EdgeColor='k';
    
    tablePatch.vertices=[MPC; LPC; PTC];
    tablePatch.faces=[1 2 3];
    
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