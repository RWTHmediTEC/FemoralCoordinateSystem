function [TFM2AFCS, LMIdx, FHC, LM, TFM] = automaticFemoralCS(femur, side, varargin)
%AUTOMATICFEMORALCS Calculate an femoral coordinate system (AFCS)
%
% REQUIRED INPUT:
%   femur: A "clean" mesh as struct with the fields vertices and faces
%       "Clean" means a closed outer surface without holes or tunnels,
%       normals are oriented outwards, no duplicate vertices, ...
%   side: 'L' or 'R': left or right femur
%
% OPTIONAL INPUT:
%   'FHC': 1x3 vector: Coordinates of the femoral head / hip joint center
%        in the coordinate system (CS) of the input femur mesh
%   'definition': The definition to construct the femoral CS
%       'Wu2002' - 2002 - Wu et al. - ISB recommendation on definitions
%           of joint coordinate systems of various joints for the reporting
%           of human joint motion - part I: ankle, hip, and spine (default)
%       'Bergmann2016' - 2016 - Bergmann et al. - Standardized Loads Acting
%           in Hip Implants
%       'Tabletop' - Defined by the table top plane
%       'MediTEC' - Defined by the mechanical axis and table top plane
%   'visualization': true (default) or false
%   'subject': Char: Identification of the subject. Default is '?'.
%
% OUTPUT:
%   TFM2AFCS: Transformation of the femur into the femoral CS
%   LMIdx: Landmark indices into the vertices of the femur
%   FHC: Coordinates of the femoral head center
%
% TODO:
%   - Add option to select of anatomical orientation of the femur CS
%
% AUTHOR: Maximilian C. M. Fischer
% 	mediTEC - Chair of Medical Engineering, RWTH Aachen University
% VERSION: 1.0.3
% DATE: 2018-08-24
% COPYRIGHT (C) 2020 Maximilian C. M. Fischer
% LICENSE: EUPL v1.2
% 

addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'src']));
addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'res']));

%% Parse inputs
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addRequired(p,'femur',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addRequired(p,'side',@(x) any(validatestring(upper(x(1)),{'R','L'})));
isPoint3d = @(x) validateattributes(x,{'numeric'},...
    {'nonempty','nonnan','real','finite','size',[1,3]});
addParameter(p,'FHC',nan, isPoint3d);
isLine3d = @(x) validateattributes(x,{'numeric'},...
    {'nonempty','nonnan','real','finite','size',[1,6]});
addParameter(p,'NeckAxis',nan, isLine3d);
validCSStrings={'Wu2002','Bergmann2016','Tabletop','MediTEC'};
addParameter(p,'definition','Wu2002',@(x) any(validatestring(x,validCSStrings)));
addParameter(p,'visualization',true,logParValidFunc);
addParameter(p,'verbose',false,logParValidFunc);
addParameter(p,'subject','?',@(x) validateattributes(x,{'char'},{'nonempty'}));
addParameter(p,'debugVisualization',false,logParValidFunc);

parse(p,femur,side,varargin{:});
femur = p.Results.femur;
side = upper(p.Results.side(1));
FHC = p.Results.FHC;
NeckAxis = p.Results.NeckAxis;
definition = p.Results.definition;
visu = logical(p.Results.visualization);
verb = logical(p.Results.verbose);
subject = p.Results.subject;
debugVisu = logical(p.Results.debugVisualization);

%% Algorithm
% Get inertia transformation for the input femur (subject)
femurProps = inertiaInfo(femur);
inertiaTFM = inv(femurProps.inverseInertiaTFM);

% Transform the vertices into the temporary inertia coordinate system
femurInertia = transformPoint3d(femur, inertiaTFM);

% For left femurs reflect along the x-axis after inertia transform
xReflection = eye(4);
if strcmp(side, 'L')
    xReflection(1,1) = -1;
end
femurInertia.vertices = transformPoint3d(femurInertia.vertices, xReflection);

if strcmp(side, 'L')
    % Orient the normals outwards after reflection
    femurInertia.faces = fliplr(femurInertia.faces);
end

if debugVisu
    NOP = 7;
    BW = 0.005;
    BSX = (1-(NOP+1)*BW)/NOP;
    set(0,'defaultAxesFontSize',14)
    dFigH1 = figure('Units','pixels','renderer','opengl', 'Color', 'w');
    MonitorsPos = get(0,'MonitorPositions');
    if     size(MonitorsPos,1) == 1
        set(dFigH1,'OuterPosition', MonitorsPos(1,:));
    elseif size(MonitorsPos,1) == 2
        set(dFigH1,'OuterPosition', MonitorsPos(2,:));
    end
    dFigH1.Name = [subject ': nICP registration (Debug Figure)'];
    dFigH1.NumberTitle = 'off';
    dFigH1.WindowState = 'maximized';
    
    tPH = uipanel('Title','Subject''s principal axes transf. ','FontSize',14,'BorderWidth',2,...
        'BackgroundColor','w','Position',[BW+(BSX+BW)*0 0.01 BSX 0.99]);
    dAxH1(1) = axes('Parent', tPH, 'Visible','off', 'Color','w');
    
    % The subject in the inertia CS
    visualizeMeshes(dAxH1(1), femurInertia);
    dFig1View = [90,-90];
    view(dAxH1(1),dFig1View); axis(dAxH1(1), 'tight');
end

% Load template
load('template.mat','template')

% Length of the template
templateLength = max(template.vertices(:,1))-min(template.vertices(:,1));
% Length of the subject
femurInertiaLength = max(femurInertia.vertices(:,1))-min(femurInertia.vertices(:,1));

% Scale subject in x-direction to the length of the template
xScale = templateLength/femurInertiaLength;
TFM2xScaling = eye(4); TFM2xScaling(1,1) = xScale;
femurXScaling = transformPoint3d(femurInertia, TFM2xScaling);
if debugVisu
    tPH = uipanel('Title','Template & length adj. subject','FontSize',14,'BorderWidth',2,...
        'BackgroundColor','w','Position',[BW+(BSX+BW)*1 0.01 BSX 0.99]);
    dAxH1(2) = axes('Parent', tPH, 'Visible','off', 'Color','w');
    view(dAxH1(2),dFig1View); axis(dAxH1(2), 'tight');
    % The scaled femur in the inertia CS
    visualizeMeshes(dAxH1(2), femurXScaling);
    templateProps.EdgeColor = 'none';
    templateProps.FaceColor = 'b';
    templateProps.FaceLighting = 'gouraud';
    % The template in the template CS
    patch(dAxH1(2), template, templateProps);
end

% Rough pre-registration of the subject to the template
femurPreReg.vertices = roughPreRegistration(template.vertices, femurXScaling.vertices);
femurPreReg.faces = femurInertia.faces;
if debugVisu
    tPH = uipanel('Title','Rough pre-registration','FontSize',14,'BorderWidth',2,...
        'BackgroundColor','w','Position',[BW+(BSX+BW)*2 0.01 BSX 0.99]);
    dAxH1(3) = axes('Parent', tPH, 'Visible','off', 'Color','w');
    view(dAxH1(3),dFig1View); axis(dAxH1(3), 'tight');
    % Subject after rough pre-registration
    visualizeMeshes(dAxH1(3), femurPreReg);
    patch(template, templateProps)
end

% Register femoral condyles of the subject to the template
femurCondReg.vertices = regFemoralCondyles(template.vertices, femurPreReg.vertices);
femurCondReg.faces = femurPreReg.faces;
if debugVisu
    tPH = uipanel('Title','Condyle registration','FontSize',14,'BorderWidth',2,...
        'BackgroundColor','w','Position',[BW+(BSX+BW)*3 0.01 BSX 0.99]);
    dAxH1(4) = axes('Parent', tPH, 'Visible','off', 'Color','w');
    view(dAxH1(4),dFig1View); axis(dAxH1(4), 'tight');
    % The subject after condyle registration
    visualizeMeshes(dAxH1(4), femurCondReg);
    patch(template, templateProps)
end

% Adapt femoral version, neck-shaft angle and neck length of the template
templatePreReg.vertices = adjustTemplateFemoralVersion(template.vertices, femurCondReg.vertices);
templatePreReg.faces = template.faces;
if debugVisu
    tPH = uipanel('Title','Neck & head registration','FontSize',14,'BorderWidth',2,...
        'BackgroundColor','w','Position',[BW+(BSX+BW)*4 0.01 BSX 0.99]);
    dAxH1(5) = axes('Parent', tPH, 'Visible','off', 'Color','w');
    view(dAxH1(5),dFig1View); axis(dAxH1(5), 'tight');
    % The adapted template
    visualizeMeshes(dAxH1(5), femurCondReg);
    patch(templatePreReg, templateProps)
end

% Non-rigid ICP registration
disp('____________ Morphing of the template mesh to the target _____________')
NRICP_ALPHA = [1e10 1e9 1e8 1e7 1e5 1e3 10 0.1 0.001]';
templateNICP = nonRigidICP(templatePreReg, femurCondReg, ...
    'alpha', NRICP_ALPHA,...
    'verbosity', double(verb));
if debugVisu
    tPH = uipanel('Title','Nonrigid ICP registration','FontSize',14,'BorderWidth',2,...
        'BackgroundColor','w','Position',[BW+(BSX+BW)*5 0.01 BSX 0.99]);
    dAxH1(6) = axes('Parent', tPH, 'Visible','off', 'Color','w');
    view(dAxH1(6),dFig1View); axis(dAxH1(6), 'tight');
    % The femur after NICP registration
    visualizeMeshes(dAxH1(6), femurCondReg);
    patch(templateNICP, templateProps)
end

% Mapping of landmarks and areas of the template to the subject
femurCondRegKDTree=createns(femurCondReg.vertices);
% Landmarks
load('template_landmarks.mat','landmarks')
if debugVisu
    drawPoint3d(dAxH1(6), templateNICP.vertices(cell2mat(struct2cell(landmarks)),:),...
        'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',3)
end
% Search the nearest neighbour of the landmarks on the NICP registered
% template to the vertices of the pre-registered femur. Return vertex
% indices.
landmarksIdx = knnsearch(femurCondRegKDTree, ...
    templateNICP.vertices(cell2mat(struct2cell(landmarks)),:));
LMNames = fieldnames(landmarks);
for lm=1:length(landmarksIdx)
    LMIdx.(LMNames{lm})=landmarksIdx(lm);
end
if debugVisu
    tPH = uipanel('Title','Mapping to the subject','FontSize',14,'BorderWidth',2,...
        'BackgroundColor','w','Position',[BW+(BSX+BW)*6 0.01 BSX 0.99]);
    dAxH1(7) = axes('Parent', tPH, 'Visible','off', 'Color','w');
    view(dAxH1(7),dFig1View); axis(dAxH1(7), 'tight');
    visualizeMeshes(dAxH1(7), femurInertia);
    drawPoint3d(dAxH1(7), femurInertia.vertices(landmarksIdx,:),...
        'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',3)
end

% Areas
load('template_areas.mat','areas')
% Delete not required area
areas(3,:)=[]; % Shaft
areas=[areas,cell(size(areas,1),1)];
for a=1:size(areas,1)
    tempFaces = template.faces(areas{a,2},:);
    [unqVertIds, ~, ~] = unique(tempFaces);
    tempVertices = templateNICP.vertices(unqVertIds,:);
    % Search the nearest neighbour of the landmarks on the NICP registered
    % template to the vertices of the pre-registered femur. Return vertex
    % indices.
    areas{a,3} = knnsearch(femurCondRegKDTree, tempVertices);
    if debugVisu
        drawPoint3d(dAxH1(6), tempVertices,...
            'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',3)
        drawPoint3d(dAxH1(7), femurInertia.vertices(areas{a,3},:),...
            'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',3)
    end
end

if debugVisu
%     % For publication
%     XLIM = [-270,240]; YLIM = [-50,45];
%     arrayfun(@(x) xlim(x,XLIM), dAxH1)
%     arrayfun(@(x) ylim(x,YLIM), dAxH1);
%     set(dFigH1, 'GraphicsSmoothing','off')
%     export_fig('Figure2', '-tif', '-r300')
end


%% Extract parameter
% Hip joint center
if isnan(FHC)
    % Fit sphere to the head of the femur, if FHC is not already available
    Head = femur.vertices(areas{ismember(areas(:,1),'Head'),3},:);
    headSphere = fitSphere(Head); 
    FHC = headSphere(1:3);
    % Radius = headSphere(4);
    if debugVisu
        patchProps.FaceColor = [223, 206, 161]/255;
        patchProps.FaceAlpha = 0.5;
        [~, dAxH2, dFigH2] = visualizeMeshes(femur, patchProps);
        dFigH2.Name=[subject ': initial FHC, neck axis, ... (Debug Figure)'];
        dFigH2.NumberTitle='Off';
        mouseControl3d(dAxH2)
        % drawSphere(dAxH2, [FHC, Radius])
        drawPoint3d(dAxH2, FHC,'MarkerFaceColor','k','MarkerEdgeColor','k');
        fontSize = 16;
        drawLabels3d(dAxH2, FHC,'FHC', 'FontSize',fontSize, 'VerticalAlignment','Top')
    end
end

if isnan(NeckAxis)
    % Neck axis
    Neck = femur.vertices(areas{ismember(areas(:,1),'Neck'),3},:);
    % Fit ellipse to the neck
    [NeckEllipse, NeckEllipseTFM] = fitEllipse3d(Neck);
    % Neck axis is defined by center and normal of the ellipse
    NeckAxis = [NeckEllipse(1:3), normalizeVector3d(transformVector3d([0 0 1], NeckEllipseTFM))];
end
NeckPlane=createPlane(NeckAxis(1:3), NeckAxis(4:6));
% Neck axis should point in lateral direction
if ~isBelowPlane(FHC,NeckPlane)
    NeckAxis=reverseLine3d(NeckAxis);
end
if debugVisu
    % drawEllipse3d(NeckEllipse)
end
% Use vertex indices of the mesh to define the neck axis
LMIdx.NeckAxis = lineToVertexIndices(NeckAxis, femur);
if debugVisu
    % debugNeckAxis = createLine3d(...
    %     femur.vertices(LMIdx.NeckAxis(1),:),...
    %     femur.vertices(LMIdx.NeckAxis(2),:));
    % debugNeckAxis(4:6) = normalizeVector3d(debugNeckAxis(4:6));
    % drawVector3d(dAxH2, debugNeckAxis(1:3),debugNeckAxis(4:6)*100,'-.r');
end

% Shaft Axis
[ShaftAxis, LMIdx.ShaftAxis] = detectShaftAxis(femur, FHC, 'debug',debugVisu);
if debugVisu
    % debugShaftAxis = createLine3d(...
    %     femur.vertices(LMIdx.ShaftAxis(1),:),...
    %     femur.vertices(LMIdx.ShaftAxis(2),:));
    % debugShaftAxis(4:6) = normalizeVector3d(debugShaftAxis(4:6));
    % drawVector3d(dAxH2, debugShaftAxis(1:3),debugShaftAxis(4:6)*100,'-.g');
    drawLine3d(dAxH2, ShaftAxis, 'Color','r', 'LineWidth',2)
end


%% Refinement of the neck axis
disp('_____________________ Refinement of the neck axis ____________________')
try
    NeckAxis = femoralNeckAxis(femur, side, NeckAxis, ShaftAxis, ...
        'visu',debugVisu, 'verbose',verb, 'subject',subject);
catch
    % In case the morphing of the neck did not work properly.
    FHC2ShaftAxis = projPointOnLine3d(FHC, ShaftAxis);
    NeckAxis = createLine3d(midPoint3d(FHC2ShaftAxis, FHC), FHC2ShaftAxis);
    OrthogonalAxis = [midPoint3d(FHC2ShaftAxis, FHC), crossProduct3d(NeckAxis(4:6), ShaftAxis(4:6))];
    DEF_NECK_SHAFT_ANGLE = 135;
    NeckAxis = transformLine3d(NeckAxis, ...
        createRotation3dLineAngle(OrthogonalAxis, deg2rad(DEF_NECK_SHAFT_ANGLE-90)));
    NeckAxis = femoralNeckAxis(femur, side, NeckAxis, ShaftAxis, ...
        'visu',debugVisu, 'verbose',verb,'subject',subject, 'PlaneVariationRange',12);
end
LMIdx.NeckAxis = lineToVertexIndices(NeckAxis, femur);
if debugVisu
    % For publication
    % set(gcf, 'GraphicsSmoothing','off')
    % export_fig('Figure4', '-tif', '-r300')
    
    % debugNeckAxis = createLine3d(...
    %     femur.vertices(LMIdx.NeckAxis(1),:),...
    %     femur.vertices(LMIdx.NeckAxis(2),:));
    % debugNeckAxis(4:6) = normalizeVector3d(debugNeckAxis(4:6));
    % drawVector3d(dAxH2, debugNeckAxis(1:3),debugNeckAxis(4:6)*100, 'r');
    NeckEdge = [NeckAxis(1:3) + NeckAxis(4:6)*50; NeckAxis(1:3) - NeckAxis(4:6)*50];
    drawEdge3d(dAxH2, NeckEdge, 'Color','g', 'LineWidth',2)
end


%% Detection of the tabletop plane
LMIdx = detectTabletopPlane(femur, side, FHC, NeckAxis, LMIdx, 'visu', debugVisu);


%% Refinement of the epicondyles
disp('___________________ Refinement of the epicondyles ____________________')
% Axis through the epicondyles in the inertia system
CondyleAxisInertia = createLine3d(...
    femurInertia.vertices(LMIdx.MedialEpicondyle,:),...
    femurInertia.vertices(LMIdx.LateralEpicondyle,:));
% Shaft axis in the inertia system
ShaftAxisInertia=createLine3d(...
    femurInertia.vertices(LMIdx.ShaftAxis(2),:),...
    femurInertia.vertices(LMIdx.ShaftAxis(1),:));
ShaftAxisInertia(4:6)=normalizeVector3d(ShaftAxisInertia(4:6));

% Cut the distal femur in the inertia system
if CondyleAxisInertia(1)>0; cutDir=1; else; cutDir=-1; end
distalCutPlaneInertia = createPlane([cutDir*1/4*femurInertiaLength 0 0], -ShaftAxisInertia(4:6));
distalFemurInertia = cutMeshByPlane(femurInertia, distalCutPlaneInertia);

% Transform into the inital USP CS
uspPreTFM = anatomicalOrientationTFM('RAS','ASR') * ...
    Tabletop(femurInertia, 'R', FHC, LMIdx, 'visu', false);
% The inertia femur is always a right femur because left femurs are mirrored
[USP_TFM, PFEA, CEA] = USP(distalFemurInertia, 'R', ...
    rotation3dToEulerAngles(uspPreTFM(1:3,1:3), 'ZYX'), ...
    'visu',debugVisu, ...
    'verbose',verb, ...
    'subject',subject);

% Transform distal femur into the USP CS
distalFemurUSP = transformPoint3d(distalFemurInertia, USP_TFM);
% Get the mapped points of the epicodyles in USP CS
MEC_map_USP = transformPoint3d(femurInertia.vertices(LMIdx.MedialEpicondyle,:), USP_TFM);
LEC_map_USP = transformPoint3d(femurInertia.vertices(LMIdx.LateralEpicondyle,:), USP_TFM);
% Refinement of the epicondyles (beta)
[MEC_USP, LEC_USP] = epicondyleRefinement(distalFemurUSP, CEA, ...
    MEC_map_USP, LEC_map_USP, 'visu',debugVisu);

% Get the epicondyle indices for the full femur
EC_Inertia = transformPoint3d([MEC_USP; LEC_USP], inv(USP_TFM));
[~, EC_Idx] = pdist2(femurInertia.vertices, EC_Inertia, 'euclidean', 'Smallest',1);
if side == 'L'; flipud(EC_Idx); end
LMIdx.MedialEpicondyle = EC_Idx(1); LMIdx.LateralEpicondyle = EC_Idx(2);

if debugVisu
    MEC = femur.vertices(LMIdx.MedialEpicondyle,:);
    LEC = femur.vertices(LMIdx.LateralEpicondyle,:);
    drawPoint3d(dAxH2, [MEC; LEC],'MarkerFaceColor','k','MarkerEdgeColor','k');
    drawLabels3d(dAxH2, MEC,'MEC', 'FontSize',fontSize, 'VerticalAlignment','Top')
    drawLabels3d(dAxH2, LEC,'LEC', 'FontSize',fontSize, 'VerticalAlignment','Top')
end


%% Refinement of the Intercondylar Notch (ICN)
LMIdx.ICN_mapped = LMIdx.IntercondylarNotch;
extremePoints = distalFemoralExtremePoints(distalFemurUSP, 'R', PFEA, 'visu', debugVisu, 'debug',0);
extremePointsInertia = structfun(@(x) transformPoint3d(x, inv(USP_TFM)), extremePoints,'uni',0);
[~, LMIdx.IntercondylarNotch] = pdist2(femurInertia.vertices, ...
    extremePointsInertia.Intercondylar, 'euclidean','Smallest',1);
[~, LMIdx.MedialProximoposteriorCondyle] = pdist2(femurInertia.vertices, ...
    extremePointsInertia.Medial, 'euclidean','Smallest',1);
[~, LMIdx.LateralProximoposteriorCondyle] = pdist2(femurInertia.vertices, ...
    extremePointsInertia.Lateral, 'euclidean','Smallest',1);

if debugVisu
    ICN = femur.vertices(LMIdx.IntercondylarNotch,:);
    drawPoint3d(dAxH2, ICN,'MarkerFaceColor','k','MarkerEdgeColor','k');
    drawLabels3d(dAxH2, ICN,'ICN', 'FontSize',fontSize, 'VerticalAlignment','Top')
    % MPPC = femur.vertices(LMIdx.MedialProximoposteriorCondyle,:);
    % drawPoint3d(dAxH2, MPPC,'MarkerFaceColor','k','MarkerEdgeColor','k');
    % drawLabels3d(dAxH2, MPPC,'MPPC')
    % LPPC = femur.vertices(LMIdx.LateralProximoposteriorCondyle,:);
    % drawPoint3d(dAxH2, LPPC,'MarkerFaceColor','k','MarkerEdgeColor','k');
    % drawLabels3d(dAxH2, LPPC,'LPPC')
    MPC = femur.vertices(LMIdx.MedialPosteriorCondyle,:);
    drawPoint3d(dAxH2, MPC,'MarkerFaceColor','k','MarkerEdgeColor','k');
    drawLabels3d(dAxH2, MPC,'MPC', 'FontSize',fontSize, 'VerticalAlignment','Top')
    LPC = femur.vertices(LMIdx.LateralPosteriorCondyle,:);
    drawPoint3d(dAxH2, LPC,'MarkerFaceColor','k','MarkerEdgeColor','k');
    drawLabels3d(dAxH2, LPC,'LPC', 'FontSize',fontSize, 'VerticalAlignment','Top')
    PTC = femur.vertices(LMIdx.PosteriorTrochantericCrest,:);
    drawPoint3d(dAxH2, PTC,'MarkerFaceColor','k','MarkerEdgeColor','k');
    drawLabels3d(dAxH2, PTC,'PTC', 'FontSize',fontSize, 'VerticalAlignment','Top')
    ttpProps.FaceColor='blue'; ttpProps.FaceAlpha=0.5; ttpProps.EdgeColor='none';
    ttpPatch.vertices=[MPC; LPC; PTC];
    ttpPatch.faces=[1 2 3];
    patch(dAxH2, ttpPatch, ttpProps)
end


%% PFEA and CEA
LM.PFEA = transformLine3d(PFEA, femurProps.inverseInertiaTFM*xReflection*inv(USP_TFM));  %#ok<MINV>
LM.CEA = transformLine3d(CEA, femurProps.inverseInertiaTFM*xReflection*inv(USP_TFM)); %#ok<MINV>
LMIdx.PFEA = lineToVertexIndices(LM.PFEA, femur);
LMIdx.CEA = lineToVertexIndices(LM.CEA, femur);


%% Construct the femoral CS
if visu
    CSIdx = ismember(validCSStrings,definition);
else
    CSIdx = false(4,1);
end
TFM.Wu2002 = Wu2002(femur,side,FHC,LMIdx, 'visu',CSIdx(1));
TFM.Bergmann2016 = Bergmann2016(femur, side, FHC, LMIdx, 'visu',CSIdx(2));
TFM.Tabletop = Tabletop(femur, side, FHC, LMIdx, 'visu',CSIdx(3));
TFM.MediTEC = MediTEC(femur, side, FHC, LMIdx, 'visu',CSIdx(4));
TFM.USP = USP_TFM*xReflection*inertiaTFM; %#ok<MINV>

TFM2AFCS = TFM.(definition);


%% Save landmarks in cartesian coordinates in input femur CS
LM.FemoralHeadCenter = FHC;
LM.NeckAxis = NeckAxis;
LM.ShaftAxis = ShaftAxis;
LM.MedialEpicondyle = femur.vertices(LMIdx.MedialEpicondyle,:);
LM.LateralEpicondyle = femur.vertices(LMIdx.LateralEpicondyle,:);
LM.IntercondylarNotch = femur.vertices(LMIdx.IntercondylarNotch,:);
LM.MedialPosteriorCondyle = femur.vertices(LMIdx.MedialPosteriorCondyle,:);
LM.LateralPosteriorCondyle = femur.vertices(LMIdx.LateralPosteriorCondyle,:);
LM.GreaterTrochanter = femur.vertices(LMIdx.GreaterTrochanter,:);
LM.LesserTrochanter = femur.vertices(LMIdx.LesserTrochanter,:);
LM.MedialProximoposteriorCondyle = femur.vertices(LMIdx.MedialProximoposteriorCondyle,:);
LM.LateralProximoposteriorCondyle = femur.vertices(LMIdx.LateralProximoposteriorCondyle,:);
LM.PosteriorTrochantericCrest = femur.vertices(LMIdx.PosteriorTrochantericCrest,:);

end