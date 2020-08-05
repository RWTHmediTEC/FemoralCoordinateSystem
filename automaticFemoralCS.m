function [fwTFM2AFCS, LMIdx, HJC, LM] = automaticFemoralCS(femur, side, varargin)
%AUTOMATICFEMORALCS Calculate an femoral coordinate system (AFCS)
%
% REQUIRED INPUT:
%   femur: A "clean" mesh as struct with the fields vertices and faces
%       "Clean" means a closed outer surface without holes or tunnels,
%       normals are oriented outwards, no duplicate vertices, ...
%   side: 'L' or 'R': left or right femur
%
% OPTIONAL INPUT:

%   'HJC': 1x3 vector: Coordinates of the hip joint center in the
%       coordinate system (CS) of the input femur mesh
%   'definition': The definition to construct the femoral CS
%       'Wu2002' - 2002 - Wu et al. - ISB recommendation on definitions
%           of joint coordinate systems of various joints for the reporting
%           of human joint motion - part I: ankle, hip, and spine (default)
%       'Bergmann2016' - 2016 - Bergmann et al. - Standardized Loads Acting
%           in Hip Implants
%       'WuBergmannComb' - A combination of Wu2002 and Bergmann2016
%       'Tabletop' - Defined by the table top plane
%       'MediTEC' - Defined by the mechanical axis and table top plane
%   'visualization': true (default) or false
%   'subject': Char: Identification of the subject. Default is 'anonymous'.
%
% OUTPUT:
%   fwTFM2AFCS: Forward transformation of the femur into the femoral CS
%   LMIdx: Landmark indices into the vertices of the femur
%
% TODO:
%   - Add option to select of anatomical orientation of the femur CS
%
% AUTHOR: Maximilian C. M. Fischer
% 	mediTEC - Chair of Medical Engineering, RWTH Aachen University
% VERSION: 1.0.3
% DATE: 2018-08-24

addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'src']));
addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'res']));

%% Parse inputs
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addRequired(p,'femur',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addRequired(p,'side',@(x) any(validatestring(upper(x(1)),{'R','L'})));
isPoint3d = @(x) validateattributes(x,{'numeric'},...
    {'nonempty','nonnan','real','finite','size',[1,3]});
addParameter(p,'HJC',nan, isPoint3d);
isLine3d = @(x) validateattributes(x,{'numeric'},...
    {'nonempty','nonnan','real','finite','size',[1,6]});
addParameter(p,'NeckAxis',nan, isLine3d);
validStrings={'Wu2002','Bergmann2016','WuBergmannComb','Tabletop','MediTEC'};
addParameter(p,'definition','Wu2002',@(x) any(validatestring(x,validStrings)));
addParameter(p,'visualization',true,logParValidFunc);
addParameter(p,'verbose',false,logParValidFunc);
addParameter(p,'subject','anonymous',@(x) validateattributes(x,{'char'},{'nonempty'}));
addParameter(p,'debugVisualization',false,logParValidFunc);

parse(p,femur,side,varargin{:});
femur = p.Results.femur;
side = upper(p.Results.side(1));
HJC = p.Results.HJC;
NeckAxis = p.Results.NeckAxis;
definition = lower(p.Results.definition);
visu = logical(p.Results.visualization);
verb = logical(p.Results.verbose);
subject = p.Results.subject;
debugVisu = logical(p.Results.debugVisualization);

%% Algorithm
% Get inertia transformation for the input femur (subject)
femurProps = inertiaInfo(femur);

% Transform the vertices into the temporary inertia coordinate system
femurInertia = transformPoint3d(femur, inv(femurProps.inverseInertiaTFM));

% For left femurs reflect along the x-axis after inertia transform
if strcmp(side, 'L')
    xReflection = eye(4);
    xReflection(1,1) = -1;
    femurInertia.vertices = transformPoint3d(femurInertia.vertices, xReflection);
    % Orient the normals outwards after reflection
    femurInertia.faces = fliplr(femurInertia.faces);
end

if debugVisu
    NOP = 7;
    BW = 0.005;
    BSX = (1-(NOP+1)*BW)/NOP;
    set(0,'defaultAxesFontSize',14)
    debugFigH1 = figure('Units','pixels','renderer','opengl', 'Color', 'w');
    MonitorsPos = get(0,'MonitorPositions');
    if     size(MonitorsPos,1) == 1
        set(debugFigH1,'OuterPosition', MonitorsPos(1,:));
    elseif size(MonitorsPos,1) == 2
        set(debugFigH1,'OuterPosition', MonitorsPos(2,:));
    end
    debugFigH1.Name = [subject ': nICP registration (Debug Figure)'];
    debugFigH1.NumberTitle = 'off';
    debugFigH1.WindowState = 'maximized';
    
    tPH = uipanel('Title','Subject''s principal axes transf. ','FontSize',14,'BorderWidth',2,...
        'BackgroundColor','w','Position',[BW+(BSX+BW)*0 0.01 BSX 0.99]);
    tAH = axes('Parent', tPH, 'Visible','off', 'Color','w');
    
    % The subject in the inertia CS
    visualizeMeshes(tAH, femurInertia);
    dFig1View = [90,-90];
    view(tAH,dFig1View); axis(tAH, 'tight');
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
    tAH = axes('Parent', tPH, 'Visible','off', 'Color','w');
    view(tAH,dFig1View); axis(tAH, 'tight');
    % The scaled femur in the inertia CS
    visualizeMeshes(tAH, femurXScaling);
    templateProps.EdgeColor = 'none';
    templateProps.FaceColor = 'b';
    templateProps.FaceLighting = 'gouraud';
    % The template in the template CS
    patch(tAH, template, templateProps);
end

% Rough pre-registration of the subject to the template
femurPreReg.vertices = roughPreRegistration(template.vertices, femurXScaling.vertices);
femurPreReg.faces = femurInertia.faces;
if debugVisu
    tPH = uipanel('Title','Rough pre-registration','FontSize',14,'BorderWidth',2,...
        'BackgroundColor','w','Position',[BW+(BSX+BW)*2 0.01 BSX 0.99]);
    tAH = axes('Parent', tPH, 'Visible','off', 'Color','w');
    view(tAH,dFig1View); axis(tAH, 'tight');
    % Subject after rough pre-registration
    visualizeMeshes(tAH, femurPreReg);
    patch(template, templateProps)
end

% Register femoral condyles of the subject to the template
femurCondReg.vertices = regFemoralCondyles(template.vertices, femurPreReg.vertices);
femurCondReg.faces = femurPreReg.faces;
if debugVisu
    tPH = uipanel('Title','Condyle registration','FontSize',14,'BorderWidth',2,...
        'BackgroundColor','w','Position',[BW+(BSX+BW)*3 0.01 BSX 0.99]);
    tAH = axes('Parent', tPH, 'Visible','off', 'Color','w');
    view(tAH,dFig1View); axis(tAH, 'tight');
    % The subject after condyle registration
    visualizeMeshes(tAH, femurCondReg);
    patch(template, templateProps)
end

% Adapt femoral version, neck-shaft angle and neck length of the template
templatePreReg.vertices = adjustTemplateFemoralVersion(template.vertices, femurCondReg.vertices);
templatePreReg.faces = template.faces;
if debugVisu
    tPH = uipanel('Title','Neck & head registration','FontSize',14,'BorderWidth',2,...
        'BackgroundColor','w','Position',[BW+(BSX+BW)*4 0.01 BSX 0.99]);
    tAH = axes('Parent', tPH, 'Visible','off', 'Color','w');
    view(tAH,dFig1View); axis(tAH, 'tight');
    % The adapted template
    visualizeMeshes(tAH, femurCondReg);
    patch(templatePreReg, templateProps)
end

% Non-rigid ICP registration
disp('____________ Morphing of the template mesh to the target _____________')
NRICP_ALPHA = [1e10 1e9 1e8 1e7 1e5 1e3 10 0.1 0.001]';
templateNICP = nonRigidICP(templatePreReg, femurCondReg, ...
    'alpha', NRICP_ALPHA,...
    'verbosity', double(verb));
if debugVisu
    tPH = uipanel('Title','Non-rigid ICP registration','FontSize',14,'BorderWidth',2,...
        'BackgroundColor','w','Position',[BW+(BSX+BW)*5 0.01 BSX 0.99]);
    tAH5 = axes('Parent', tPH, 'Visible','off', 'Color','w');
    view(tAH5,dFig1View); axis(tAH5, 'tight');
    % The femur after NICP registration
    visualizeMeshes(tAH5, femurCondReg);
    patch(templateNICP, templateProps)
end

% Mapping of landmarks and areas of the template to the subject
femurCondRegKDTree=createns(femurCondReg.vertices);
% Landmarks
load('template_landmarks.mat','landmarks')
if debugVisu
    drawPoint3d(tAH5, templateNICP.vertices(cell2mat(struct2cell(landmarks)),:),...
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
    tAH6 = axes('Parent', tPH, 'Visible','off', 'Color','w');
    view(tAH6,dFig1View); axis(tAH6, 'tight');
    visualizeMeshes(tAH6, femurInertia);
    drawPoint3d(tAH6, femurInertia.vertices(landmarksIdx,:),...
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
        drawPoint3d(tAH5, tempVertices,...
            'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',3)
        drawPoint3d(tAH6, femurInertia.vertices(areas{a,3},:),...
            'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',3)
    end
end

if debugVisu
    % % For publication
    % set(gcf, 'GraphicsSmoothing','off')
    % export_fig('Figure2', '-tif', '-r300')
end

%% Extract parameter
% Hip joint center
if isnan(HJC)
    % Fit sphere to the head of the femur, if HJC is not already available
    Head = femur.vertices(areas{ismember(areas(:,1),'Head'),3},:);
    headSphere = fitSphere(Head); 
    HJC = headSphere(1:3);
    Radius = headSphere(4);
    if debugVisu
        patchProps.FaceColor = [223, 206, 161]/255;
        patchProps.FaceAlpha = 0.5;
        [~, debugAxH2, debugFigH2] = visualizeMeshes(femur, patchProps);
        debugFigH2.Name=[subject ': initial HJC, neck axis, ... (Debug Figure)'];
        debugFigH2.NumberTitle='Off';
        mouseControl3d(debugAxH2)
        drawSphere(debugAxH2, [HJC, Radius])
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
if ~isBelowPlane(HJC,NeckPlane)
    NeckAxis=reverseLine3d(NeckAxis);
end
if debugVisu
    drawEllipse3d(NeckEllipse)
end
% Use vertex indices of the mesh to define the neck axis
LMIdx.NeckAxis = lineToVertexIndices(NeckAxis, femur);
if debugVisu
    debugNeckAxis = createLine3d(...
        femur.vertices(LMIdx.NeckAxis(1),:),...
        femur.vertices(LMIdx.NeckAxis(2),:));
    debugNeckAxis(4:6) = normalizeVector3d(debugNeckAxis(4:6));
    drawVector3d(debugAxH2, debugNeckAxis(1:3),debugNeckAxis(4:6)*100,'-.r');
end

% Shaft Axis
[ShaftAxis, LMIdx.ShaftAxis] = detectShaftAxis(femur, HJC, 'debug',debugVisu);
if debugVisu
    debugShaftAxis = createLine3d(...
        femur.vertices(LMIdx.ShaftAxis(1),:),...
        femur.vertices(LMIdx.ShaftAxis(2),:));
    debugShaftAxis(4:6) = normalizeVector3d(debugShaftAxis(4:6));
    drawVector3d(debugAxH2, debugShaftAxis(1:3),debugShaftAxis(4:6)*100,'-.g');
end

% Neck Orthogonal
NeckOrthogonal(1:3) = NeckAxis(1:3);
NeckOrthogonal(4:6) = crossProduct3d(NeckAxis(4:6), ShaftAxis(4:6));
% Use vertex indices of the mesh to define the neck orthogonal
if strcmp(side, 'L'); NeckOrthogonal(4:6)=-NeckOrthogonal(4:6); end
LMIdx.NeckOrthogonal = lineToVertexIndices(NeckOrthogonal, femur);
if debugVisu
    debugNeckOrthogonal = createLine3d(...
        femur.vertices(LMIdx.NeckOrthogonal(1),:),...
        femur.vertices(LMIdx.NeckOrthogonal(2),:));
    debugNeckOrthogonal(4:6) = normalizeVector3d(debugNeckOrthogonal(4:6));
    drawVector3d(debugAxH2, debugNeckOrthogonal(1:3),debugNeckOrthogonal(4:6)*100,'-.b');
end


%% Refinement of the neck axis
disp('_______________ Refinement of the anatomical neck axis ________________')
NeckAxis = ANA(femur.vertices, femur.faces, side, ...
    LMIdx.NeckAxis, LMIdx.ShaftAxis, LMIdx.NeckOrthogonal,...
    'visu', debugVisu,'verbose',verb,'subject', subject);
LMIdx.NeckAxis = lineToVertexIndices(NeckAxis, femur);
if debugVisu
    % For publication
    % set(gcf, 'GraphicsSmoothing','off')
    % export_fig('Figure4', '-tif', '-r300')
    
    debugNeckAxis = createLine3d(...
        femur.vertices(LMIdx.NeckAxis(1),:),...
        femur.vertices(LMIdx.NeckAxis(2),:));
    debugNeckAxis(4:6) = normalizeVector3d(debugNeckAxis(4:6));
    drawVector3d(debugAxH2, debugNeckAxis(1:3),debugNeckAxis(4:6)*100,'r');
end
% Neck Orthogonal
NeckOrthogonal(1:3) = NeckAxis(1:3);
NeckOrthogonal(4:6) = crossProduct3d(NeckAxis(4:6), ShaftAxis(4:6));
% Use vertex indices of the mesh to define the neck orthogonal
if strcmp(side, 'L'); NeckOrthogonal(4:6)=-NeckOrthogonal(4:6); end
LMIdx.NeckOrthogonal = lineToVertexIndices(NeckOrthogonal, femur);


%% Detection of the tabletop plane
LMIdx = detectTabletopPlane(femur, side, HJC, LMIdx, 'visu', debugVisu);


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
    Tabletop(femurInertia, 'R', HJC, LMIdx, 'visu', false);
% The inertia femur is always a right femur because left femurs are mirrored
[USP_TFM, PFEA, CEA] = USP(distalFemurInertia.vertices, distalFemurInertia.faces, ...
    'R', ...
    rotation3dToEulerAngles(uspPreTFM(1:3,1:3), 'ZYX'), ...
    'visu',debugVisu, ...
    'verbose',verb, ...
    'subject',subject);
if debugVisu
    % % For publication
    % set(gcf, 'GraphicsSmoothing','off')
    % export_fig('Figure6', '-tif', '-r300')
end

% Transform distal femur into the USP CS
distalFemurUSP = transformPoint3d(distalFemurInertia, USP_TFM);
% Get the mapped points of the epicodyles in USP CS
MEC_map_USP = transformPoint3d(femurInertia.vertices(LMIdx.MedialEpicondyle,:), USP_TFM);
LEC_map_USP = transformPoint3d(femurInertia.vertices(LMIdx.LateralEpicondyle,:), USP_TFM);
% Refinement of the epicondyles (beta)
[MEC_USP, LEC_USP] = epicondyleRefinement(distalFemurUSP, CEA, ...
    MEC_map_USP, LEC_map_USP, 'visu',debugVisu);

% Get the epicondyle indices for the full femur
ecInertia = transformPoint3d([MEC_USP; LEC_USP],inv(USP_TFM));
[~, ecIdx] = pdist2(femurInertia.vertices,ecInertia,'euclidean','Smallest',1);

if side == 'L'; flipud(ecIdx); end
LMIdx.MedialEpicondyle=ecIdx(1); LMIdx.LateralEpicondyle=ecIdx(2);

if debugVisu
    MEC = femur.vertices(LMIdx.MedialEpicondyle,:);
    LEC = femur.vertices(LMIdx.LateralEpicondyle,:);
    drawPoint3d(debugAxH2, [MEC; LEC],'MarkerFaceColor','r','MarkerEdgeColor','r');
    text(debugAxH2, MEC(1),MEC(2),MEC(3),'MEC')
    text(debugAxH2, LEC(1),LEC(2),LEC(3),'LEC')
end


%% Refinement of the Intercondylar Notch (ICN)
extremePoints = distalFemurExtremePoints(distalFemurUSP, 'R', PFEA, 'visu', debugVisu, 'debug',0);
extremePointsInertia = structfun(@(x) transformPoint3d(x, inv(USP_TFM)), extremePoints,'uni',0);
[~, LMIdx.IntercondylarNotch] = pdist2(femurInertia.vertices, ...
    extremePointsInertia.Intercondylar, 'euclidean','Smallest',1);
[~, LMIdx.ProximoposteriorMedialCondyle] = pdist2(femurInertia.vertices, ...
    extremePointsInertia.Medial, 'euclidean','Smallest',1);
[~, LMIdx.ProximoposteriorLateralCondyle] = pdist2(femurInertia.vertices, ...
    extremePointsInertia.Lateral, 'euclidean','Smallest',1);

if debugVisu
    ICN = femur.vertices(LMIdx.IntercondylarNotch,:);
    drawPoint3d(debugAxH2, ICN,'MarkerFaceColor','r','MarkerEdgeColor','r');
    text(debugAxH2, ICN(1),ICN(2),ICN(3),'ICN')
    PPMC = femur.vertices(LMIdx.ProximoposteriorMedialCondyle,:);
    drawPoint3d(debugAxH2, PPMC,'MarkerFaceColor','r','MarkerEdgeColor','r');
    text(debugAxH2, PPMC(1),PPMC(2),PPMC(3),'PPMC')
    PPLC = femur.vertices(LMIdx.ProximoposteriorLateralCondyle,:);
    drawPoint3d(debugAxH2, PPLC,'MarkerFaceColor','r','MarkerEdgeColor','r');
    text(debugAxH2, PPLC(1),PPLC(2),PPLC(3),'PPLC')
end


%% Construct the femoral CS
switch definition
    case 'wu2002'
        fwTFM2AFCS = Wu2002(femur,side,HJC,LMIdx, 'visu',visu);
    case 'wubergmanncomb'
        fwTFM2AFCS = WuBergmannComb(femur,side,HJC,LMIdx, 'visu',visu);
    case 'bergmann2016'
        fwTFM2AFCS = Bergmann2016(femur, side, HJC, LMIdx, 'visu',visu);
    case 'tabletop'
        fwTFM2AFCS = Tabletop(femur, side, HJC, LMIdx, 'visu',visu);
    case 'meditec'
        fwTFM2AFCS = MediTEC(femur, side, HJC, LMIdx, 'visu',visu);
end

%% Save landmarks in cartesian coordinates in input femur CS
LM.FemoralHeadCenter = HJC;
LM.NeckAxis = NeckAxis;
LM.ShaftAxis = ShaftAxis;
LM.MedialEpicondyle = femur.vertices(LMIdx.MedialEpicondyle,:);
LM.LateralEpicondyle = femur.vertices(LMIdx.LateralEpicondyle,:);
LM.IntercondylarNotch = femur.vertices(LMIdx.IntercondylarNotch,:);
LM.MedialPosteriorCondyle = femur.vertices(LMIdx.MedialPosteriorCondyle,:);
LM.LateralPosteriorCondyle = femur.vertices(LMIdx.LateralPosteriorCondyle,:);
LM.GreaterTrochanter = femur.vertices(LMIdx.GreaterTrochanter,:);
LM.LesserTrochanter = femur.vertices(LMIdx.LesserTrochanter,:);
LM.ProximoposteriorMedialCondyle = femur.vertices(LMIdx.ProximoposteriorMedialCondyle,:);
LM.ProximoposteriorLateralCondyle = femur.vertices(LMIdx.ProximoposteriorLateralCondyle,:);
LM.PosteriorTrochantericCrest = femur.vertices(LMIdx.PosteriorTrochantericCrest,:);

end