function [fwTFM2AFCS, LMIdx, HJC] = automaticFemoralCS(femur, side, varargin)
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
%       'Tabletop' - A combination of Murphy1987 and Bergmann2016
%       'TabletopMediTEC' - Table top method taking into account the 
%           mechanical axis 
%   'visualization': true (default) or false
%   'subject': Char: Identification of the subject. Default is 'anonymous'.
% 
% OUTPUT:
%   fwTFM2AFCS: Forward transformation of the femur into the femoral CS
%   LMIdx: Landmark indices into the vertices of the femur
%
% TODO:
%   Add sanity check for MEC & LEC refinement for osteophytic distal femurs
%   Add intercondylar notch refinement (ICN)
%
% AUTHOR: Maximilian C. M. Fischer
% 	mediTEC - Chair of Medical Engineering, RWTH Aachen University
% VERSION: 1.0.2
% DATE: 2018-05-11

addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'src']));
addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'res']));

%% Parse inputs
p = inputParser;
addRequired(p,'femur',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addRequired(p,'side',@(x) any(validatestring(upper(x(1)),{'R','L'})));
isPoint3d = @(x) validateattributes(x,{'numeric'},...
    {'nonempty','nonnan','real','finite','size',[1,3]});
addParameter(p,'HJC',nan, isPoint3d);
isLine3d = @(x) validateattributes(x,{'numeric'},...
    {'nonempty','nonnan','real','finite','size',[1,6]});
addParameter(p,'NeckAxis',nan, isLine3d);
validStrings={'Wu2002','Bergmann2016','WuBergmannComb','Tabletop','TabletopMediTEC'};
addParameter(p,'definition','Wu2002',@(x) any(validatestring(x,validStrings)));
addParameter(p,'visualization',true,@islogical);
addParameter(p,'verbose',false,@islogical);
addParameter(p,'subject','anonymous',@(x) validateattributes(x,{'char'},{'nonempty'}));
addParameter(p,'debugVisualization',false,@islogical);

parse(p,femur,side,varargin{:});
femur = p.Results.femur;
side = upper(p.Results.side(1));
HJC = p.Results.HJC;
NeckAxis = p.Results.NeckAxis;
definition = p.Results.definition;
visu = p.Results.visualization;
verb = p.Results.verbose;
subject = p.Results.subject;
debugVisu = p.Results.debugVisualization;

%% Algorithm
% Get inertia transformation
femurProps = inertiaInfo(femur);

% Transform the vertices into the temporary inertia coordinate system
femurInertia = transformPoint3d(femur, inv(femurProps.inverseInertiaTFM));

% For left femurs reflect along the x-axis after inertia transform
if strcmp(side, 'L')
    xReflection=eye(4);
    xReflection(1,1)=-1;
    femurInertia.vertices=transformPoint3d(femurInertia.vertices, xReflection);
    % Orient the normals outwards after reflection
    femurInertia.faces=fliplr(femurInertia.faces);
end

if debugVisu
    % Patch properties
    patchProps.EdgeColor = [0.5 0.5 0.5];
    patchProps.FaceColor = [.75 .75 .75];
    patchProps.FaceAlpha = 1;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    
    % The femur in the inertia CS
    [~, ~, debugFigH1]=visualizeMeshes(femurInertia, patchProps);
    set(debugFigH1, 'Name', [subject ': nICP registration (Debug Figure)'], 'NumberTitle', 'Off')
    hold on
end

% Load template mesh
load('template.mat','template')
if debugVisu
    % The template in the template CS
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = 'y';
    patchProps.FaceAlpha = 0.75;
    patch(template, patchProps)
end

% Length of the template femur
templateLength=max(template.vertices(:,1))-min(template.vertices(:,1));
% Length of the input femur
femurInertiaLength=max(femurInertia.vertices(:,1))-min(femurInertia.vertices(:,1));

% Scale input femur in x-direction
xScale=templateLength/femurInertiaLength;
TFM2xScaling=eye(4); TFM2xScaling(1,1)=xScale;
femurxScaling=transformPoint3d(femurInertia, TFM2xScaling);
if debugVisu
    % The scaled femur in the inertia CS
    patchProps.FaceColor = 'g';
    patch(femurxScaling, patchProps)
end

% Rough pre-registration
femurPreReg.vertices = roughPreRegistration(template.vertices, femurxScaling.vertices);
femurPreReg.faces = femurInertia.faces;
if debugVisu
    % The femur in the pre-registration CS
    patchProps.FaceColor = 'b';
    patch(femurPreReg, patchProps)
end

% Register femoral condyles
femurCondReg.vertices = regFemoralCondyles(template.vertices, femurPreReg.vertices);
femurCondReg.faces = femurPreReg.faces;
if debugVisu
    % The femur in the ... CS
    patchProps.FaceColor = 'w';
    patch(femurCondReg, patchProps)
end


% Adapt femoral version of the template
templatePreReg.vertices = adjustTemplateFemoralVersion(template.vertices, femurCondReg.vertices);
templatePreReg.faces = template.faces;
if debugVisu
    % The template in the ... CS
    patchProps.FaceColor = 'm';
    patch(templatePreReg, patchProps)
end

% non-rigid ICP registration - mediTEC implementation
disp('____________ Morphing of the template mesh to the target _____________')
NRICP_ALPHA = [1e10 1e9 1e8 1e7 1e5 1e3 10 0.1 0.001]';
templateNICP = nonRigidICP(templatePreReg, femurCondReg, ...
    'alpha', NRICP_ALPHA,...
    'verbosity', double(verb));
if debugVisu
    % The femur after NICP registration
    patchProps.FaceColor = 'c';
    patch(templateNICP, patchProps)
end

% Mapping of landmarks and areas of the template to the source
femurCondRegKDTree=createns(femurCondReg.vertices);
% Landmarks
load('template_landmarks.mat','landmarks')
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
    drawPoint3d(femurInertia.vertices(landmarksIdx,:),...
        'MarkerFaceColor','k','MarkerEdgeColor','k')
end

% Areas
load('template_areas.mat','areas')
areas=[areas,cell(size(areas,1),1)]; %#ok<NODEF>
for a=1:size(areas,1)
    tempFaces = template.faces(areas{a,2},:);
    [unqVertIds, ~, ~] = unique(tempFaces);
    tempVertices = templateNICP.vertices(unqVertIds,:);
    % Search the nearest neighbour of the landmarks on the NICP registered
    % template to the vertices of the pre-registered femur. Return vertex
    % indices.
    areas{a,3} = knnsearch(femurCondRegKDTree, tempVertices);
    if debugVisu
        drawPoint3d(femurInertia.vertices(areas{a,3},:),...
            'MarkerFaceColor',[0.4 .08 .08],'MarkerEdgeColor',[0.4 .08 .08])
    end
end

if debugVisu
    legend({'Source Inertia','Template','Source Scaled',...
        'Source Pre-Registered','Source Condyle-Registered',...
        'Template Version-Registered','Template nICP'})
    mouseControl3d
end

%% Extract parameter
% Hip joint center
if isnan(HJC)
    % fit sphere to the head of the femur, if HJC is not already available
    Head = femur.vertices(areas{ismember(areas(:,1),'Head'),3},:);
    [HJC, Radius] = spherefit(Head); HJC=HJC';
    if debugVisu
        [~, debugAxH2, debugFigH2] = visualizeMeshes(femur);
        set(debugFigH2, 'Name', [subject ': initial HJC, neck axis, ... (Debug Figure)'],...
            'NumberTitle','Off')
        mouseControl3d
        hold on
        drawSphere([HJC, Radius])
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
    NeckAxis2=createLine3d(...
        femur.vertices(LMIdx.NeckAxis(1),:),...
        femur.vertices(LMIdx.NeckAxis(2),:));
    NeckAxis2(4:6)=normalizeVector3d(NeckAxis2(4:6));
    drawVector3d(NeckAxis2(1:3),NeckAxis2(4:6)*100,'r');
end

% Shaft Axis
Shaft = femur.vertices(areas{ismember(areas(:,1),'Shaft'),3},:);
% Fit ellipsoid to the shaft
ShaftEllipsoid = inertiaEllipsoid(Shaft);
% Construct the main shaft axis from the shaft ellipsoid
ShaftVector = transformVector3d([1 0 0], eulerAnglesToRotation3d(ShaftEllipsoid(7:9)));
ShaftAxis = [ShaftEllipsoid(1:3) ShaftVector];
ShaftPlane = createPlane(ShaftAxis(1:3), ShaftAxis(4:6));
if ~isBelowPlane(HJC,ShaftPlane)
    ShaftAxis=reverseLine3d(ShaftAxis);
end
% Use vertex indices of the mesh to define the shaft axis
LMIdx.ShaftAxis = lineToVertexIndices(ShaftAxis, femur);
if debugVisu
    ShaftAxis2=createLine3d(...
        femur.vertices(LMIdx.ShaftAxis(1),:),...
        femur.vertices(LMIdx.ShaftAxis(2),:));
    ShaftAxis2(4:6)=normalizeVector3d(ShaftAxis2(4:6));
    drawVector3d(ShaftAxis2(1:3),ShaftAxis2(4:6)*100,'g');
end

% Neck Orthogonal
NeckOrthogonal(1:3) = NeckAxis(1:3);
NeckOrthogonal(4:6) = crossProduct3d(NeckAxis(4:6), ShaftAxis(4:6));
% Use vertex indices of the mesh to define the neck orthogonal
if strcmp(side, 'L'); NeckOrthogonal(4:6)=-NeckOrthogonal(4:6); end
LMIdx.NeckOrthogonal = lineToVertexIndices(NeckOrthogonal, femur);
if debugVisu
    NeckOrthogonal2=createLine3d(...
        femur.vertices(LMIdx.NeckOrthogonal(1),:),...
        femur.vertices(LMIdx.NeckOrthogonal(2),:));
    NeckOrthogonal2(4:6)=normalizeVector3d(NeckOrthogonal2(4:6));
    drawVector3d(NeckOrthogonal2(1:3),NeckOrthogonal2(4:6)*100,'b');
end

%% Refinement of the neck axis
disp('_______________ Detection of the anatomical neck axis ________________')
NeckAxis = ANA(femur.vertices, femur.faces, side, ...
    LMIdx.NeckAxis, LMIdx.ShaftAxis, LMIdx.NeckOrthogonal,...
    'visu', debugVisu,'verbose',verb,'subject', subject);
LMIdx.NeckAxis = lineToVertexIndices(NeckAxis, femur);

%% Refinement of the epicondyles
disp('________________ Refinement of the epicondyles (beta) ________________')
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
distalCutPlaneInertia=createPlane([cutDir*1/4*femurInertiaLength 0 0], -ShaftAxisInertia(4:6));
distalFemurInertia = cutMeshByPlane(femurInertia, distalCutPlaneInertia);

% Transform into the USP CS
Y = normalizeVector3d(ShaftAxisInertia(4:6));
X = normalizeVector3d(crossProduct3d(ShaftAxisInertia(4:6), CondyleAxisInertia(4:6)));
Z = normalizeVector3d(crossProduct3d(X, Y));
uspInitialTFM = inv([[inv([X; Y; Z]), HJC']; [0 0 0 1]]);
uspInitialRot = rotation3dToEulerAngles(uspInitialTFM);

USP_TFM=USP(distalFemurInertia.vertices, distalFemurInertia.faces, side, uspInitialRot,...
    'visu', debugVisu,'verbose',verb, 'subject', subject);

% Get indices of the epicodyles
distalFemurUSP=transformPoint3d(distalFemurInertia, USP_TFM);
[~, minIdx]=min(distalFemurUSP.vertices(:,3));
[~, maxIdx]=max(distalFemurUSP.vertices(:,3));
ecInertia=transformPoint3d(distalFemurUSP.vertices([minIdx,maxIdx],:),inv(USP_TFM));
[~, ecIdx] = pdist2(femurInertia.vertices,ecInertia,'euclidean','Smallest',1);

if strcmp(side, 'L'); flipud(ecIdx); end
LMIdx.MedialEpicondyle=ecIdx(1); LMIdx.LateralEpicondyle=ecIdx(2);

if debugVisu
    drawPoint3d(debugAxH2, ...
        femur.vertices([LMIdx.MedialEpicondyle,LMIdx.LateralEpicondyle],:),...
        'MarkerFaceColor','r','MarkerEdgeColor','r');
end

%% Construct the femoral CS
switch definition
    case 'Wu2002'
        fwTFM2AFCS = Wu2002(femur,side,HJC,...
            femur.vertices(LMIdx.MedialEpicondyle,:),...
            femur.vertices(LMIdx.LateralEpicondyle,:),...
            'visu', visu);
    case 'WuBergmannComb'
        [fwTFM2AFCS, ...
            LMIdx.MedialPosteriorCondyle, ...
            LMIdx.LateralPosteriorCondyle]...
            = WuBergmannComb(femur,side,HJC,...
            femur.vertices(LMIdx.MedialPosteriorCondyle,:),...
            femur.vertices(LMIdx.LateralPosteriorCondyle,:),...
            femur.vertices(LMIdx.IntercondylarNotch,:),...
            'visu', visu);
    case 'Bergmann2016'
        [fwTFM2AFCS, LMIdx] = Bergmann2016(femur, side, HJC, LMIdx, 'visu', visu);
    case 'Tabletop'
        [fwTFM2AFCS, ...
            LMIdx.MedialPosteriorCondyle, ...
            LMIdx.LateralPosteriorCondyle, ...
            LMIdx.PosteriorTrochantericCrest]...
            = Tabletop(femur, side, HJC, ...
            femur.vertices(LMIdx.MedialPosteriorCondyle,:),...
            femur.vertices(LMIdx.LateralPosteriorCondyle,:),...
            femur.vertices(LMIdx.IntercondylarNotch,:),...
            NeckAxis, ShaftAxis, 'visu', visu);
    case 'TabletopMediTEC'
        [fwTFM2AFCS, ...
            LMIdx.MedialPosteriorCondyle, ...
            LMIdx.LateralPosteriorCondyle, ...
            LMIdx.PosteriorTrochantericCrest]...
            = TabletopMediTEC(femur, side, HJC, ...
            femur.vertices(LMIdx.MedialPosteriorCondyle,:),...
            femur.vertices(LMIdx.LateralPosteriorCondyle,:),...
            femur.vertices(LMIdx.IntercondylarNotch,:),...
            NeckAxis, 'visu', visu);
end

end