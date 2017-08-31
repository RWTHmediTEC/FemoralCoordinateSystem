function [fwTFM2AFCS, LMIdx, HJC] = automaticFemoralCS(femur, side, varargin)
%AUTOMATICFEMORALCS Calculate an femoral coordinate system (AFCS)
%
% REQUIRED INPUT:
%   femur: A "clean" mesh as struct with the fields vertices and faces
%   side: 'L' or 'R': left or right femur
%   
% OPTIONAL INPUT:
%   'HJC': 1x3 vector: Coordinates of the hip joint center in the 
%       coordinate system (CS) of the input femur mesh
%   'definition': The definition to construct the femoral CS
%       'Wu2002' - 2002 - Wu et al. - ISB recommendation on definitions
%                     of joint coordinate systems of various 
%                     joints for the reporting of human joint motion -
%                     part I: ankle, hip, and spine (default)
%       'WuBergmannComb' - A combination of Wu2002 and Bergmann2016
%   'iterations': number of iterations of the non-rigid ICP
%   'visualization': true (default) or false
% 
% OUTPUT:
%   fwTFM2AFCS: Forward transformation of the femur into the femoral CS
%   LMIdx: Landmark indices into the vertices of the femur
%
% TODO:
%   - Add other definitions: Bergmann2016, Murphy1987, etc.
%
% AUTHOR: Maximilian C. M. Fischer
% 	mediTEC - Chair of Medical Engineering, RWTH Aachen University
% VERSION: 1.0
% DATE: 2017-07-07

addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'src']));
addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'res']));


%% Parse inputs
p = inputParser;
addRequired(p,'femur',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addRequired(p,'side',@(x) any(validatestring(x,{'R','L'})));
isPoint = @(x) validateattributes(x,{'numeric'},...
    {'nonempty','nonnan','real','finite','size',[1,3]});
addParameter(p,'HJC',nan, isPoint);
addParameter(p,'definition','Wu2002',@(x) any(validatestring(x,{'Wu2002','WuBergmannComb'})));
addParameter(p,'iterations',10,...
    @(x)validateattributes(x,{'numeric'},{'scalar', '>',1, '<', 100}));
addParameter(p,'visualization',true,@islogical);

parse(p,femur,side,varargin{:});
femur = p.Results.femur;
side = p.Results.side;
HJC = p.Results.HJC;
definition = p.Results.definition;
nICPiter = p.Results.iterations;
visu = p.Results.visualization;

visuDebug = false;

%% Algorithm
% Get inertia transformation
femurProps = inertiaInfo(femur);

% Transform the vertices into the temporary inertia coordinate system
femurInertia.vertices = ...
    transformPoint3d(femur.vertices, inv(femurProps.inverseInertiaTFM));
femurInertia.faces = femur.faces;

% For left femurs reflect along the x-axis after inertia transform
if strcmp(side, 'L')
    xReflection=eye(4);
    xReflection(1,1)=-1;
    femurInertia.vertices=transformPoint3d(femurInertia.vertices, xReflection);
end

if visuDebug
    % New figure
    monitorsPosition = get(0,'MonitorPositions');
    FigHandle = figure('Units','pixels','renderer','opengl', 'Color', 'w','ToolBar','figure',...
        'WindowScrollWheelFcn',@zoomWithWheel,'WindowButtonDownFcn',@rotateWithLeftMouse);
    if     size(monitorsPosition,1) == 1
        set(FigHandle,'OuterPosition',monitorsPosition(1,:));
    elseif size(monitorsPosition,1) == 2
        set(FigHandle,'OuterPosition',monitorsPosition(2,:));
    end
    hold on
    % title({'The femur in the automatic femoral coordinate system (AFCS)';...
    title({'Left mouse - Rotate | Mouse wheel - Zoom'})
    cameratoolbar('SetCoordSys','none')
    axis equal; axis on; xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]');
    lightHandle(1) = light; light('Position', -1*(get(lightHandle(1),'Position')));
    view(90,0)
    hold on
    
    % Patch properties
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = 'r';
    patchProps.FaceAlpha = 0.75;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    
    % The femur in the inertia CS
    patch(femurInertia, patchProps)
end

% Load target mesh
load('target.mat','target')
if visuDebug
    % The target in the target CS
    patchProps.FaceColor = 'y';
    patch(target, patchProps)
end

% Length of the target femur
targetLength=max(target.vertices(:,1))-min(target.vertices(:,1));
% Length of the input femur
femurInertiaLength=max(femurInertia.vertices(:,1))-min(femurInertia.vertices(:,1));

% Scale input femur in x-direction
xScale=targetLength/femurInertiaLength;
TFM2xScaling=eye(4); TFM2xScaling(1,1)=xScale;
femurxScaling.vertices=transformPoint3d(femurInertia.vertices, TFM2xScaling);
femurxScaling.faces=femurInertia.faces;
if visuDebug
    % The scaled femur in the inertia CS
    patchProps.FaceColor = 'g';
    patch(femurxScaling, patchProps)
end

% Rough pre registration
femurPreReg.vertices = roughPreRegistration(target.vertices, femurxScaling.vertices);
femurPreReg.faces = femurInertia.faces;
if visuDebug
    % The femur in the pre-registration CS
    patchProps.FaceColor = 'b';
    patch(femurPreReg, patchProps)
end

% non-rigid ICP registration
[femurNICPReg.vertices,~,~]=....
    nonrigidICP(target, femurPreReg,...
    'preall',true,'vis',false,'iter',nICPiter,'verb',true);
femurNICPReg.faces=femurPreReg.faces;
if visuDebug
    % The femur after NICP registration
    patchProps.FaceColor = 'c';
    patch(femurNICPReg, patchProps)
end

% Mapping of landmarks and areas of the target to the source
femurNICPRegKDTree=createns(femurNICPReg.vertices);
% Landmarks
load('target.mat','landmarks')
landmarksIdx = knnsearch(femurNICPRegKDTree, cell2mat(landmarks(:,2)));
for lm=1:size(landmarks,1)
    LMIdx.(landmarks{lm,1})=landmarksIdx(lm);
end
if visuDebug
    drawPoint3d(femurInertia.vertices(landmarksIdx,:),...
        'MarkerFaceColor','k','MarkerEdgeColor','k')
end

% Areas
load('target.mat','areas')
areas=[areas,cell(size(areas,1),1)];
for a=1:size(areas,1)
    tempFaces = target.faces(areas{a,2},:);
    [unqVertIds, ~, ~] = unique(tempFaces);
    tempVertices = target.vertices(unqVertIds,:);
    areas{a,3} = knnsearch(femurNICPRegKDTree, tempVertices);
    if visuDebug
        drawPoint3d(femurInertia.vertices(areas{a,3},:),...
            'MarkerFaceColor','y','MarkerEdgeColor','y')
    end
end

% Construct the femoral CS
if isnan(HJC)
    HJC = femur.vertices(areas{ismember(areas(:,1),'Head'),3},:);
    % fit sphere to the head of the femur, if HJC is not a single point
    [HJC, ~] = spherefit(HJC); HJC=HJC';
end
switch definition
    case 'Wu2002'
        fwTFM2AFCS = Wu2002(femur,side,HJC,...
            femur.vertices(LMIdx.MedialEpicondyle,:),...
            femur.vertices(LMIdx.LateralEpicondyle,:),...
            'visu', visu);
    case 'WuBergmannComb'
        [fwTFM2AFCS, LMIdx.MedialPosteriorCondyle, LMIdx.LateralPosteriorCondyle]...
            = WuBergmannComb(femur,side,HJC,...
            femur.vertices(LMIdx.MedialPosteriorCondyle,:),...
            femur.vertices(LMIdx.LateralPosteriorCondyle,:),...
            femur.vertices(LMIdx.IntercondylarNotch,:),...
            'visu', visu);
end

end

function props = inertiaInfo(Mesh)

% Get Volume (V), Center of Mass (CoM), Inertia Tensor (J) of the Bone
[props.V, props.CoM, props.J] = VolumeIntegrate(Mesh.vertices, Mesh.faces);

% Get Principal Axes (pAxes) & Principal Moments of Inertia (Jii)
[props.pAxes, props.Jii] = eig(props.J); % Sign of the Eigenvectors can change (In agreement with their general definition)

% Keep the determinant positive
if det(props.pAxes) < 0
    props.pAxes = -1*props.pAxes;
end

% Create a affine transformation to move the Bone into his own Inertia System
props.inverseInertiaTFM = [ [props.pAxes props.CoM]; [0 0 0 1] ];

end

function zoomWithWheel(~,evnt)
if evnt.VerticalScrollCount > 0
    CVA_old = get(gca,'CameraViewAngle');
    CVA_new = CVA_old + 1;
    draw
elseif evnt.VerticalScrollCount < 0
    CVA_old = get(gca,'CameraViewAngle');
    CVA_new = CVA_old - 1;
    draw
end
    function draw
        set(gca,'CameraViewAngle',CVA_new)
        drawnow
    end
end

function rotateWithLeftMouse(src,~)
if strcmp(get(src,'SelectionType'),'normal')
    cameratoolbar('SetMode','orbit')
end
end