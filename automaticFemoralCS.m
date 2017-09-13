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
%           of joint coordinate systems of various joints for the reporting
%           of human joint motion - part I: ankle, hip, and spine (default)
%       'Bergmann2016' - 2016 - Bergmann et al. - Standardized Loads Acting
%           in Hip Implants
%       'WuBergmannComb' - A combination of Wu2002 and Bergmann2016
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
validStrings={'Wu2002','Bergmann2016','WuBergmannComb'};
addParameter(p,'definition','Wu2002',@(x) any(validatestring(x,validStrings)));
addParameter(p,'visualization',true,@islogical);

parse(p,femur,side,varargin{:});
femur = p.Results.femur;
side = p.Results.side;
HJC = p.Results.HJC;
definition = p.Results.definition;
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
    % Patch properties
    patchProps.EdgeColor = [0.5 0.5 0.5];
    patchProps.FaceColor = [.75 .75 .75];
    patchProps.FaceAlpha = 1;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    
    % The femur in the inertia CS
    visualizeMeshes(femurInertia, patchProps)
    hold on
end

% Load template mesh
load('template.mat','template')
if visuDebug
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
femurxScaling.vertices=transformPoint3d(femurInertia.vertices, TFM2xScaling);
femurxScaling.faces=femurInertia.faces;
if visuDebug
    % The scaled femur in the inertia CS
    patchProps.FaceColor = 'g';
    patch(femurxScaling, patchProps)
end

% Rough pre registration
femurPreReg.vertices = roughPreRegistration(template.vertices, femurxScaling.vertices);
femurPreReg.faces = femurInertia.faces;
if visuDebug
    % The femur in the pre-registration CS
    patchProps.FaceColor = 'b';
    patch(femurPreReg, patchProps)
end

% non-rigid ICP registration - mediTEC implementation
templateNICP = nonRigidICP(template, femurPreReg, 'alpha', [1e10 1e9 1e8 1e7 1e5 1e3 10 0.1 0.001]');
if visuDebug
    % The femur after NICP registration
    patchProps.FaceColor = 'c';
    patch(templateNICP, patchProps)
end

% Mapping of landmarks and areas of the template to the source
femurPreRegKDTree=createns(femurPreReg.vertices);
% Landmarks
load('template_landmarks.mat','landmarks')
% Search the nearest neighbour of the landmarks on the NICP registered 
% template to the vertices of the pre-registered femur. Return vertex
% indices.
landmarksIdx = knnsearch(femurPreRegKDTree, ...
    templateNICP.vertices(cell2mat(struct2cell(landmarks)),:));
LMNames = fieldnames(landmarks);
for lm=1:length(landmarksIdx)
    LMIdx.(LMNames{lm})=landmarksIdx(lm);
end
if visuDebug
    drawPoint3d(femurInertia.vertices(landmarksIdx,:),...
        'MarkerFaceColor','k','MarkerEdgeColor','k')
end

% Areas
load('template_areas.mat','areas')
areas=[areas,cell(size(areas,1),1)];
for a=1:size(areas,1)
    tempFaces = template.faces(areas{a,2},:);
    [unqVertIds, ~, ~] = unique(tempFaces);
    tempVertices = templateNICP.vertices(unqVertIds,:);
    % Search the nearest neighbour of the landmarks on the NICP registered
    % template to the vertices of the pre-registered femur. Return vertex
    % indices.
    areas{a,3} = knnsearch(femurPreRegKDTree, tempVertices);
    if visuDebug
        drawPoint3d(femurInertia.vertices(areas{a,3},:),...
            'MarkerFaceColor',[0.4 .08 .08],'MarkerEdgeColor',[0.4 .08 .08])
    end
end

if visuDebug
    legend({'Source Inertia','Template','Source Scaled',...
        'Source Pre-Registered','Template nICP'})
    mouseControl3d
end

% Construct the femoral CS
if isnan(HJC)
    % fit sphere to the head of the femur, if HJC is not already available
    Head = femur.vertices(areas{ismember(areas(:,1),'Head'),3},:);
    [HJC, Radius] = spherefit(Head); HJC=HJC';
    if visuDebug
        visualizeMeshes(femur)
        mouseControl3d
        hold on
        drawSphere([HJC, Radius])
    end
end
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
        % Neck axis
        Neck = femur.vertices(areas{ismember(areas(:,1),'Neck'),3},:);
        % Fit ellipse to the neck
        [NeckEllipse, NeckEllipseTFM] = fitEllipse3d(Neck);
        % Neck axis is defined by center and normal of the ellipse
        NeckAxis = [NeckEllipse(1:3), transformVector3d([0 0 1], NeckEllipseTFM)];
        % Use vertex indices of the mesh to define the neck axis
        NeckAxisPoints = intersectLineMesh3d(NeckAxis, femur.vertices, femur.faces);
        NeckAxisPoints = unique(NeckAxisPoints,'rows','stable');
        [~, NeckAxis_Idx] = pdist2(femur.vertices,NeckAxisPoints,'euclidean','Smallest',1);
        LMIdx.NeckAxis = [NeckAxis_Idx(1); NeckAxis_Idx(end)];
        if visuDebug
            drawLine3d(NeckAxis);
            drawEllipse3d(NeckEllipse)
        end
        % Shaft
        Shaft = femur.vertices(areas{ismember(areas(:,1),'Shaft'),3},:);
        [fwTFM2AFCS, LMIdx.MedialPosteriorCondyle, LMIdx.LateralPosteriorCondyle]...
            = Bergmann2016(femur, side, HJC, ...
            femur.vertices(LMIdx.MedialPosteriorCondyle,:),...
            femur.vertices(LMIdx.LateralPosteriorCondyle,:),...
            femur.vertices(LMIdx.IntercondylarNotch,:),...
            LMIdx.NeckAxis, Shaft, 'visu', visu);
end

end