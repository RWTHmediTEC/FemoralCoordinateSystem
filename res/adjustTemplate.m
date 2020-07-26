clearvars; close all; opengl hardware

addpath('..\..\..\..\General\Code\#public')
addpath(genpath('..\src'))

load('template_mesh.mat')

%% Switches
manualSelectionSwitch = 0;
inertiaSwitch = 0;
controlSwitch = 0;
weightSwitch = 0;

%% Manual landmark selection
if manualSelectionSwitch
    % \\sagnix\Programmierung\public\Matlab\ManualLandmarkSelection
    addpath('..\..\..\..\General\Code\ManualLandmarkSelection')
    
    if exist('template_landmarks.mat','file')
        load('template_landmarks.mat','landmarks')
        LM=[fieldnames(landmarks),...
            mat2cell(template_mesh.vertices(cell2mat(struct2cell(landmarks)),:),...
            ones(1,length(fieldnames(landmarks))),3)];
    else
        LM = {'GreaterTrochanter';'LesserTrochanter';...
            'MedialEpicondyle';'LateralEpicondyle';...
            'MedialPosteriorCondyle';'LateralPosteriorCondyle';...
            'IntercondylarNotch'};
    end
    LM = selectLandmarks(template_mesh, LM);
    template_mesh_KDTree=createns(template_mesh.vertices);
    Idx = knnsearch(template_mesh_KDTree, cell2mat(LM(:,2)));
    for lm=1:size(LM,1)
        landmarks.(LM{lm,1})=Idx(lm);
    end
    save('template_landmarks.mat','landmarks')
    
    %% Manual face selection
    addpath('..\..\..\..\General\Code\ManualFaceSelection')
    
    if exist('template_areas.mat','file')
        load('template_areas.mat','areas')
    else
        areas = {'Head';'Neck';'Shaft'};
    end
    areas = selectFaces(template_mesh, areas);
    
    save('template_areas.mat','areas')
end

%% Inertia Alignment
if inertiaSwitch
    % Get inertia transformation
    inertiaProps = inertiaInfo(template_mesh);
    inertiaTFM = inv(inertiaProps.inverseInertiaTFM);
    % Transform the vertices into the temporary inertia coordinate system
    meshInertia.vertices = ...
        transformPoint3d(template_mesh.vertices, inertiaTFM);
    meshInertia.faces = template_mesh.faces;
    
    MechanicalAxis = [0 0 0 1 0 0];
    
    MEC=meshInertia.vertices(landmarks.MedialEpicondyle,:);
    LEC=meshInertia.vertices(landmarks.LateralEpicondyle,:);
    EpicondyleAxis = createLine3d(MEC, LEC);
    
    Z = normalizeVector3d(crossProduct3d(EpicondyleAxis(4:6),MechanicalAxis(4:6)));
    
    X = normalizeVector3d(-MechanicalAxis(4:6));
    
    Y = normalizeVector3d(crossProduct3d(Z, X));
    
    TFM = inv([[inv([X; Y; Z]), [0; 0; 0]]; [0 0 0 1]]);
    
    template.vertices = transformPoint3d(meshInertia.vertices, TFM);
    template.faces=meshInertia.faces;
    
    visualizeMeshes(template)
    mouseControl3d
    
    save('template.mat','template')
end

%% Construct controls
if controlSwitch
    addpath(genpath('..'))
    addpath('D:\Biomechanics\General\Code\#public')
    addpath(genpath('D:\Biomechanics\General\Code\#external\matGeom\matGeom'))
    load('template.mat','template')
    
    % Get HJC, neck axis & shaft axis
    [TFM2AFCS, LMIdx, HJC] = automaticFemoralCS(template, 'R',...
        'definition','Bergmann2016','vis', true, 'verbose', true);
    save('template_controls.mat', 'TFM2AFCS', 'LMIdx', 'HJC')
    
    load('template_controls.mat')
    
    neckAxis=createLine3d(template.vertices(LMIdx.NeckAxis(1),:), template.vertices(LMIdx.NeckAxis(2),:));
    shaftAxis=createLine3d(template.vertices(LMIdx.ShaftAxis(1),:), template.vertices(LMIdx.ShaftAxis(2),:));
    [~, P1, ] = distanceLines3d(neckAxis, shaftAxis);
    
    C(1,:)=HJC;
    C(2,:)=P1;
    C(3,:)=template.vertices(LMIdx.IntercondylarNotch,:);
    C(4,:)=template.vertices(LMIdx.MedialPosteriorCondyle,:);
    C(5,:)=template.vertices(LMIdx.LateralPosteriorCondyle,:);
    
    BE=[1,2; 2,3; 4,5];
    
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [223, 206, 161]/255;
    patchProps.FaceAlpha = 0.5;
    patchProps.FaceLighting = 'gouraud';
    visualizeMeshes(template, patchProps)
    
    for b=1:size(BE,1)
        drawEdge3d(C(BE(b,1),:),C(BE(b,2),:),'Color','k');
    end
    
    pointProps.Marker='o';
    pointProps.MarkerFaceColor='k';
    pointProps.MarkerEdgeColor='y';
    pointProps.MarkerSize=7;
    pointProps.LineStyle='none';
    drawPoint3d(C,pointProps)
    
    anatomicalViewButtons('SRA')
    
    save('template_controls.mat', 'C', 'BE', '-append')
end

%% Calculate weights
if weightSwitch
    addpath(genpath('D:\Biomechanics\General\Code\#external\#Mesh\gptoolbox'))
    
    load('template.mat','template')
    load('template_controls.mat', 'C')
    
    P=1:size(C,1);
    
    % Compute boundary conditions
    [bVertices, bConditions] = boundary_conditions(template.vertices,template.faces, C, P);
    % Compute weights
    weights = biharmonic_bounded(template.vertices, template.faces, bVertices, bConditions, 'OptType','quad');
    % Normalize weights
    weights = weights./repmat(sum(weights,2),1,size(weights,2));
    
    save('template_weights.mat', 'weights')
end

%% Test skinning
addpath(genpath('D:\Biomechanics\General\Code\#external\#Mesh\gptoolbox'))
addpath(genpath('D:\Biomechanics\General\Code\#external\matGeom\matGeom'))

load('template.mat','template')
load('template_controls.mat', 'C','TFM2AFCS')
P=1:size(C,1);
load('template_weights.mat', 'weights')

patchProps.EdgeColor = 'none';
patchProps.FaceColor = [223, 206, 161]/255;
patchProps.FaceAlpha = 0.75;
patchProps.FaceLighting = 'gouraud';
visualizeMeshes(transformPoint3d(template, TFM2AFCS), patchProps)

pointProps.Marker='o';
pointProps.MarkerFaceColor='k';
pointProps.MarkerEdgeColor='y';
pointProps.MarkerSize=7;
pointProps.LineStyle='none';
drawPoint3d(transformPoint3d(C, TFM2AFCS),pointProps)
neckLine = transformLine3d(createLine3d(C(1,:), C(2,:)), TFM2AFCS);
drawLine3d(neckLine)

patchProps.FaceColor = 'r';
patchProps.FaceAlpha = 0.5;

pointProps.MarkerFaceColor='r';
pointProps.MarkerEdgeColor='k';

refVersion = measureFemoralVersionBergmann2016(C(1,:), C(2,:), C(3,:), C(4,:), C(5,:));

Angle = deg2rad(linspace(-25-refVersion,40-refVersion,14)); % neg. ~= retroversion, pos. ~= anteversion
for a=1:length(Angle)
    skinnedTemplate = template;
    new_C = C;
    % Change position of control point by a rotation around the shaft axis
    ROT=createRotation3dLineAngle(createLine3d(C(2,:),C(3,:)), Angle(a));
    new_C(1,:)=transformPoint3d(C(1,:), ROT);
    
    [T,AX,AN,Sm,O] = skinning_transformations(C,P,[],new_C);
    
    % number of handles
    m = numel(P);%+size(BE,1);
    % dimension (2 or 3)
    dim = size(C,2);
    % Extract scale
    TR = zeros(dim,dim+1,m);
    TR(1:dim,1:dim,:) = Sm;
    Sm = reshape(Sm,[dim dim*m])';
    TR(1:dim,dim+1,:) = permute(O-stacktimes(Sm,O),[2 3 1]);
    % Perform scale as linear blend skinning, before translations and rotations
    [new_V] = lbs(template.vertices,TR,weights);
    Q = axisangle2quat(AX,AN);
    % quattrans2udq expect 3D translations, so pad with zeros
    T = [T zeros(size(T,1),1)];
    % convert quaternions and translations into dualquaternions
    DQ = quattrans2udq(Q,T);
    
    skinnedTemplate.vertices = dualquatlbs(new_V,DQ,weights);
    
    vIdx2rem = ismembertol(skinnedTemplate.vertices, template.vertices, 'ByRows', true, 'DataScale', 1e10);
    skinnedTemplate = removeMeshVertices(skinnedTemplate, vIdx2rem);
    
    patch(transformPoint3d(skinnedTemplate, TFM2AFCS), patchProps)
    drawPoint3d(transformPoint3d(new_C(1,:), TFM2AFCS), pointProps)
    
    measureFemoralVersionBergmann2016(new_C(1,:), new_C(2,:), new_C(3,:), new_C(4,:), new_C(5,:))
    
    % Calculate femoral version
    newNeckLine = transformLine3d(createLine3d(new_C(1,:), new_C(2,:)), TFM2AFCS);
    drawLine3d(newNeckLine)
    FV=table(nan(1,1),'VariableNames',{'femoralVersion'},'RowNames',{'Angle [°]'});
    FV{1,1} = measureFemoralVersionBergmann2016(new_C(1,:), new_C(2,:), new_C(3,:), new_C(4,:), new_C(5,:));
    disp(FV)
end

anatomicalViewButtons('RAS')

title('Dual-quaternion skinning with bounded biharmonic weights')