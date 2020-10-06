function preallignedTemplate = registerProximalFemur(template, subject)
%REGISTERPROXIMALFEMUR adjusts the proximal femur of the template to allign
% it with the neck and head of the subject.
%
%   WORKFLOW
%   1. Adjust the the femoral version (aka torsion)
%   2. Adjust the neck shaft angle (NSA, aka CCD angle)
%   3. Adjust the neck length (NL)
%
% AUTHOR: Maximilian C. M. Fischer
% COPYRIGHT (C) 2020 Maximilian C. M. Fischer
% LICENSE: EUPL v1.2
%

% Create kd-tree for the subject points
subjectKDTree = createns(subject);

load('template_controls.mat', 'C')
P=1:size(C,1);
load('template_weights.mat', 'weights')

refVersion = measureFemoralVersionBergmann2016(C(1,:), C(2,:), C(3,:), C(4,:), C(5,:));
Angle = deg2rad(-25-refVersion:2.5:40-refVersion); % neg. ~= retroversion, pos. ~= anteversion
rotAxis_Version = createLine3d(C(2,:),C(3,:));

MM = nan(1,length(Angle));
tempTemplate = cell(1,length(Angle));
C_ver = nan([size(C),length(Angle)]);
for a=1:length(Angle)
    C_ver(:,:,a) = C;
    % Change position of control point by a rotation around the shaft axis
    ROT = createRotation3dLineAngle(rotAxis_Version, Angle(a));
    C_ver(1,:,a) = transformPoint3d(C(1,:), ROT);
    
    [T,AX,AN,Sm,O] = skinning_transformations(C,P,[],C_ver(:,:,a));
    
    % number of handles
    m = numel(P);
    % dimension (2 or 3)
    dim = size(C,2);
    % Extract scale
    TR = zeros(dim,dim+1,m);
    TR(1:dim,1:dim,:) = Sm;
    Sm = reshape(Sm,[dim dim*m])';
    TR(1:dim,dim+1,:) = permute(O-stacktimes(Sm,O),[2 3 1]);
    % Perform scale as linear blend skinning, before translations and rotations
    [new_V] = lbs(template,TR,weights);
    Q = axisangle2quat(AX,AN);
    % convert quaternions and translations into dualquaternions
    DQ = quattrans2udq(Q,T);
    
    tempTemplate{a} = dualquatlbs(new_V,DQ,weights);
    [~, DD] = knnsearch(subjectKDTree, tempTemplate{a});
    MM(a) = sum(DD);
end
% Find skinned template with minimum deviation to the subject
[~, minIdx] = min(MM);


%% Adjust neck shaft angle (NSA)
C_NSA = C_ver(:,:,minIdx);
refNSA = rad2deg(vectorAngle3d(C_NSA(1,:)-C_NSA(2,:),C_NSA(3,:)-C_NSA(2,:)));
Angle = deg2rad(120-refNSA:2.5:140-refNSA);
rotAxis_NSA = [C_NSA(2,:), crossProduct3d(C_NSA(1,:)-C_NSA(2,:),C_NSA(2,:)-C_NSA(3,:))];

MM = nan(1,length(Angle));
tempTemplate = cell(1,length(Angle));
C_NSA = nan([size(C),length(Angle)]);
for a=1:length(Angle)
    C_NSA(:,:,a) = C_ver(:,:,minIdx);
    % Change position of the FHC by a rotation around the orthogonal of
    % shaft and neck axis
    ROT = createRotation3dLineAngle(rotAxis_NSA, Angle(a));
    C_NSA(1,:,a) = transformPoint3d(C_NSA(1,:,a), ROT);
    
    [T,AX,AN,Sm,O] = skinning_transformations(C,P,[],C_NSA(:,:,a));
    
    % number of handles
    m = numel(P);
    % dimension (2 or 3)
    dim = size(C,2);
    % Extract scale
    TR = zeros(dim,dim+1,m);
    TR(1:dim,1:dim,:) = Sm;
    Sm = reshape(Sm,[dim dim*m])';
    TR(1:dim,dim+1,:) = permute(O-stacktimes(Sm,O),[2 3 1]);
    % Perform scale as linear blend skinning, before translations and rotations
    [new_V] = lbs(template,TR,weights);
    Q = axisangle2quat(AX,AN);
    % convert quaternions and translations into dualquaternions
    DQ = quattrans2udq(Q,T);
    
    tempTemplate{a} = dualquatlbs(new_V,DQ,weights);
    [~, DD] = knnsearch(subjectKDTree, tempTemplate{a});
    MM(a) = sum(DD);
end
% Find skinned template with minimum deviation to the subject
[~, minIdx] = min(MM);

%% Adjust neck length (NL)
C_NL = C_NSA(:,:,minIdx);
refNL = distancePoints3d(C_NL(1,:), C_NL(2,:));
Dist = 40-refNL:2:60-refNL;

MM = nan(1,length(Dist));
tempTemplate = cell(1,length(Dist));
C_NL = nan([size(C),length(Dist)]);
for a=1:length(Dist)
    C_NL(:,:,a) = C_NSA(:,:,minIdx);
    C_NL(1,:,a) = C_NL(1,:,a) - normalizeVector3d(C_NL(2,:,a)-C_NL(1,:,a))*Dist(a);
    
    [T,AX,AN,Sm,O] = skinning_transformations(C,P,[],C_NL(:,:,a));
    
    % number of handles
    m = numel(P);
    % dimension (2 or 3)
    dim = size(C,2);
    % Extract scale
    TR = zeros(dim,dim+1,m);
    TR(1:dim,1:dim,:) = Sm;
    Sm = reshape(Sm,[dim dim*m])';
    TR(1:dim,dim+1,:) = permute(O-stacktimes(Sm,O),[2 3 1]);
    % Perform scale as linear blend skinning, before translations and rotations
    [new_V] = lbs(template,TR,weights);
    Q = axisangle2quat(AX,AN);
    % convert quaternions and translations into dualquaternions
    DQ = quattrans2udq(Q,T);
    
    tempTemplate{a} = dualquatlbs(new_V,DQ,weights);
    [~, DD] = knnsearch(subjectKDTree, tempTemplate{a});
    MM(a) = sum(DD);
end
% Find skinned template with minimum deviation to the subject
[~, minIdx] = min(MM);

preallignedTemplate = tempTemplate{minIdx};

end