function preallignedTarget = adjustTemplateFemoralVersion(target, source)

% Create Kd-tree for the source points
sourceKDTRee=createns(source);

load('template_controls.mat', 'C')
P=1:size(C,1); 
load('template_weights.mat', 'weights')

refVersion = measureFemoralVersionBergmann2016(C(1,:), C(2,:), C(3,:), C(4,:), C(5,:));
Angle = deg2rad(-25-refVersion:2.5:40-refVersion); % neg. ~= retroversion, pos. ~= anteversion
rotAxis_Version = createLine3d(C(2,:),C(3,:));

MM=nan(1,length(Angle));
tempTargets=cell(1,length(Angle));
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
    [new_V] = lbs(target,TR,weights);
    Q = axisangle2quat(AX,AN);
    % quattrans2udq expect 3D translations, so pad with zeros
    T = [T zeros(size(T,1),1)];
    % convert quaternions and translations into dualquaternions
    DQ = quattrans2udq(Q,T);
    
    tempTargets{a} = dualquatlbs(new_V,DQ,weights);
    [~, DD] = knnsearch(sourceKDTRee, tempTargets{a});
    MM(a) = sum(DD);
end
% Find skinned target with minimum deviation to the source
[~, minIdx]=min(MM);


% Adjust CCD angle
C_CCD = C_ver(:,:,minIdx);
refCCD = rad2deg(vectorAngle3d(C_CCD(1,:)-C_CCD(2,:),C_CCD(3,:)-C_CCD(2,:)));
Angle = deg2rad(120-refCCD:2.5:140-refCCD);
rotAxis_CCD = [C_CCD(2,:), crossProduct3d(C_CCD(1,:)-C_CCD(2,:),C_CCD(2,:)-C_CCD(3,:))];

MM=nan(1,length(Angle));
tempTargets=cell(1,length(Angle));
C_CCD = nan([size(C),length(Angle)]);
for a=1:length(Angle)
    C_CCD(:,:,a) = C_ver(:,:,minIdx);
    % Change position of the FHC by a rotation around the orthogonal of 
    % shaft and neck axis
    ROT = createRotation3dLineAngle(rotAxis_CCD, Angle(a));
    C_CCD(1,:,a) = transformPoint3d(C_CCD(1,:,a), ROT);
    
    [T,AX,AN,Sm,O] = skinning_transformations(C,P,[],C_CCD(:,:,a));
    
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
    [new_V] = lbs(target,TR,weights);
    Q = axisangle2quat(AX,AN);
    % quattrans2udq expect 3D translations, so pad with zeros
    T = [T zeros(size(T,1),1)];
    % convert quaternions and translations into dualquaternions
    DQ = quattrans2udq(Q,T);
    
    tempTargets{a} = dualquatlbs(new_V,DQ,weights);
    [~, DD] = knnsearch(sourceKDTRee, tempTargets{a});
    MM(a) = sum(DD);
end
% Find skinned target with minimum deviation to the source
[~, minIdx]=min(MM);

% Adjust neck length
C_len = C_CCD(:,:,minIdx);
refLength = distancePoints3d(C_len(1,:), C_len(2,:));
Dist = 40-refLength:2:60-refLength;

MM=nan(1,length(Dist));
tempTargets=cell(1,length(Dist));
C_len = nan([size(C),length(Dist)]);
for a=1:length(Dist)
    C_len(:,:,a) = C_CCD(:,:,minIdx);
    C_len(1,:,a) = C_len(1,:,a) - normalizeVector3d(C_len(2,:,a)-C_len(1,:,a))*Dist(a);
    
    [T,AX,AN,Sm,O] = skinning_transformations(C,P,[],C_len(:,:,a));
    
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
    [new_V] = lbs(target,TR,weights);
    Q = axisangle2quat(AX,AN);
    % quattrans2udq expect 3D translations, so pad with zeros
    T = [T zeros(size(T,1),1)];
    % convert quaternions and translations into dualquaternions
    DQ = quattrans2udq(Q,T);
    
    tempTargets{a} = dualquatlbs(new_V,DQ,weights);
    [~, DD] = knnsearch(sourceKDTRee, tempTargets{a});
    MM(a) = sum(DD);
end
% Find skinned target with minimum deviation to the source
[~, minIdx]=min(MM);

preallignedTarget=tempTargets{minIdx};

end