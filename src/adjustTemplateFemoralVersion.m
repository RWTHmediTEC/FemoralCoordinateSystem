function preallignedTarget = adjustTemplateFemoralVersion(target, source)

% Create Kd-tree for the source points
sourceKDTRee=createns(source);

load('template_controls.mat', 'C')
P=1:size(C,1); %#ok<NODEF>
load('template_weights.mat', 'weights')

NoA=8;
Angle = deg2rad(linspace(-25,37,NoA)); % neg. ~= retroversion, pos. ~= anteversion
MM=nan(1,NoA);
tempTargets=cell(1,NoA);
for a=1:length(Angle)
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
    [new_V] = lbs(target,TR,weights);
    Q = axisangle2quat(AX,AN);
    % quattrans2udq expect 3D translations, so pad with zeros
    T = [T zeros(size(T,1),1)];
    % convert quaternions and translations into dualquaternions
    DQ = quattrans2udq(Q,T);
    
    tempTargets{a}=dualquatlbs(new_V,DQ,weights);
    
    [~, DD] = knnsearch(sourceKDTRee, tempTargets{a});
    MM(a) = sum(DD);
end
% Find skinned target with minimum deviation to the source
[~, minIdx]=min(MM);
preallignedTarget=tempTargets{minIdx};

end