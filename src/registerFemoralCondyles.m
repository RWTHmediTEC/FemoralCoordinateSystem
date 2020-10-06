function preallignedSubject = registerFemoralCondyles(template, subject)
%REGISTERFEMORALCONDYLES registers the subject's femoral condyles with the 
% template's femoral condyles
%
% AUTHOR: Maximilian C. M. Fischer
% COPYRIGHT (C) 2020 Maximilian C. M. Fischer
% LICENSE: EUPL v1.2
%

% Keep only the distal part of the bones
distalTemplate = template(template(:,1)<(1/2*min(template(:,1))),:);
distalSubject = subject(subject(:,1)<(1/2*min(template(:,1))),:);

% Create Kd-tree for the target points
templateKDTree = createns(distalTemplate);

% Step size of the rotation around a specific axis in degree
StepSize = 5;
Steps = 0:StepSize:360-StepSize;
NoS = length(Steps);
R = cell(1,NoS);
MM = nan(1,NoS);
for i = 1:length(R)
    % Rotation around the x-axis
    R{i} = createRotationOx(deg2rad(Steps(i)));
    [~, DD] = knnsearch(templateKDTree, transformPoint3d(distalSubject, R{i}));
    MM(i) = sum(DD);
end
[~, minIdx] = min(MM);

preallignedSubject = transformPoint3d(subject, R{minIdx});

end
 