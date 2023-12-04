function preallignedSubject = roughPreRegistration(template, subject)
%ROUGHPREREGISTRATION registers roughly subject femur to template femur
%
% AUTHOR: Maximilian C. M. Fischer
% COPYRIGHT (C) 2020-2023 Maximilian C. M. Fischer
% LICENSE: EUPL v1.2
%

% Create kd-tree for the template points
templateKDTree = createns(template);

% Position source to the template in x direction
xPosition = max(template(:,1)) - max(subject(:,1));
TFM2xPosition = eye(4); TFM2xPosition(1,4) = xPosition;
tempSubject1 = transformPoint3d(subject, TFM2xPosition);

% A 180° rotation around the y-axis (upside down)
R180Y = [cos(pi) 0 sin(pi); 0 1 0; -sin(pi) 0 cos(pi)];
tempSubject2 = transformPoint3d(subject, R180Y);
% Position subject to the template in x direction
xPosition = max(template(:,1)) - max(tempSubject2(:,1));
TFM2xPosition = eye(4); TFM2xPosition(1,4) = xPosition;
tempSubject2 = transformPoint3d(tempSubject2, TFM2xPosition);

% Step size of the rotation around a specific axis in degree
StepSize = 10;
Steps = 0:StepSize:360-StepSize;
NoS = length(Steps);
RX = cell(1,NoS);
MM1 = nan(1,NoS);
MM2 = nan(1,NoS);

for i = 1:length(RX)
    % Rotation around the x-axis
    RX{i} = createRotationOx(deg2rad(Steps(i)));
    [~, DD1] = knnsearch(templateKDTree, transformPoint3d(tempSubject1, RX{i}));
    MM1(i) = sum(DD1);
    % Also try upside down
    [~, DD2] = knnsearch(templateKDTree, transformPoint3d(tempSubject2, RX{i}));
    MM2(i) = sum(DD2);
end
[min1, I1] = min(MM1);
[min2, I2] = min(MM2);

if min1<min2
    preallignedSubject = transformPoint3d(tempSubject1, RX{I1});
else
    preallignedSubject = transformPoint3d(tempSubject2, RX{I2});
end