function Prealligned_source = roughPreRegistration(target, source)

% Create Kd-tree for the target points
targetKDTRee=createns(target);

% Position source to the target in x direction
xPosition = max(target(:,1))-max(source(:,1));
TFM2xPosition=eye(4); TFM2xPosition(1,4)=xPosition;
source1=transformPoint3d(source, TFM2xPosition);

% A 180° rotation around the y-axis (upside down)
R180Y = [cos(pi) 0 sin(pi); 0 1 0; -sin(pi) 0 cos(pi)];
source2=transformPoint3d(source, R180Y);
% Position source to the target in x direction
xPosition = max(target(:,1))-max(source2(:,1));
TFM2xPosition=eye(4); TFM2xPosition(1,4)=xPosition;
source2=transformPoint3d(source2, TFM2xPosition);

% Step size of the rotation around a specific axis in degree
StepSize = 10;
Steps = 0:StepSize:360-StepSize;
NoS = length(Steps);
RX = cell(1,NoS);
MM1=nan(1,NoS);
MM2=nan(1,NoS);

for i = 1:length(RX)
    theta = degtorad(Steps(i));
    % switch axis  % create rotation
    % case 'X' % around x-axis
    RX{i} = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
    % case 'Y' % around y-axis
    % R{i} = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    % case 'Z' % around z-axis
    % R{i} = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
    % end
    [~, DD1] = knnsearch(targetKDTRee,source1*RX{i});
    MM1(i) = sum(DD1);
    % Also try upside down
    [~, DD2] = knnsearch(targetKDTRee,source2*RX{i});
    MM2(i) = sum(DD2);
end
[min1, I1]=min(MM1);
[min2, I2]=min(MM2);

if min1<min2
    Prealligned_source = source1*RX{I1};
else
    Prealligned_source = source2*RX{I2};
end
 