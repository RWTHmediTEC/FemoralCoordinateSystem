function Prealligned_source = regFemoralCondyles(target, source)

% Keep only the distal part of the bones
distTarget=target(target(:,1)<(1/2*min(target(:,1))),:);
distSource=source(source(:,1)<(1/2*min(target(:,1))),:);

% Create Kd-tree for the target points
targetKDTRee=createns(distTarget);

% Step size of the rotation around a specific axis in degree
StepSize = 5;
Steps = 0:StepSize:360-StepSize;
NoS = length(Steps);
R = cell(1,NoS);
MM=nan(1,NoS);
for i = 1:length(R)
    theta = deg2rad(Steps(i));
    % switch axis  % create rotation
    % case 'X' % around x-axis
    R{i} = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
    % case 'Y' % around y-axis
    % R{i} = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    % case 'Z' % around z-axis
    % R{i} = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
    % end
    [~, DD] = knnsearch(targetKDTRee,distSource*R{i});
    MM(i) = sum(DD);
end
[~, I]=min(MM);

Prealligned_source = source*R{I};

end
 