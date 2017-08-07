function Prealligned_source = roughPreRegistration(target, source)

% Create Kd-tree for the target points
targetKDTRee=createns(target);

% A 180° rotation around the y-axis
R180Y = [cos(pi) 0 sin(pi); 0 1 0; -sin(pi) 0 cos(pi)];

% Step size of the rotation around a specific axis in degree
StepSize = 10;
Steps = 0:StepSize:360-StepSize;
NoS = length(Steps);
R1 = cell(1,NoS);
MM1=nan(1,NoS);
MM2=nan(1,NoS);

for i = 1:length(R1)
    theta = degtorad(Steps(i));
    % switch axis  % create rotation
    % case 'X' % around x-axis
    R1{i} = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
    % case 'Y' % around y-axis
    % R{i} = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    % case 'Z' % around z-axis
    % R{i} = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
    % end
    [~, DD] = knnsearch(targetKDTRee,source*R1{i});
    MM1(i) = sum(DD);
    % Also try upside down
    [~, DD] = knnsearch(targetKDTRee,source*R1{i}*R180Y);
    MM2(i) = sum(DD);
end
[min1, I1]=min(MM1);
[min2, I2]=min(MM2);

if min1<min2
    Prealligned_source = source*R1{I1};
else
    Prealligned_source = source*(R1{I2}*R180Y);
end
 