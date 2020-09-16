clearvars; close all; opengl hardware

% Clone example data
if ~exist('VSD', 'dir')
    try
        !git clone https://github.com/RWTHmediTEC/VSDFullBodyBoneModels VSD
        rmdir('VSD/.git', 's')
    catch
        warning([newline 'Clone (or copy) the example data from: ' ...
            'https://github.com/RWTHmediTEC/VSDFullBodyBoneModels' newline 'to: ' ...
            fileparts([mfilename('fullpath'), '.m']) '\VSD' ...
            ' and try again!' newline])
        return
    end
end

% Load subject names
load('VSD\MATLAB\res\VSD_Subjects.mat', 'Subjects')
Subjects = table2cell(Subjects);
Subjects(1:2:20,4) = {'L'}; Subjects(2:2:20,4) = {'R'};

for s=1%:size(Subjects, 1)
    load(['VSD\Bones\' Subjects{s,1} '.mat'], 'B');
    
    femur = B(ismember({B.name}, ['Femur_' Subjects{s,4}])).mesh;
    [fwTFM2AFCS, LMIdx, HJC, LM] = automaticFemoralCS(femur, Subjects{s,4},...
        'definition','MediTEC', 'visu',1, 'verb',0, 'debug',1,'Subject', Subjects{s,1});
    
end

% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';