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
sides = {'R','L'};

for s=1%:size(Subjects, 1)
    load(['VSD\Bones\' Subjects.Number{s} '.mat'],'B');
    
    sideIdx = randi(2);
    
    [fwTFM2AFCS, LMIdx, HJC, LM] = automaticFemoralCS(B(sideIdx+3).mesh, sides{sideIdx},...
        'definition','MediTEC', 'visu',1, 'verb',0, 'debug',1);
    
end

% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';