clearvars; close all

addpath(genpath('src'))

% Clone example data
if ~exist('VSD', 'dir')
    try
        !git clone https://github.com/MCM-Fischer/VSDFullBodyBoneModels VSD
        rmdir('VSD/.git', 's')
    catch
        warning([newline 'Clone (or copy) the example data from: ' ...
            'https://github.com/RWTHmediTEC/VSDFullBodyBoneModels' newline 'to: ' ...
            fileparts([mfilename('fullpath'), '.m']) '\VSD' ...
            ' and try again!' newline])
        return
    end
end

% Select subjects of the VSD
subjectXLSX = 'VSD\MATLAB\res\VSD_Subjects.xlsx';
Subjects = readtable(subjectXLSX);
Subjects{2:2:height(Subjects),7} = 'R';
Subjects{1:2:height(Subjects),7} = 'L'; 

for s=1%:size(Subjects, 1)
    id = Subjects{s,1}{1};
    side = Subjects{s,7};
    
    load(['VSD\Bones\' id '.mat'], 'B');
    
    femur = B(ismember({B.name}, ['Femur_' side])).mesh;
    % Remove sesamoid bones if present
    femur = splitMesh(femur, 'mostVertices');
    if s == 22 && strcmp(side, 'R')
        % Subject with a hinged total knee arthroplasty (TKA)
        [TFM2FCS, LM, LMIdx, TFM] = automaticFemoralCS(femur, side,...
            'definition','Wu2002', 'visu',1, 'verb',1, 'debug',0, 'Subject',id, ...
            'minimalRefinement', 1); % Activated "minimalRefinement"
    else
        % All other subjects
        [TFM2FCS, LM, LMIdx, TFM] = automaticFemoralCS(femur, side,...
            'definition','mediTEC', 'visu',1, 'verb',1, 'debug',0, 'Subject',id);
    end
end

% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';