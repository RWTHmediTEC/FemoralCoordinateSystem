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

% Select subjects of the VSD
Subjects = [1 9 13 19 23 24 27 35 36 42 46 49 50 55 56 57 61 62 64 66];
Subjects = arrayfun(@(x) ['z' num2str(x, '%03i')], Subjects', 'uni',0);
Subjects(1:2:20,2) = {'L'}; Subjects(2:2:20,2) = {'R'};

for s=1%:size(Subjects, 1)
    name = Subjects{s,1};
    side = Subjects{s,2};
    
    load(['VSD\Bones\' name '.mat'], 'B');
    
    femur = B(ismember({B.name}, ['Femur_' side])).mesh;
    [TFM2FCS, LM, LMIdx, TFM] = automaticFemoralCS(femur, side,...
        'definition','MediTEC', 'visu',1, 'verb',1, 'debug',0, 'Subject',name);
    
end

% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';