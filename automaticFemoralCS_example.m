clearvars; close all; opengl hardware

% Load example data
load('data\F.mat')

% Select example
Idx = 1;

[fwTFM2AFCS, LMIdx, HJC, LM] = automaticFemoralCS(F(Idx).mesh, F(Idx).side,...
    'definition','MediTEC', 'visu',1, 'verb',0, 'debug',1);

% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';