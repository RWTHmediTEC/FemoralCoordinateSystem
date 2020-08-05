clearvars; close all; opengl hardware

% Load example data
load('data\F.mat')

% Select example
Idx=4;

[fwTFM2AFCS, LMIdx, HJC, LM] = automaticFemoralCS(F(Idx).mesh, F(Idx).side,...
    'definition','MediTEC','vis',1, 'verb',0, 'debug',0);

% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';