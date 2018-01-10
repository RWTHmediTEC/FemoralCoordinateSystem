clearvars; close all; opengl hardware
% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';

% Load example data
load('data\F.mat')

% Select example
Idx=4;

[fwTFM2AFCS, LMIdx] = automaticFemoralCS(F(Idx).mesh, F(Idx).side,...
    'definition','TabletopMediTEC','vis', true);