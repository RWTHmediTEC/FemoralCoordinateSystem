clearvars; close all; opengl hardware

% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';

load('data\F.mat')

Idx=1;

[fwTFM2AFCS, LMIdx] = automaticFemoralCS(F(Idx).mesh, F(Idx).side,...
    'definition','Bergmann2016','vis', true);