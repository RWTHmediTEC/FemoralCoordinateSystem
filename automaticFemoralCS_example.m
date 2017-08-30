clearvars; close all; opengl hardware

[List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
List.f = List.f'; List.p = List.p';

load('data\F.mat')

Idx=4;

[fwTFM2AFCS, LMIdx] = automaticFemoralCS(F(Idx).mesh, F(Idx).side,...
    'definition','WuBergmannComb','vis', true,'iter',3);