function TFM = Bergmann2016(femur, side, HJC, LMIdx, varargin)
%BERGMANN2016 calculates the femoral coordinate system used in Bergmann2016
%
% 2016 - Bergmann et al. - Standardized Loads Acting in Hip Implants:
%   "The origin of this coordinate system is located at the centre of the 
%   femoral head. The +z axis points upward and is defined by the line 
%   connecting the two points where the curved femoral mid-line intersected 
%   with the neck axis (P1) and where it passes the intercondylar notch 
%   (P2). The +x axis points laterally and is oriented parallel to the 
%   proximal contour of the condyles. The +y axis points in the anterior 
%   direction."
%
% AUTHOR: Maximilian C. M. Fischer
% COPYRIGHT (C) 2020 Maximilian C. M. Fischer
% LICENSE: EUPL v1.2
%

% inputs
p = inputParser;
addRequired(p,'femur',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addRequired(p,'side',@(x) any(validatestring(x,{'R','L'})));
addOptional(p,'visualization',true,@islogical);
parse(p,femur,side,varargin{:});

femur = p.Results.femur;
visu = p.Results.visualization;

% Landmarks
MPC = femur.vertices(LMIdx.MedialPosteriorCondyle,:);
LPC = femur.vertices(LMIdx.LateralPosteriorCondyle,:);
ICN = femur.vertices(LMIdx.IntercondylarNotch,:);
NeckAxis = createLine3d(femur.vertices(LMIdx.NeckAxis(1),:),femur.vertices(LMIdx.NeckAxis(2),:));
ShaftAxis = createLine3d(femur.vertices(LMIdx.ShaftAxis(1),:),femur.vertices(LMIdx.ShaftAxis(2),:));

%% Construction of P1
[~, P1, ~] = distanceLines3d(NeckAxis, ShaftAxis);
StraightFemurAxis = createLine3d(ICN, P1);

%% inital transformation
% Connection of the most posterior points of the condyles
PosteriorCondyleAxis = createLine3d(MPC, LPC);

Z = normalizeVector3d(StraightFemurAxis(4:6));
Y = normalizeVector3d(crossProduct3d(StraightFemurAxis(4:6), PosteriorCondyleAxis(4:6)));
X = normalizeVector3d(crossProduct3d(Y, Z));

TFM = [[X;Y;Z],[0 0 0]'; [0 0 0 1]]*createTranslation3d(-HJC);
% If it is a left femur, rotate 180° around the Z axis
if strcmp(side, 'L')
    TFM = createRotationOz(pi)*TFM;
end

%% visualization
if visu
    % The femur in the AFCS
    femurCS = transformPoint3d(femur, TFM);
    
    % Patch properties
    patchProps.FaceAlpha = 0.75;
    [~, axH, figH] = visualizeMeshes(femurCS, patchProps);
    figH.NumberTitle = 'off';
    figH.Name = 'Bergmann2016 coordinate system';
    drawAxis3d(axH,35,1.5)
    
    % Landmarks
    LM_Idx = struct2cell(LMIdx);
    LM_Idx = cell2mat(LM_Idx(structfun(@(x) length(x) == 1, LMIdx)));
    drawPoint3d(femurCS.vertices(LM_Idx, :),...
        'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
    drawPoint3d(transformPoint3d([P1; ICN], TFM),...
        'MarkerFaceColor','k','MarkerEdgeColor','k')
    % Axes
    edgeProps.LineStyle='-';
    edgeProps.Color='k';
    
    drawEdge3d(axH,clipLine3d(transformLine3d(StraightFemurAxis, TFM),...
        boundingBox3d(femurCS.vertices)),edgeProps)
    drawEdge3d(axH,...
        femurCS.vertices(LMIdx.NeckAxis(1),:),...
        femurCS.vertices(LMIdx.NeckAxis(2),:), edgeProps);
    edgeProps.Marker='o';
    edgeProps.MarkerEdgeColor='k';
    edgeProps.MarkerFaceColor='k';
    
    PC=[femurCS.vertices(LMIdx.MedialPosteriorCondyle,:);...
        femurCS.vertices(LMIdx.LateralPosteriorCondyle,:)];
    drawEdge3d(axH,PC(1,:),PC(2,:), edgeProps)
    
    text(axH,PC(:,1),PC(:,2),PC(:,3),{'MPC';'LPC'})
    
    anatomicalViewButtons(axH,'RAS')
end

end