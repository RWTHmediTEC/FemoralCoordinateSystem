function [shaftAxis, shaftAxisIdx] = detectShaftAxis(femur, FHC, varargin)
%DETECTSHAFTAXIS identifies the shaft axis of the femur
% 
%   WORKFLOW
%   - Find the two vertices with the maximum distance between each other.
%	- Based on the length of the vector connecting these two vertices   
%     resect 1/4 at of the femur at the proximal and distal end.
%   - Fit an ellipsoid to the remaining vertices of the shaft and take the 
%     major axis as the shaft axis.
%
%   INPUT
%     femur: femur mesh
%     FHC: (approximate) femoral head center (FHC). Required to 
%     distinguish between proximal and distal end.
%
%   OUTPUT
%     shaftAxis: the femoral shaft axis defined by the origin
%     (approximately in the middle of the shaft) and a vector pointing in
%     the distal direction.
%     shaftAxisIdx: the shaft axis defined by two vertex indices of the
%     femur.
%
% AUTHOR: Maximilian C. M. Fischer
% COPYRIGHT (C) 2020-2023 Maximilian C. M. Fischer
% LICENSE: EUPL v1.2
%

p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addRequired(p,'femur',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
isPoint3d = @(x) validateattributes(x,{'numeric'},...
    {'nonempty','nonnan','real','finite','size',[1,3]});
addRequired(p,'FHC', isPoint3d);
addParameter(p,'debugVisualization',false,logParValidFunc);

parse(p,femur,FHC,varargin{:});
femur = p.Results.femur;
FHC = p.Results.FHC;
debugVisu = logical(p.Results.debugVisualization);

% Get femur points with maximal distance between each other
[maxDist,maxDistIdx]=max(pdist(femur.vertices));
[maxDistPts(1), maxDistPts(2)] = distanceVectorIndToSub(size(femur.vertices,1), maxDistIdx);
maxDistPts = femur.vertices(maxDistPts,:);
% Check which point is closest to the FHC
DistFHC = distancePoints3d(maxDistPts,FHC);
% Create a line pointing in distal direction
if DistFHC(1)>DistFHC(2)
    maxDistPts = flipud(maxDistPts);
end
maxDistLine = createLine3d(maxDistPts(1,:),maxDistPts(2,:));

maxDistVector = normalizeVector3d(maxDistLine(4:6));
% Create cutting planes 
PROX_CUT = 1/4;
proxCutPlane = createPlane(...
    maxDistPts(1,:)+maxDistVector*maxDist*PROX_CUT, maxDistVector);
DIST_CUT = PROX_CUT + 1/2;
distCutPlane = createPlane(...
    maxDistPts(1,:)+maxDistVector*maxDist*DIST_CUT, -maxDistVector);

% The shaft is defined by the mesh between the cutting planes
[shaft, ~, proxPart] = cutMeshByPlane(femur, proxCutPlane);
[shaft, ~, distPart] = cutMeshByPlane(shaft,distCutPlane);

% Fit ellipsoid to the shaft
shaftEllipsoid = equivalentEllipsoid(shaft.vertices);
% Construct the main shaft axis from the shaft ellipsoid
shaftVector = transformVector3d([1 0 0], eulerAnglesToRotation3d(shaftEllipsoid(7:9)));
shaftAxis = [shaftEllipsoid(1:3) shaftVector];
shaftPlane = createPlane(shaftAxis(1:3), shaftAxis(4:6));
% Shaft axis should point in distal direction
if ~isBelowPlane(FHC,shaftPlane)
    shaftAxis=reverseLine3d(shaftAxis);
end
% Use vertex indices of the mesh to define the shaft axis
shaftAxisIdx = lineToVertexIndices(shaftAxis, femur);

if debugVisu
    patchProps.FaceAlpha = 0.5;
    [~, debugAx, debugFig] = visualizeMeshes([proxPart,distPart], patchProps);
    set(debugFig, 'Name', 'Shaft axis (debug figure)','NumberTitle','Off')
    axis(debugAx, 'off','tight')
    mouseControl3d(debugAx)
    % Draw max. distance line
    drawPoint3d(debugAx, maxDistPts, 'MarkerFaceColor','c', 'MarkerEdgeColor','c')
    drawEdge3d(maxDistPts, 'Marker','o', 'MarkerEdgeColor','k', 'MarkerFaceColor','k',...
        'Color','k', 'LineWidth',2);
    planeProps.FaceColor = 'k';
    planeProps.EdgeColor = 'k';
    planeProps.FaceAlpha = 0.5;
    drawPlatform(debugAx,[intersectLinePlane(shaftAxis, proxCutPlane), proxCutPlane(4:9)], 70, planeProps)
    drawPlatform(debugAx,[intersectLinePlane(shaftAxis, distCutPlane), distCutPlane(4:9)], 70, planeProps)
    % Draw shaft
    patchProps.FaceColor = 'r';
    patchProps.EdgeColor = 'none';
    patch(debugAx, shaft, patchProps)
    % Draw shaft axis
    drawLine3d(debugAx, shaftAxis, 'Color','r', 'LineWidth',2)
    drawArrow3d(debugAx, shaftAxis(1:3), normalizeVector3d(shaftAxis(4:6))*50,'r')
    % Draw FHC
    drawPoint3d(FHC, 'Marker','o', 'MarkerEdgeColor','g', 'MarkerFaceColor','g', 'MarkerSize',10)
    
    % For publication
    % OffsetPoint = projPointOnLine3d(FHC, shaftAxis);
    % OffsetAxis = normalizeLine3d([OffsetPoint FHC-OffsetPoint]);
    % camTar = mean(femur.vertices);
    % set(gca,'CameraTarget',camTar);
    % camVec = normalizeVector3d(crossProduct3d(OffsetAxis(4:6),normalizeVector3d(shaftAxis(4:6))));
    % set(gca,'CameraPosition',camTar+[0 0 100]+camVec*1000);
    % set(gca,'CameraUpVector',OffsetAxis(4:6));
    % set(gca,'CameraViewAngle',15)
    % set(gcf,'GraphicsSmoothing','off')
    % export_fig('Figure3', '-tif', '-r300')
end

end