function [shaftAxis, shaftAxisIdx] = detectShaftAxis(femur, HJC, varargin)

p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addRequired(p,'femur',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
isPoint3d = @(x) validateattributes(x,{'numeric'},...
    {'nonempty','nonnan','real','finite','size',[1,3]});
addRequired(p,'HJC', isPoint3d);
addParameter(p,'debugVisualization',false,logParValidFunc);

parse(p,femur,HJC,varargin{:});
femur = p.Results.femur;
HJC = p.Results.HJC;
debugVisu = logical(p.Results.debugVisualization);

% Get femur points with maximal distance between each other
[maxDist,maxDistIdx]=max(pdist(femur.vertices));
[maxDistPts(1), maxDistPts(2)] = distanceVectorIndToSub(size(femur.vertices,1), maxDistIdx);
maxDistPts = femur.vertices(maxDistPts,:);
% Check which point is closest to the HJC
DistHJC = distancePoints3d(maxDistPts,HJC);
% Create a line pointing in distal direction
if DistHJC(1)>DistHJC(2)
    maxDistPts=flipud(maxDistPts);
end
maxDistLine = createLine3d(maxDistPts(1,:),maxDistPts(2,:));

maxDistVector = normalizeVector3d(maxDistLine(4:6));
% Create cutting planes 
PROX_CUT = 1/4;
proximalCutPlane = createPlane(...
    maxDistPts(1,:)+maxDistVector*maxDist*PROX_CUT, maxDistVector);
DIST_CUT = PROX_CUT + 1/2;
distalCutPlane = createPlane(...
    maxDistPts(1,:)+maxDistVector*maxDist*DIST_CUT, -maxDistVector);

% The shaft is defined by the mesh between the cutting planes
shaft = cutMeshByPlane(cutMeshByPlane(femur, proximalCutPlane),distalCutPlane);

% Fit ellipsoid to the shaft
shaftEllipsoid = equivalentEllipsoid(shaft.vertices);
% Construct the main shaft axis from the shaft ellipsoid
shaftVector = transformVector3d([1 0 0], eulerAnglesToRotation3d(shaftEllipsoid(7:9)));
shaftAxis = [shaftEllipsoid(1:3) shaftVector];
shaftPlane = createPlane(shaftAxis(1:3), shaftAxis(4:6));
if ~isBelowPlane(HJC,shaftPlane)
    shaftAxis=reverseLine3d(shaftAxis);
end
% Use vertex indices of the mesh to define the shaft axis
shaftAxisIdx = lineToVertexIndices(shaftAxis, femur);
if debugVisu
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [223, 206, 161]/255;
    patchProps.FaceAlpha = 0.5;
    patchProps.FaceLighting = 'gouraud';
    [~, debugAx, debugFig] = visualizeMeshes(femur,patchProps);
    set(debugFig, 'Name', 'Debug figure: shaft axis','NumberTitle','Off')
    mouseControl3d(debugAx)
    % max distance
    drawPoint3d(debugAx, maxDistPts,'MarkerFaceColor','c','MarkerEdgeColor','c')
    drawVector3d(maxDistLine(1:3),maxDistLine(4:6),'c');
    drawPlane3d(debugAx,proximalCutPlane)
    drawPlane3d(debugAx,distalCutPlane)
    % shaft
    patchProps.FaceColor = 'r';
    patch(debugAx, shaft, patchProps)
    % shaft axis
    shaftAxisIdxBased=createLine3d(...
        femur.vertices(shaftAxisIdx(1),:),...
        femur.vertices(shaftAxisIdx(2),:));
    shaftAxisIdxBased(4:6) = normalizeVector3d(shaftAxisIdxBased(4:6));
    drawVector3d(shaftAxisIdxBased(1:3),shaftAxisIdxBased(4:6)*100,'g');
    drawLine3d(debugAx, shaftAxisIdxBased,'k');
end

end