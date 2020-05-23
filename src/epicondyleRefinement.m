function [MEC_USP, LEC_USP] = epicondyleRefinement(distalFemurUSP, CEA, ...
    MEC_morph_USP, LEC_morph_USP, varargin)

%% Parse inputs
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addRequired(p,'distalFemurUSP',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addParameter(p,'visualization',false,logParValidFunc);
parse(p,distalFemurUSP,varargin{:});
visu = logical(p.Results.visualization);

% Get min/max indices of the epicodyles in USP system
[~, ecMinMaxIdxUSP(1)]=min(distalFemurUSP.vertices(:,3)); % MEC
[~, ecMinMaxIdxUSP(2)]=max(distalFemurUSP.vertices(:,3)); % LEC
MEC_max_USP=distalFemurUSP.vertices(ecMinMaxIdxUSP(1),:);
LEC_max_USP=distalFemurUSP.vertices(ecMinMaxIdxUSP(2),:);
% Get the CEA indices of the epicodyles in USP system
CEA_intersectPts=intersectLineMesh3d(CEA, distalFemurUSP.vertices, distalFemurUSP.faces);
[~, ecCEA_IdxUSP(1)]=min(CEA_intersectPts(:,3)); % MEC
[~, ecCEA_IdxUSP(2)]=max(CEA_intersectPts(:,3)); % LEC
MEC_CEA_USP=CEA_intersectPts(ecCEA_IdxUSP(1),:);
LEC_CEA_USP=CEA_intersectPts(ecCEA_IdxUSP(2),:);
% Sanity check in case of osteophytes
MEC_CEA_morph_Dist = distancePoints3d(MEC_morph_USP, MEC_CEA_USP);
if distancePoints3d(MEC_max_USP, MEC_CEA_USP) > MEC_CEA_morph_Dist
    MEC_morph_sphere = [MEC_morph_USP, MEC_CEA_morph_Dist];
    MEC_CEA_sphere = [MEC_CEA_USP, MEC_CEA_morph_Dist];
    MEC_cand=clipPoints3d(distalFemurUSP.vertices, MEC_morph_sphere, 'shape', 'sphere');
    MEC_cand= [MEC_cand; clipPoints3d(distalFemurUSP.vertices, MEC_CEA_sphere, 'shape', 'sphere')];
    [~, tempCandIdx]=min(MEC_cand(:,3));
    MEC_USP=MEC_cand(tempCandIdx,:);
else
    % Keep max for MEC
    MEC_USP=MEC_max_USP;
end
LEC_CEA_morph_Dist = distancePoints3d(LEC_morph_USP, LEC_CEA_USP);
if distancePoints3d(LEC_max_USP, LEC_CEA_USP) > LEC_CEA_morph_Dist
    LEC_morph_sphere = [LEC_morph_USP, LEC_CEA_morph_Dist];
    LEC_CEA_sphere = [LEC_CEA_USP, LEC_CEA_morph_Dist];
    LEC_cand=clipPoints3d(distalFemurUSP.vertices, LEC_morph_sphere, 'shape', 'sphere');
    LEC_cand= [LEC_cand; clipPoints3d(distalFemurUSP.vertices, LEC_CEA_sphere, 'shape', 'sphere')];
    [~, tempCandIdx]=max(LEC_cand(:,3));
    LEC_USP=LEC_cand(tempCandIdx,:);
else
    % Keep max for MEC
    LEC_USP=LEC_max_USP;
end

if visu
    [~, axH,]=visualizeMeshes(distalFemurUSP);
    drawLine3d(axH, CEA)
    drawPoint3d(axH, [MEC_morph_USP; LEC_morph_USP],...
        'MarkerFaceColor','r','MarkerEdgeColor','r');
    drawPoint3d(axH, [MEC_max_USP; LEC_max_USP],...
        'MarkerFaceColor','g','MarkerEdgeColor','g');
    drawPoint3d(axH, [MEC_CEA_USP; LEC_CEA_USP],...
        'MarkerFaceColor','b','MarkerEdgeColor','b');
    if exist('MEC_cand', 'var')
        drawPoint3d(axH, MEC_cand,'MarkerEdgeColor','y');
    end
    if exist('LEC_cand', 'var')
        drawPoint3d(axH, LEC_cand,'MarkerEdgeColor','y');
    end
    drawPoint3d(axH, [MEC_USP; LEC_USP],...
        'MarkerFaceColor','k','MarkerEdgeColor','k');
    anatomicalViewButtons(axH, 'ASR')
end

end