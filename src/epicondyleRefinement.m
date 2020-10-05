function [MEC_USP, LEC_USP] = epicondyleRefinement(distalFemurUSP, CEA, MEC_map, LEC_map, varargin)
%EPICONDYLEREFINEMENT refines the position of the epicondyles
%
% INPUT
%   distalFemurUSP: Distal part of the femur in the unified sagittal plane
%                   coordinate system (USP CS, see USP.m).
%   CEA: The central elliptic axis in the USP CS.
%   MEC_map & LEC_map: Mapped epicondyles from the template in the USP CS.
%
% OUTPUT
%   MEC_USP & LEC_USP: Final epicondyles in the USP CS.
%
% AUTHOR: Maximilian C. M. Fischer
% COPYRIGHT (C) 2020 Maximilian C. M. Fischer
% LICENSE: EUPL v1.2
%

%% Parse inputs
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addRequired(p,'distalFemurUSP',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addParameter(p,'visualization',false,logParValidFunc);
parse(p,distalFemurUSP,varargin{:});
visu = logical(p.Results.visualization);

% Get the CEA intersections with the epicodyles in USP system
CEA_intersectPts = intersectLineMesh3d(CEA, distalFemurUSP.vertices, distalFemurUSP.faces);

%% MEC
% Most medial vertex of the distal femur
[~, MEC_max_Idx] = min(distalFemurUSP.vertices(:,3)); % MEC
MEC_max = distalFemurUSP.vertices(MEC_max_Idx,:);
% Medial intersection of the CEA with distal femur
[~, MEC_CEA_Idx] = min(CEA_intersectPts(:,3));
MEC_CEA = CEA_intersectPts(MEC_CEA_Idx,:);
% Sanity check in case of osteophytes
MEC_CEA_map_Dist = distancePoints3d(MEC_map, MEC_CEA);
MEC_CEA_max_Dist = distancePoints3d(MEC_max, MEC_CEA);
if MEC_CEA_max_Dist > MEC_CEA_map_Dist
    MEC_map_sphere = [MEC_map, MEC_CEA_map_Dist];
    MEC_CEA_sphere = [MEC_CEA, MEC_CEA_map_Dist];
    MEC_area = clipPoints3d(distalFemurUSP.vertices, MEC_map_sphere, 'shape','sphere');
    MEC_area = [MEC_area; clipPoints3d(distalFemurUSP.vertices, MEC_CEA_sphere, 'shape','sphere')];
    [~, MEC_USP_Idx] = min(MEC_area(:,3));
    MEC_USP = MEC_area(MEC_USP_Idx,:);
else
    % Take MEC_max
    MEC_USP = MEC_max;
end

%% LEC
% Get the most lateral vertex of the distal femur in the in USP CS
[~, LEC_max_Idx] = max(distalFemurUSP.vertices(:,3));
LEC_max = distalFemurUSP.vertices(LEC_max_Idx,:);
% Lateral intersection of the CEA with distal femur
[~, LEC_CEA_Idx] = max(CEA_intersectPts(:,3));
LEC_CEA = CEA_intersectPts(LEC_CEA_Idx,:);
% Sanity check in case of osteophytes
LEC_CEA_map_Dist = distancePoints3d(LEC_map, LEC_CEA);
LEC_CEA_max_Dist = distancePoints3d(LEC_max, LEC_CEA);
if LEC_CEA_max_Dist > 20 % [mm]
    LEC_USP = LEC_map;
elseif LEC_CEA_max_Dist > LEC_CEA_map_Dist
    LEC_map_sphere = [LEC_map, LEC_CEA_map_Dist];
    LEC_CEA_sphere = [LEC_CEA, LEC_CEA_map_Dist];
    LEC_area = clipPoints3d(distalFemurUSP.vertices, LEC_map_sphere, 'shape','sphere');
    LEC_area = [LEC_area; clipPoints3d(distalFemurUSP.vertices, LEC_CEA_sphere, 'shape','sphere')];
    [~, LEC_USP_Idx] = max(LEC_area(:,3));
    LEC_USP = LEC_area(LEC_USP_Idx,:);
else
    % Take LEC_max
    LEC_USP = LEC_max;
end

if visu
    [~, axH, figH] = visualizeMeshes(distalFemurUSP);
    figH.NumberTitle = 'off';
    figH.Name = 'Epicondyle refinement (debug figure)';
    drawLine3d(axH, CEA)
    drawPoint3d(axH, [MEC_map; LEC_map], 'MarkerFaceColor','r','MarkerEdgeColor','r');
    drawLabels3d(axH, MEC_map, 'MEC_{map}', 'Color','r','HorizontalAlignment','Left')
    drawLabels3d(axH, LEC_map, 'LEC_{map}', 'Color','r','HorizontalAlignment','Right')
    drawPoint3d(axH, [MEC_max; LEC_max], 'MarkerFaceColor','g','MarkerEdgeColor','g');
    drawLabels3d(axH, MEC_max, 'MEC_{max}', 'Color','g','HorizontalAlignment','Left')
    drawLabels3d(axH, LEC_max, 'LEC_{max}', 'Color','g','HorizontalAlignment','Right')
    drawPoint3d(axH, [MEC_CEA; LEC_CEA], 'MarkerFaceColor','b','MarkerEdgeColor','b');
    drawLabels3d(axH, MEC_CEA, 'MEC_{CEA}', 'Color','b','HorizontalAlignment','Left')
    drawLabels3d(axH, LEC_CEA, 'LEC_{CEA}', 'Color','b','HorizontalAlignment','Right')
    if exist('MEC_area', 'var')
        drawPoint3d(axH, MEC_area,'MarkerEdgeColor','y');
    end
    if exist('LEC_area', 'var')
        drawPoint3d(axH, LEC_area,'MarkerEdgeColor','y');
    end
    drawPoint3d(axH, [MEC_USP; LEC_USP],'MarkerFaceColor','k','MarkerEdgeColor','k');
    drawLabels3d(axH, MEC_USP, 'MEC', 'Color','k','HorizontalAlignment','Left')
    drawLabels3d(axH, LEC_USP, 'LEC', 'Color','k','HorizontalAlignment','Right')
    
    anatomicalViewButtons(axH, 'ASR')
end

end