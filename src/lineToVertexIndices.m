function Idx = lineToVertexIndices(line, mesh)
%LINETOVERTEXINDICES converts the representation of a line from origin and
% direction vector to two vertex indices of a the mesh.
%
% AUTHOR: Maximilian C. M. Fischer
% COPYRIGHT (C) 2020-2023 Maximilian C. M. Fischer
% LICENSE: EUPL v1.2
%

[linePts, linePos] = intersectLineMesh3d(line, mesh.vertices, mesh.faces);
[linePts, linePtsIdx] = unique(linePts,'rows','stable');
[~, linePosIdx] = sort(linePos(linePtsIdx));
linePts = linePts(linePosIdx,:);
[~, Idx] = pdist2(mesh.vertices,[linePts(1,:);linePts(end,:)],'euclidean','Smallest',1);

end