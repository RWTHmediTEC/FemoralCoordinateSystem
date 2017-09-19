function varargout = visualizeMeshes(Mesh, PatchProps)

if nargin == 1
    PatchProps(1).EdgeColor = 'none';
    PatchProps(1).FaceColor = [223, 206, 161]/255;
    PatchProps(1).FaceAlpha = 1;
    PatchProps(1).FaceLighting = 'gouraud';
end

% New figure
MonitorsPos = get(0,'MonitorPositions');
figHandle = figure('Units','pixels','renderer','opengl', 'Color', 'w');
% figHandle.ToolBar='none';
% figHandle.MenuBar='none';
% figHandle.WindowScrollWheelFcn=@M_CB_Zoom;
% FigHandle.WindowButtonDownFcn=@M_CB_RotateWithLeftMouse;
if     size(MonitorsPos,1) == 1
    set(figHandle,'OuterPosition',[1 50 MonitorsPos(1,3)-1 MonitorsPos(1,4)-50]);
elseif size(MonitorsPos,1) == 2
    set(figHandle,'OuterPosition',[1+MonitorsPos(1,3) 50 MonitorsPos(2,3)-1 MonitorsPos(2,4)-50]);
end

meshHandle=zeros(length(Mesh),1);
for i=1:length(Mesh)
    meshHandle(i) = patch(Mesh(i), PatchProps);
end

H_Light(1) = light; light('Position', -1*(get(H_Light(1),'Position')));
% cameratoolbar('SetCoordSys','none')
axis on; axis equal; hold on
xlabel x; ylabel y; zlabel z;

if nargout==1
    varargout = {meshHandle};
elseif nargout==2
    varargout{1} = {meshHandle};
    varargout{2} = {figHandle};
end

end