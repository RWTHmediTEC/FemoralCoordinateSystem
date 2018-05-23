function viewButtonsRAS(varargin)

switch length(varargin)
    case 0
        hAx = gca;
        mode = 'RAS';
    case 1
        if numel(varargin{1}) == 1 && ishandle(varargin{1})
            hAx = varargin{1};
        else
            hAx = gca;
            mode = varargin{1};
        end
    case 2
        hAx = varargin{1};
        mode = varargin{2};
end

mouseControl3d(hAx)

% uicontrol Button Size
BSX = 0.1; BSY = 0.023;

%Font properies
FontPropsA.FontUnits = 'normalized';
FontPropsA.FontSize = 0.8;
% Rotate-buttons
switch mode
    case 'RAS'
        uicontrol('Units','normalized','Position',[0.5-BSX*3/2 0.01+BSY BSX BSY],FontPropsA,...
            'String','Left','Callback','mouseControl3d(gca, [0 1 0 0; -1 0 0 0; 0 0 1 0; 0 0 0 1])');
        uicontrol('Units','normalized','Position',[0.5-BSX*3/2     0.01 BSX BSY],FontPropsA,...
            'String','Right','Callback','mouseControl3d(gca, [0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1])');
        uicontrol('Units','normalized','Position',[0.5-BSX*1/2 0.01+BSY BSX BSY],FontPropsA,...
            'String','Anterior','Callback','mouseControl3d(gca, [1 0 0 0;0 1 0 0; 0 0 1 0; 0 0 0 1])');
        uicontrol('Units','normalized','Position',[0.5-BSX*1/2     0.01 BSX BSY],FontPropsA,...
            'String','Posterior','Callback','mouseControl3d(gca, [-1 0 0 0;0 -1 0 0; 0 0 1 0; 0 0 0 1])');
        uicontrol('Units','normalized','Position',[0.5+BSX*1/2 0.01+BSY BSX BSY],FontPropsA,...
            'String','Superior','Callback','mouseControl3d(gca, [1 0 0 0;0 0 1 0; 0 -1 0 0; 0 0 0 1])');
        uicontrol('Units','normalized','Position',[0.5+BSX*1/2     0.01 BSX BSY],FontPropsA,...
            'String','Inferior','Callback','mouseControl3d(gca, [-1 0 0 0;0 0 -1 0; 0 -1 0 0; 0 0 0 1])');
    case 'ASR'
        uicontrol('Units','normalized','Position',[0.5-BSX*3/2     0.01 BSX BSY],FontPropsA,...
            'String','Posterior','Callback','mouseControl3d(gca, [0 0 -1 0; -1 0 0 0; 0 1 0 0; 0 0 0 1])');
        uicontrol('Units','normalized','Position',[0.5-BSX*3/2 0.01+BSY BSX BSY],FontPropsA,...
            'String','Anterior','Callback','mouseControl3d(gca, [0 0 1 0; 1 0 0 0; 0 1 0 0; 0 0 0 1])');
        uicontrol('Units','normalized','Position',[0.5-BSX*1/2 0.01+BSY BSX BSY],FontPropsA,...
            'String','Superior','Callback','mouseControl3d(gca, [0 0 1 0;0 1 0 0; -1 0 0 0; 0 0 0 1])');
        uicontrol('Units','normalized','Position',[0.5-BSX*1/2     0.01 BSX BSY],FontPropsA,...
            'String','Inferior','Callback','mouseControl3d(gca, [0 0 1 0;0 -1 0 0; 1 0 0 0; 0 0 0 1])');
        uicontrol('Units','normalized','Position',[0.5+BSX*1/2 0.01+BSY BSX BSY],FontPropsA,...
            'String','Left','Callback','mouseControl3d(gca, [1 0 0 0;0 0 -1 0; 0 1 0 0; 0 0 0 1])');
        uicontrol('Units','normalized','Position',[0.5+BSX*1/2     0.01 BSX BSY],FontPropsA,...
            'String','Right','Callback','mouseControl3d(gca, [-1 0 0 0;0 0 1 0; 0 1 0 0; 0 0 0 1])');
end

end