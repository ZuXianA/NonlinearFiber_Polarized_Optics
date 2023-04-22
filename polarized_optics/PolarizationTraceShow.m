function PolarizationTraceShow(dat, S)
    name = strcat(num2str(S(1)), ',', num2str(S(2)), ',', num2str(S(3)));
    filename = strcat('hwp', name);
    
    % 绘制Poincaré球上的轨迹
    fig = figure;
    ax = axes('Parent', fig, 'Projection', 'orthographic');
    plot3(ax, dat(:, 1), dat(:, 2), dat(:, 3), 'b-', 'LineWidth', 2, 'DisplayName', name);
    hold(ax, 'on');
    scatter3(ax, dat(1, 1), dat(1, 2), dat(1, 3), 'r', 'SizeData', 30, 'DisplayName', 'Start');
    scatter3(ax, dat(end, 1), dat(end, 2), dat(end, 3), 'g', 'SizeData', 30, 'DisplayName', 'End');
    xlim(ax, [-1, 1]);
    ylim(ax, [-1, 1]);
    zlim(ax, [-1, 1]);
    xlabel(ax, 'S1');
    ylabel(ax, 'S2');
    zlabel(ax, 'S3');
    legend(ax, 'Location', 'best');
    title(ax, filename);
    grid(ax, 'on');
    box(ax, 'on');
    pbaspect(ax, [1 1 1]);
    view(ax, [-38 30]);
    set(fig, 'Renderer', 'OpenGL');
    set(fig, 'WindowStyle', 'normal');
    set(fig, 'Units', 'pixels');
    set(fig, 'Position', [0 0 800 600]);
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'InvertHardcopy', 'off');
    set(fig, 'Color', [1 1 1]);
    set(fig, 'ToolBar', 'figure');
    set(fig, 'MenuBar', 'figure');
    drawnow;
    
    % 保存图像
    print(fig, filename, '-dpng', '-r300');
end
