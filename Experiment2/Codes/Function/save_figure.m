%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : VO_2 neuromorphic dynamics modulate fish social behavior through robotic embodiment
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function save_figure(fig, figureDir, fileBase)

pngPath = fullfile(figureDir, [fileBase '.png']);
pdfPath = fullfile(figureDir, [fileBase '.pdf']);
figPath = fullfile(figureDir, [fileBase '.fig']);

try
    exportgraphics(fig, pngPath, 'Resolution', 300);
catch
    saveas(fig, pngPath);
end
try
    exportgraphics(fig, pdfPath, 'ContentType', 'vector');
catch
    try
        saveas(fig, pdfPath);
    catch
    end
end
try
    savefig(fig, figPath);
catch
end

close(fig);

end
