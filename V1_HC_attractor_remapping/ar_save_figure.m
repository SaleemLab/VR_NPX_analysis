function ar_save_figure(fig, save_path, svg)
% .fig, .png and vector .pdf (and vector .svg for large scatters, svg = true).
% A file that is locked (e.g. open in a viewer) is written as <name>_new.<ext>.
if nargin < 3, svg = false; end
write = @(f, ext) save_one(fig, [save_path ext], f);
write(@(p) savefig(fig, p), '.fig');
write(@(p) exportgraphics(fig, p, 'Resolution', 150), '.png');
write(@(p) exportgraphics(fig, p, 'ContentType', 'vector'), '.pdf');
if svg, write(@(p) print(fig, p, '-dsvg', '-vector'), '.svg'); end
end

function save_one(~, path, f)
try
    f(path);
catch err
    [d, n, e] = fileparts(path);
    alt = fullfile(d, [n '_new' e]);
    warning('Could not write %s (%s); writing %s', path, err.message, alt);
    f(alt);
end
end
