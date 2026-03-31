% sets default xline properties
function h = xline_default(x, varargin)
    h = xline(x, varargin{:});
    h.FontSize = 12;             % change font size
    h.FontName = 'Times New Roman';  % change font family
end
