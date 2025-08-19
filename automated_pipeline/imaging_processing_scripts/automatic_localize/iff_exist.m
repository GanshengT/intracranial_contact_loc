function y = iff_exist(x, dflt), if exist('x','var') && ~isempty(x), y=x; else, y=dflt; end, end
