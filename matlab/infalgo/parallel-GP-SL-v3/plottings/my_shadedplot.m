function [ha, hb, hc] = my_shadedplot(x, y1, y2, varargin)

% convert input vectors so that column vectors are also valid inputs
x = x(:)';
y1 = y1(:)';
y2 = y2(:)';
[ha, hb, hc] = shadedplot(x, y1, y2, varargin{:});

