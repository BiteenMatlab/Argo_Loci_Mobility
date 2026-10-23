function h = drawFullLine(ax,point,angle_degrees,varargin)
%drawFullLine Draw a line that spans the entire plot
%
%    drawFullLine(ax,point,angle_degrees) draws a line in the
%    specified axes that goes through the specified point at the
%    specified angle (in degrees). The line is drawn to span the
%    entire plot.
%
%    drawFullLine(___,Name,Value) passes name-value parameter pairs
%    to the line function.

% Steve Eddins


limits = axis(ax);
width = abs(limits(2) - limits(1));
height = abs(limits(4) - limits(3));
d = 2*hypot(width,height);
x1 = point(1) - d*cosd(angle_degrees);
x2 = point(1) + d*cosd(angle_degrees);
y1 = point(2) - d*sind(angle_degrees);
y2 = point(2) + d*sind(angle_degrees);
h = line(ax,'XData',[x1 x2],'YData',[y1 y2],varargin{:});
end

