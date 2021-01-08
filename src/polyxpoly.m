function [xi, yi, ii] = polyxpoly(x1, y1, x2, y2, filterCode)
%POLYXPOLY  Intersection points for lines or polygon edges.
%
%   [XI,YI] = POLYXPOLY(X1,Y1,X2,Y2)  returns the intersection points of two polylines in a planar,
%             Cartesian system, with vertices defined by x1, y1, x2 and y2
%
%   [xi,yi,ii] = polyxpoly(...) returns a two-column array of line segment
%                indices corresponding to the intersection points
%
%   [xi,yi] = polyxpoly(...,'unique') filters out duplicate intersections,
%             which may result if the input polylines are self-intersecting.

%  Written by:  TinhNN
%  $Revision: 0.1 $    $Date: 2021/01/08

    % check input
    
    % check for valid filterCode
    
    % check x and y vectors
    
    
    % compute all intersection points
    xi = [];
    yi = [];
    ii = [];
    
    % format intersection points according to type and filterCode
    
end


