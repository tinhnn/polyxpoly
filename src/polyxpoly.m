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
    narginchk(4,5)
    
    % check for valid filterCode
    if nargin < 5
        filterCode = 'all';
    else
        filterCode = validatestring(filterCode, {'all','unique'});
    end
    
    % check x and y vectors
    checkxyvector(x1,y1)
    checkxyvector(x2,y2)
    
    % compute all intersection points
    [xi, yi, ii] = IntersecPointAll(x1,y1,x2,y2);
    
    % format intersection points according to type and filterCode
    if strcmp(filterCode,'unique') && (nargout < 3)
        [~,index] = unique([xi yi],'rows');
        index = sort(index);
        xi = xi(index);
        yi = yi(index);
    else
        if strcmp(filterCode,'unique')
            warning('''%s'' flag is ignored when third output (%s) is requested.', ...
                'unique', 'II')
        end
    end
end


function [xi,yi,ii] = IntersecPointAll(x1,y1,x2,y2)
%INTERSECPOINTALL  Unfiltered line or polygon intersection points.
%   [XI,YI,II] = INTERSECPOINTALL(X1,Y1,X2,Y2) returns the unfiltered intersection 
%   points of two sets of lines or polygons, along with a two-column index
%   of line segment numbers corresponding to the intersection points.
%   Note: intersection points are ordered from lowest to hightest line 
%   segment numbers.
    xi = [];
    yi = [];
    ii = [];
end

