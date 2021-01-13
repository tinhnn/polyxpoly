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
    
    % convert to column vectors
    x1 = x1(:);
    y1 = y1(:);
    x2 = x2(:);
    y2 = y2(:);
    
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
        % Ensure that intersections are reported only once.
        if numel(xi) > 1
            [xi, yi, ii] = filterIntersections(xi,yi,ii);
        end
        if strcmp(filterCode,'unique')
            warning('''%s'' flag is ignored when third output (%s) is requested.', ...
                'unique', 'II')
        end
    end
end

%--------------------------------------------------------------------------
function [xi,yi,ii] = IntersecPointAll(x1,y1,x2,y2)
%INTERSECPOINTALL  Unfiltered line or polygon intersection points.
%   [XI,YI,II] = INTERSECPOINTALL(X1,Y1,X2,Y2) returns the unfiltered intersection 
%   points of two sets of lines or polygons, along with a two-column index
%   of line segment numbers corresponding to the intersection points.
%   Note: intersection points are ordered from lowest to hightest line 
%   segment numbers.
    err = eps*1e5;

    % form line segment matrices
    xs1 = [x1 [x1(2:end); x1(1)]];
    ys1 = [y1 [y1(2:end); y1(1)]];
    xs2 = [x2 [x2(2:end); x2(1)]];
    ys2 = [y2 [y2(2:end); y2(1)]];

    % remove last segment (for self-enclosed polygons, this is a non-segment;
    % for lines, there are n-1 line segments)
    xs1 = xs1(1:end-1,:);  ys1 = ys1(1:end-1,:);
    xs2 = xs2(1:end-1,:);  ys2 = ys2(1:end-1,:);

    % tile matrices for vectorized intersection calculations
    N1 = length(xs1(:,1));  N2 = length(xs2(:,1));
    X1 = reshape(repmat(xs1,1,N2)',2,N1*N2)';
    Y1 = reshape(repmat(ys1,1,N2)',2,N1*N2)';
    X2 = repmat(xs2,N1,1);
    Y2 = repmat(ys2,N1,1);

    % compute slopes
    w = warning;
    warning off
    m1 = (Y1(:,2) - Y1(:,1)) ./ (X1(:,2) - X1(:,1));
    m2 = (Y2(:,2) - Y2(:,1)) ./ (X2(:,2) - X2(:,1));
    m1(find(abs(m1)>1/err)) = inf;  m2(find(abs(m2)>1/err)) = inf;
    warning(w)

    % compute y-intercepts (note: imaginary values for vertical lines)
    b1 = zeros(size(m1));  b2 = zeros(size(m2));
    i1 = find(m1==inf);  if ~isempty(i1),  b1(i1) = X1(i1)*i;  end
    i2 = find(m2==inf);  if ~isempty(i2),  b2(i2) = X2(i2)*i;  end
    i1 = find(m1~=inf);  if ~isempty(i1),  b1(i1) = Y1(i1) - m1(i1).*X1(i1);  end
    i2 = find(m2~=inf);  if ~isempty(i2),  b2(i2) = Y2(i2) - m2(i2).*X2(i2);  end

    % zero intersection coordinate arrays
    sz = size(X1(:,1));  x0 = zeros(sz);  y0 = zeros(sz);

    % parallel lines (do not intersect except for similar lines)
    % for similar lines, take the low and high points
    idx = find( abs(m1-m2)<err | (isinf(m1)&isinf(m2)) );
    if ~isempty(idx)
    % non-similar lines
        sub = find(abs(b1(idx)-b2(idx))>err);  j = idx(sub);
        x0(j) = nan;  y0(j) = nan;
    % similar lines (non-vertical)
        sub = find(abs(b1(idx)-b2(idx))<err & m1(idx)~=inf);  j = idx(sub);
        Xlo = max([min(X1(j,:),[],2) min(X2(j,:),[],2)],[],2);
        Xhi = min([max(X1(j,:),[],2) max(X2(j,:),[],2)],[],2);
        if ~isempty(j)
            j0 = find(abs(Xlo-Xhi)<=err);
            j1 = find(abs(Xlo-Xhi)>err);
            x0(j(j0)) = Xlo(j0);
            y0(j(j0)) = Y1(j(j0)) + m1(j(j0)).*(Xlo(j0) - X1(j(j0)));
            x0(j(j1)) = Xlo(j1) + i*Xhi(j1);
            y0(j(j1)) = (Y1(j(j1)) + m1(j(j1)).*(Xlo(j1) - X1(j(j1)))) + ...
                         i*(Y1(j(j1)) + m1(j(j1)).*(Xhi(j1) - X1(j(j1))));
        end
    % similar lines (vertical)
        sub = find(abs(b1(idx)-b2(idx))<err & m1(idx)==inf);  j = idx(sub);
        Ylo = max([min(Y1(j,:),[],2) min(Y2(j,:),[],2)],[],2);
        Yhi = min([max(Y1(j,:),[],2) max(Y2(j,:),[],2)],[],2);
        if ~isempty(j)
            y0(j) = Ylo + i*Yhi;
            x0(j) = X1(j) + i*X1(j);
        end
    end

    % non-parallel lines
    idx = find(abs(m1-m2)>err);
    if ~isempty(idx)
    % non-vertical/non-horizontal lines
        sub = find(m1(idx)~=inf & m2(idx)~=inf & ...
                   abs(m1(idx))>eps & abs(m2(idx))>eps);
        j = idx(sub);
        x0(j) = (Y1(j) - Y2(j) + m2(j).*X2(j) - m1(j).*X1(j)) ./ ...
                (m2(j) - m1(j));
        y0(j) = Y1(j) + m1(j).*(x0(j)-X1(j));
    % first line vertical
        sub = find(m1(idx)==inf);  j = idx(sub);
        x0(j) = X1(j);
        y0(j) = Y2(j) + m2(j).*(x0(j)-X2(j));
    % second line vertical
        sub = find(m2(idx)==inf);  j = idx(sub);
        x0(j) = X2(j);
        y0(j) = Y1(j) + m1(j).*(x0(j)-X1(j));
    % first line horizontal, second line non-vertical
        sub = find(abs(m1(idx))<=eps & m2(idx)~=inf);  j = idx(sub);
        y0(j) = Y1(j);
        x0(j) = (Y1(j) - Y2(j) + m2(j).*X2(j)) ./ m2(j);
    % second line horizontal, first line non-vertical
        sub = find(abs(m2(idx))<=eps & m1(idx)~=inf);  j = idx(sub);
        y0(j) = Y2(j);
        x0(j) = (Y1(j) - y0(j) - m1(j).*X1(j)) ./ -m1(j);
    end

    % throw away points that lie outside of line segments
    dx1 = [min(X1,[],2)-x0, x0-max(X1,[],2)];
    dy1 = [min(Y1,[],2)-y0, y0-max(Y1,[],2)];
    dx2 = [min(X2,[],2)-x0, x0-max(X2,[],2)];
    dy2 = [min(Y2,[],2)-y0, y0-max(Y2,[],2)];
    [irow,icol] = find([dx1 dy1 dx2 dy2]>err);
    idx = sort(unique(irow));
    x0(idx) = nan;
    y0(idx) = nan;

    % retrieve only intersection points (no nans)
    idx = find(~isnan(x0));
    xi = x0(idx);  yi = y0(idx);

    % determine indices of line segments that intersect
    i1 = ceil(idx/N2);  i2 = rem(idx,N2);
    if ~isempty(i2),  i2(find(i2==0)) = N2;  end
    ii = [i1 i2];

    % combine all intersection points
    indx = union(find(imag(xi)),find(imag(yi)));
    for n=length(indx):-1:1
        j = indx(n);
        ii = [ii(1:j-1,:); ii(j,:); ii(j:end,:)];
        xi = [xi(1:j-1); imag(xi(j)); real(xi(j:end))];
        yi = [yi(1:j-1); imag(yi(j)); real(yi(j:end))];
    end
end

%--------------------------------------------------------------------------
function [xi,yi,ii] = filterIntersections(xi,yi,ii)
    % TODO
end