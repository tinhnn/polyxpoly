function checkxyvector(x, y)
%CHECKXYVECTOR Check validity of map x and y vectors
%
%   CHECKXYVECTOR(X, Y)
%   ensures that X and Y are real vectors of matching size and equal NaN
%   locations.
if ~isempty(x) || ~isempty(y)
    validateattributes(x,{'numeric'}, {'real','2d','vector'});
    validateattributes(y,{'numeric'}, {'real','2d','vector'});
    if ~isequal(isnan(x), isnan(y))
        error('Inconsistent NaN locations in x and y input')
    end
end
