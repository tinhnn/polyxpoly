# intersections
## polyxpoly
Line or polygon intersection points.
## Syntax
```matlab
[xi,yi] = polyxpoly(x1,y1,x2,y2)
[xi,yi,ii] = polyxpoly(...)
[xi,yi] = polyxpoly(...,'unique')
```
## Description

`[xi,yi] = polyxpoly(x1,y1,x2,y2)` returns the intersection points of two polylines in a planar, Cartesian system, with vertices defined by x1, y1, x2 and y2. The output arguments, xi and yi, contain the x- and y-coordinates of each point at which a segment of the first polyline intersects a segment of the second. In the case of overlapping, collinear segments, the intersection is actually a line segment rather than a point, and both endpoints are included in xi, yi.

----

`[xi,yi,ii] = polyxpoly(___)` returns a two-column array of line segment indices corresponding to the intersection points. The k-th row of ii indicates which polyline segments give rise to the intersection point xi(k), yi(k).

----

`[xi,yi] = polyxpoly(___,'unique')` filters out duplicate intersections, which may result if the input polylines are self-intersecting.

----

## Examples

```matlab
x1 = [1 5 5 1 1];
y1 = [1 1 5 5 1];
x2 = [0 7];
y2 = [0 6];
plot(x1,y1);
hold on;
plot(x2,y2);
```

Intersect the straight segment with the rectangle.
```matlab
[xi, yi,ii] = polyxpoly(x1,y1,x2,y2) ;
plot(xi,yi,'*r');
```


## Reference
- https://www.mathworks.com/help/map/ref/polyxpoly.html
- https://github.com/hananel/piggottDesignCode/blob/master/code/polyxpoly.m
