
#include "line_functions.h"

bool intersect(point p1, point p2, point p3, point p4) {
    float d1 = 0, d2 = 0, d3 = 0, d4 = 0;
    
    d1 = cross_product(p3, p4, p1);
    d2 = cross_product(p3, p4, p2);
    d3 = cross_product(p1, p2, p3);
    d4 = cross_product(p1, p2, p4);
    
    if(((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 >0)))
        return true;
    else if(d1 == 0 && on_line(p3, p4, p1))
        return true;
    else if(d2 == 0 && on_line(p3, p4, p2))
        return true;
    else if(d3 == 0 && on_line(p1, p2, p3))
        return true;
    else if(d4 == 0 && on_line(p1, p2, p4))
        return true;
    else
        return false;
}

float cross_product(point p0, point p1, point p2) {
    return (p1.x - p0.x)*(p2.y - p0.y) - (p2.x - p0.x)*(p1.y - p0.y);
}

bool on_line(point p0, point p1, point p2) {
    float min_x = 0, min_y = 0, max_x = 0, max_y = 0;
    
    if(p0.x < p1.x) {
        min_x = p0.x;
        max_x = p1.x;
    }
    else {
        min_x = p1.x;
        max_x = p0.x;
    }
    
    if(p0.y < p1.y) {
        min_y = p0.y;
        max_y = p1.y;
    }
    else {
        min_y = p1.y;
        max_y = p0.y;
    }
    
    if(min_x <= p2.x && p2.x <= max_x && min_y <= p2.y && p2.y <= max_y)
        return true;
    else
        return false;
}
