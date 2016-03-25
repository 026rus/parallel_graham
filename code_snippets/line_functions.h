


#ifndef LINE_FUNCTIONS_H
#define LINE_FUNCTIONS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

//Define point structure
typedef struct {
    float x, y;
} point;

float cross_product(point p0, point p1, point p2);
bool on_line(point p0, point p1, point p2);
bool intersect(point p1, point p2, point p3, point p4);

#ifdef __cplusplus
}
#endif

#endif /* LINE_FUNCTIONS_H */

