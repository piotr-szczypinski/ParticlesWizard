#ifndef VVORONOI_H
#define VVORONOI_H


    #include "./libqhull/qhull_a.h"
    #include "./libqhull/libqhull.h"



typedef struct VoronoiData
{
    unsigned int verticesno;
    unsigned int cellsno;
    unsigned int validcellsno;
    double* vertices;
    int** cells;
    int* trajectoryindex;
    double circle[3];
} VoronoiData;

extern void collect_voronoi(VoronoiData* data);

#endif // VVORONOI_H
