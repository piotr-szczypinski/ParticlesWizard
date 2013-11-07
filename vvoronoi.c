
#include "vvoronoi.h"


void getfacecenter(facetT *facet, double vert[2]) {
    int k, num;

    if (qh CENTERtype != qh_ASvoronoi && qh CENTERtype != qh_AScentrum)
        return;

    if (qh CENTERtype == qh_ASvoronoi) {
        num= qh hull_dim-1;
        if(num > 2) num = 2;
        if (!facet->normal || !facet->upperdelaunay || !qh ATinfinity) {
            if (!facet->center)
                facet->center= qh_facetcenter(facet->vertices);
            for (k=0; k < num; k++)
                vert[k] = facet->center[k];
        }else {
            for (k=0; k < num; k++)
                vert[k] = qh_INFINITE;

        }
    }
}

void collect_voronoi(VoronoiData *data)
{
    int ver;
    double verd[2];

    int k, numcenters, numvertices= 0, numneighbors, numinf, vid=1, vertex_i, vertex_n;
    facetT *facet, **facetp, *neighbor, **neighborp;
    setT *vertices;
    vertexT *vertex;
    boolT isLower;

    //facetT *facetlist;
    unsigned int numfacets= (unsigned int) qh num_facets;
    //facetlist = qh facet_list;

    vertices= qh_markvoronoi(qh facet_list, NULL, !qh_ALL, &isLower, &numcenters);
    FOREACHvertex_i_(vertices) {
        if (vertex) {
            numvertices++;
            numneighbors = numinf = 0;
            FOREACHneighbor_(vertex) {
                if (neighbor->visitid == 0)
                    numinf= 1;
                else if (neighbor->visitid < numfacets)
                    numneighbors++;
            }
            if (numinf && !numneighbors) {
                SETelem_(vertices, vertex_i)= NULL;
                numvertices--;
            }
        }
    }

    if(data->trajectoryindex != NULL) free(data->trajectoryindex);
    if(data->vertices != NULL) free(data->vertices);
    if(data->cells != NULL)
    {
        for(k = 0; k < data->cellsno; k++) free(data->cells[k]);
        free(data->cells);
    }

    data->verticesno = numcenters;
    data->trajectoryindex = (int*) malloc(sizeof(int) * numcenters);
    data->vertices = (double*) malloc(sizeof(double) * 2 * numcenters);
    data->cellsno = qh_setsize(vertices);
    data->cells = (int**) malloc(sizeof(int*) * qh_setsize(vertices));
    for(k = 0; k < data->cellsno; k++) data->cells[k] = NULL;

    ver = 0;
    for (k=qh hull_dim-1; k--; )
    {
        data->vertices[ver] = qh_INFINITE;
        ver++;
    }
    FORALLfacet_(qh facet_list) {
        if (facet->visitid && facet->visitid < numfacets)
        {
            getfacecenter(facet, verd);
            data->vertices[ver] = verd[0]; ver++;
            data->vertices[ver] = verd[1]; ver++;
        }
    }


    ver = 0;

    FOREACHvertex_i_(vertices)
    {
        numneighbors= 0;
        numinf=0;
        if (vertex) {
            if (qh hull_dim == 3)
                qh_order_vertexneighbors(vertex);
            else if (qh hull_dim >= 4)
                qsort(SETaddr_(vertex->neighbors, facetT),
                      (size_t)qh_setsize(vertex->neighbors),
                      sizeof(facetT *), qh_compare_facetvisit);
            FOREACHneighbor_(vertex) {
                if (neighbor->visitid == 0)
                    numinf= 1;
                else if (neighbor->visitid < numfacets)
                    numneighbors++;
            }
        }

        if (numinf)
            numneighbors++;

        data->cells[ver] = (int*) malloc(sizeof(int) * (numneighbors + 1));
        data->cells[ver][0] = numneighbors;
        k = 1;

        if (vertex) {
            FOREACHneighbor_(vertex) {
                if (neighbor->visitid == 0) {
                    if (numinf)
                    {
                        numinf= 0;
                        data->cells[ver][k] = neighbor->visitid;
                        k++;
                    }
                }
                else if (neighbor->visitid < numfacets)
                {
                    data->cells[ver][k] = neighbor->visitid;
                    k++;
                }
            }
        }
        ver++;
    }
    qh_settempfree(&vertices);
}
