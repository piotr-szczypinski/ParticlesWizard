/*
 * Copyright 2013  Piotr M. Szczypiński <piotr.szczypinski@p.lodz.pl>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
