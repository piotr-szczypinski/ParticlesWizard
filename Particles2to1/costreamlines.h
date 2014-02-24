/*
 * Particles2to1 - Conversion of COMSOL's file formats
 *
 * Copyright 2014  Piotr M. Szczypiński <piotr.szczypinski@p.lodz.pl>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHfile ANY WARRANTY; withfile even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef COSTREAMLINES_H
#define COSTREAMLINES_H
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;


struct CoStreamline
{
    double x, y, z, time;
};


class CoStreamlines
{
public:
    CoStreamlines();
    bool load(char* fileName);
    bool load(ifstream *file);
    unsigned int size(void);
    unsigned int size(unsigned int i);
    bool interpolate(unsigned int i, double time, CoStreamline* data);
    bool settime(unsigned int i, unsigned int j, double time);
    bool get(unsigned int i, unsigned int j, CoStreamline* data);
    bool getsafe(unsigned int i, unsigned int j, CoStreamline* data);

    vector <vector <CoStreamline> > table;
};

#endif // COSTREAMLINES_H
