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

#include "covelocities.h"

CoVelocities::CoVelocities()
{
}

bool CoVelocities::load(char* fileName)
{
    ifstream file;
    file.open(fileName);
    if (!file.is_open()) return false;
    if (!file.good()) return false;
    bool ret = load(&file);
    file.close();
    return ret;
}

bool CoVelocities::load(ifstream *file)
{
    table.clear();
    string inputstring;
    do
    {
        getline(*file, inputstring);
    }
    while(!file->eof() && inputstring.at(0) == '%');

    do
    {
        CoVelocity d;
        *file >> d.x;
        *file >> d.y;
        *file >> d.z;
        *file >> d.u;
        *file >> d.v;
        *file >> d.w;

        if(!file->eof())
        {
            table.push_back(d);
        }
    }
    while(!file->eof());

    if(table.size()>0) return true;
    return false;
}

unsigned int CoVelocities::size(void)
{
    return table.size();
}


bool CoVelocities::get(unsigned int i, CoVelocity* data)
{
    if(i >= table.size()) return false;
    *data = table[i];
    return true;
}

unsigned int CoVelocities::unzero(double epsilon)
{
    unsigned int er = 0;
    for(unsigned int i = table.size()-1; i < table.size(); i--)
    {
        CoVelocity* vv = &table[i];
        double dd = vv->u*vv->u+vv->v*vv->v+vv->w*vv->w;
        if(dd <= epsilon)
        {
            er++;
            table.erase(table.begin()+i);
        }
    }
    return er;
}

bool CoVelocities::interpolate(CoVelocity* data, unsigned int* index)
{
    double mind = -1.0;
    double xx = data->x;
    double yy = data->y;
    double zz = data->z;
    for(vector<CoVelocity>::iterator i = table.begin(); i != table.end(); ++i)
    {
        double x = i->x - data->x;
        double y = i->y - data->y;
        double z = i->z - data->z;
        double dd = x*x+y*y+z*z;

        if(dd < mind || mind < 0)
        {
            mind = dd;
            data->u = i->u;
            data->v = i->v;
            data->w = i->w;
            xx = i->x;
            yy = i->y;
            zz = i->z;

            *index = i-table.begin();
        }
    }
    if(mind >= 0.0)
    {
        data->x = xx;
        data->y = yy;
        data->z = zz;
        return true;
    }
    return false;
}
