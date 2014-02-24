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

#include "costreamlines.h"

CoStreamlines::CoStreamlines()
{
}
bool CoStreamlines::load(char* fileName)
{
    ifstream file;
    file.open(fileName);
    if (!file.is_open()) return false;
    if (!file.good()) return false;
    bool ret = load(&file);
    file.close();
    return ret;
}

bool CoStreamlines::load(ifstream *file)
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
        CoStreamline d;
        unsigned int i;
        double c;
        *file >> d.x;
        *file >> d.y;
        *file >> d.z;
        *file >> i;
        *file >> c;

        d.time = -1;
        if(!file->eof())
        {
            if(i >= table.size()) table.resize(i+1);
            table[i].push_back(d);
        }
    }
    while(!file->eof());

    if(table.size()>0) return true;
    return false;
}
bool CoStreamlines::settime(unsigned int i, unsigned int j, double time)
{
    if(table.size() <= i) return false;
    if(table[i].size() <= j) return false;
    table[i][j].time = time;
    return true;
}
bool CoStreamlines::getsafe(unsigned int i, unsigned int j, CoStreamline* data)
{
    if(table.size() <= i) return false;
    if(table[i].size() <= j) return false;
    *data = table[i][j];
    return true;
}

bool CoStreamlines::get(unsigned int i, unsigned int j, CoStreamline* data)
{
    *data = table[i][j];
    return true;
}


unsigned int CoStreamlines::size(unsigned int i)
{
    if(table.size() <= i) return 0;
    return table[i].size();
}
unsigned int CoStreamlines::size(void)
{
    return table.size();
}

bool CoStreamlines::interpolate(unsigned int i, double time, CoStreamline* data)
{
    static unsigned int lasti = -1;
    static unsigned int lasta = -1;
    unsigned int lastb = -1;
    if(table.size() <= i) return false;

    vector <CoStreamline> line = table[i];
    vector<CoStreamline>::iterator linestart;
    if(lasti == i)
    {
        if(lasta >= table[i].size()) return false;
        if(lasta == 0) linestart = line.begin();
        else linestart = line.begin()+lasta-1;
    }
    else linestart = line.begin();

    lasta = -1;
    for(vector<CoStreamline>::iterator ii = linestart; ii != line.end(); ++ii)
    {
        if(ii->time >= time)
        {
            lasta = ii-line.begin();
            break;
        }
    }

    if(lasta >= line.size()) return false;

    if(lasta > 0) lastb = lasta - 1;
    else
    {
        if(time == line[lasta].time)
        {
            *data = line[lasta];
            return true;
        }
        return false;
    }

    double fracb = (time - line[lastb].time) / (line[lasta].time - line[lastb].time);
    double fraca = (1.0 - fracb);

    data->time = time;
    data->x = line[lastb].x * fraca + line[lasta].x * fracb;
    data->y = line[lastb].y * fraca + line[lasta].y * fracb;
    data->z = line[lastb].z * fraca + line[lasta].z * fracb;
    return true;
}
