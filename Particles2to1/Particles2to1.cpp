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

#include <fstream>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "getoption.h"
#include "costreamlines.h"
#include "covelocities.h"


#define epsilon 0.0000000001

//=============================================================================
// Function for scanning arguments

char* velocity_file = NULL;
char* input_file = NULL;
char* output_file = NULL;
double time_step = 0.001;
double time_end = -1.0;
double velo_add = 0.0;
bool isprinthelp = false;
bool debug_log = false;

double dist_mul = 1.0;
double time_mul = 1.0;

ofstream debug_stream;



int scan_parameters(int argc, char* argv[])
{
    int argi = 1;
    while(argi < argc)
    {
        GET_STRING_OPTION("-v", "--velocity", velocity_file)
        else GET_STRING_OPTION("-s", "--streamlines", input_file)
        else GET_STRING_OPTION("-o", "--output", output_file)
        else GET_DOUBLE_OPTION("-t", "--time", time_step)
        else GET_DOUBLE_OPTION("-x", "--veloadd", velo_add)
        else GET_DOUBLE_OPTION("-e", "--endtime", time_end)
        else GET_DOUBLE_OPTION("-m", "--distmul", dist_mul)
        else GET_DOUBLE_OPTION("-n", "--timemul", time_mul)
        else GET_NOARG_OPTION("-d", "--debug", debug_log, true)
        else GET_NOARG_OPTION("/?", "--help", isprinthelp, true)
        else return argi;
        argi++;
    }
    return 0;
}

//=============================================================================
// Print help
void printhelp(void)
{
//------------------------------------------------------------------------------------------v
    printf("Usage: Prticles2to1 [OPTION]...\n");
    printf("Converts COMSOL's streamline and velocity files particle flow file.\n");
    printf("Ignores units information, use -n and -m to correct output values.\n");
    printf("Version 2014.02.22, by Piotr M. Szczypinski\n");
    printf("Options:\n");
    printf("  -v, --velocity <file>     Velocity input file <file>\n");
    printf("  -s, --streamlines <file>  Streamline input file <file>\n");
    printf("  -o, --output <file>       Saves results to <file>\n");
    printf("  -t, --time <time>         Time step in units of input files (dafault 0.001)\n");
    printf("  -x, --veloadd <value>     Velocity to be added (dafault 0.0)\n");
    printf("  -e, --endtime <time>      End time (dafault autodetect)\n");
    printf("  -m, --distmul <value>     Output distance value multiplier (dafault 1.0)\n");
    printf("  -n, --timemul <value>     Output time value multiplier (dafault 1.0)\n");
    printf("  -d, --debug               Internal data dump to debug.log\n");
    printf("  /?, --help                Display this help and exit\n\n");


//    Example:
//    Particles2to1 -v velocity.txt -s streamlines.txt -o output.txt -t 20 -x 0.05 -d
//
}

//=============================================================================
bool savedata(char* filename, CoStreamlines* lines, double maxtime, double timestep, double distmul, double timemul)
{
    if(timestep <= epsilon || maxtime <= epsilon) return false;

    ofstream file;
    file.open(filename);
    if (!file.is_open()) return false;
    if (!file.good()) return false;

    unsigned int steps = maxtime / timestep;
    unsigned int linenu = lines->size();

    file << "% Nodes: " << linenu << "\n";
    file << "% Expressions: " << (steps+1)*3 << "\n";
    file << "% x ";

    for(unsigned int k = 0; k <= steps; k++)
    {
        double kk = k*timestep*timemul;
        file << " qx (?) @ t=" << kk << " qy (?) @ t=" << kk << " qz (?) @ t=" << kk;
    }
    file << "\n";

    for(unsigned int l = 0; l < linenu; l++)
    {
        file << l+1 << " ";

        if(lines->size(l) > 0)
        {
            unsigned int k;
            for(k = 0; k < steps+1; k++)
            {
                CoStreamline data;
                if(!lines->interpolate(l, k*timestep, &data))
                {
                    file << " NaN NaN NaN";
                    break;
                }
                file << " " << data.x*distmul;
                file << " " << data.y*distmul;
                file << " " << data.z*distmul;
            }
            for(; k < steps+1; k++) file << " NaN NaN NaN";
            file << "\n";
        }
    }
    file.close();
    return true;
}


//=============================================================================

bool computetimeshifts(CoStreamlines* lines, CoVelocities* velocities, double* maxtime)
{
    unsigned int dotpro_all = 0;
    unsigned int dotpro_small = 0;
    unsigned int zero_shift = 0;
    unsigned int linenu = lines->size();


    if(debug_stream.is_open())
    {
        debug_stream << "Stream lines = " << linenu << "\n";
        debug_stream << "Velocity added = " << velo_add << "\n";
    }


    for(unsigned int l = 0; l < linenu; l++)
    {
        if(lines->size(l) > 0)
        {
            lines->settime(l, 0, 0.0);
            CoStreamline prevelem;
            lines->settime(l, 0, 0.0);
            lines->get(l, 0, &prevelem);


            if(debug_stream.is_open())
            {
                debug_stream << l << " 0";
                debug_stream << " nodes(" << prevelem.x << "," << prevelem.y << "," << prevelem.z << ")";
                debug_stream << " time["<< prevelem.time << "]\n";
            }

            unsigned int linelen = lines->size(l);
            for(unsigned int k = 1; k < linelen; k++)
            {
                CoStreamline elem;
                unsigned int index = 0;
                lines->get(l, k, &elem);

                CoVelocity velo;
                velo.x = (elem.x + prevelem.x)/2.0;
                velo.y = (elem.y + prevelem.y)/2.0;
                velo.z = (elem.z + prevelem.z)/2.0;

                velocities->interpolate(&velo, &index);

                double dx, dy, dz;
                dx = elem.x - prevelem.x;
                dy = elem.y - prevelem.y;
                dz = elem.z - prevelem.z;

                double dd = (dx*dx+dy*dy+dz*dz);
                double dv = (velo.u*velo.u + velo.v*velo.v + velo.w*velo.w);
                dd = sqrt(dd);
                dv = sqrt(dv);

                if(dd < epsilon)
                {
                    zero_shift++;
                }

                if(dv < epsilon)
                {
                    cerr << "Error velocity zero: " << l << " " << k << " velocity index="<< index<<"\n";
                    return false;
                }

                double dotpro = (dx*velo.u+dy*velo.v+dz*velo.w);
                dotpro /= (dd*dv);
                dotpro_all++;
                if(dotpro < 0.94)
                {
                    dotpro_small ++;
                }
                dv *= dotpro;
                dd /= (dv+velo_add); //time
                dd += prevelem.time;

                lines->settime(l, k, dd);
                lines->get(l, k, &prevelem);
                if(debug_stream.is_open())
                {
                    debug_stream << l << " " << k << " dotprod["<< dotpro << "]";
                    debug_stream << " nodes(" << prevelem.x << "," << prevelem.y << "," << prevelem.z << ")->";
                    debug_stream << "(" << elem.x << "," << elem.y << "," << elem.z << ")";
                    debug_stream << " field(" << velo.x << "," << velo.y << "," << velo.z << ")";
                    debug_stream << "<" << velo.u << "," << velo.v << "," << velo.w << ">";
                    debug_stream << " time["<< dd << "]\n";
                }
                if(*maxtime < dd) *maxtime = dd;

                lines->get(l, k, &prevelem);
            }
        }
    }
    if(dotpro_small) cerr << "Warning: " << dotpro_small <<"/"<<dotpro_all<< " vectors do not align by more then 20 deg.\n";
    if(zero_shift) cerr << "Warning: " << zero_shift <<" duplicate nodes on stream lines\n";
    return true;
}

//=============================================================================
// Main function

int main(int argc, char* argv[])
{
    int ret = scan_parameters(argc, argv);

    if(isprinthelp || argc <= 1)
    {
        printhelp();
        return ret;
    }
    if(ret != 0)
    {
        if(ret < argc) fprintf(stderr, "Incorrect operand: %s\n", argv[ret]);
        fprintf(stderr, "Try LdaConsole --help for more information.\n");
        return ret;
    }
    if(debug_log) debug_stream.open("debug.log");

    if(velocity_file == NULL) {fprintf(stderr, "Unspecified velocities (-v)\n"); return -1;}
    if(input_file == NULL) {fprintf(stderr, "Unspecified streamlines (-i)\n"); return -2;}
    if(output_file == NULL) {fprintf(stderr, "Unspecified output file (-o)\n"); return -3;}
    if(time_step < 0.000000001) {fprintf(stderr, "Incorrect timestep (-t)\n"); return -4;}

    CoStreamlines lines;
    CoVelocities velocities;
    double maxtime = -1.0;
    if(!lines.load(input_file)) {fprintf(stderr, "Cannot load sreamlines\n"); return -5;}
    if(!velocities.load(velocity_file)) {fprintf(stderr, "Cannot load velocities\n"); return -6;}
    if(debug_stream.is_open())
    {
        CoVelocity v;
        unsigned int imax = velocities.size();
        for(unsigned int i = 0; i < imax; i++)
        {
            velocities.get(i, &v);
            debug_stream << i << " (" << v.x << "," << v.y << "," << v.z <<")<" << v.u << "," << v.v << "," << v.w <<">\n";
        }
    }

    if(velo_add < 0.0) velo_add = 0.0;
    unsigned int allvel = velocities.size();
    unsigned int unzeroed = velocities.unzero(epsilon);
    if(unzeroed > 0) fprintf(stderr, "Warning: Ignores %i/%i zero velocity vectors\n", unzeroed, allvel);

    if(!computetimeshifts(&lines, &velocities, &maxtime)) {fprintf(stderr, "Cannot compute timeshifts\n"); return -8;}
    if(time_end > 0.0)
    {
        int ss = int(time_end/time_step)+1;
        cout << "Maximum time detected = " << maxtime << " maximum time used = " << ss*time_step << " time step = " << time_step << " steps = " << ss << "\n";
        maxtime = time_end;
    }
    else
    {
        int ss = int(maxtime/time_step)+1;
        cout << "Maximum time detected = " << maxtime << " maximum time used = " << ss*time_step << " time step = " << time_step << " steps = " << ss << "\n";
    }
    if(!savedata(output_file, &lines, maxtime, time_step, dist_mul, time_mul)) {fprintf(stderr, "Cannot save output file\n"); return -9;}

    if(debug_log) debug_stream.close();
    return 0;
}
