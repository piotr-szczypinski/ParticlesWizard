/*
 * Particles helper functions
 * 
 * Copyright 2013  Piotr M. Szczypiński <piotr.szczypinski@p.lodz.pl>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either 
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public 
 * License along with this library.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>

unsigned int number_of_trajectories = 0;
double node_repetition_time = -1.0;
double** trajectory_nodes = NULL;
double* randomization_time = NULL;
int* last_in_tube_particle_index_at_start = NULL; // = (unsigned int)((nodes_on_trajectory[trajectory]*node_repetition_time - randomization_time[trajectory]) / particle_repetition_time[trajectory]);

unsigned int* nodes_on_trajectory = NULL;
double* particle_repetition_time = NULL;
double* particle_volume = NULL;
double minmaxxyz[3][2];


// Loads file with tesselation data and trajectories. Returns 1 on success.
int czLoadTrajectories(const char* filename)
{
    int e, i;
    float f;
    int first = 1;
    FILE *plik = fopen(filename, "rt");
    if(plik == NULL) return 0;

    if(randomization_time != NULL) free(randomization_time); randomization_time = NULL;
    if(particle_repetition_time != NULL) free(particle_repetition_time); particle_repetition_time = NULL;
    if(particle_volume != NULL) free(particle_volume); particle_volume = NULL;
    if(nodes_on_trajectory != NULL) free(nodes_on_trajectory); nodes_on_trajectory = NULL;
    if(last_in_tube_particle_index_at_start != NULL) free(last_in_tube_particle_index_at_start); last_in_tube_particle_index_at_start = NULL;
    if(trajectory_nodes != NULL)
    {
        for(unsigned int trajectory = 0; trajectory < number_of_trajectories; trajectory++) free(trajectory_nodes[trajectory]);
        free(trajectory_nodes);
        trajectory_nodes = NULL;
    }

    e = fscanf(plik, "number_of_trajectories %i ", &i);
    if(e <= 0) {fclose(plik); return 0;}
    number_of_trajectories = i;
    e = fscanf(plik, "node_repetition_time %e ", &f);
    if(e <= 0) {fclose(plik); return 0;}
    node_repetition_time = f;

    randomization_time = (double*) malloc(sizeof(double) * number_of_trajectories);
    if(randomization_time == NULL) {fclose(plik); return 0;}

    particle_repetition_time = (double*) malloc(sizeof(double) * number_of_trajectories);
    if(particle_repetition_time == NULL) {fclose(plik); return 0;}

    particle_volume = (double*) malloc(sizeof(double) * number_of_trajectories);
    if(particle_volume == NULL) {fclose(plik); return 0;}

    nodes_on_trajectory = (unsigned int*) malloc(sizeof(unsigned int) * number_of_trajectories);
    if(nodes_on_trajectory == NULL) {fclose(plik); return 0;}

    last_in_tube_particle_index_at_start = (int*) malloc(sizeof(int) * number_of_trajectories);
    if(last_in_tube_particle_index_at_start == NULL) {fclose(plik); return 0;}

    trajectory_nodes = (double**) malloc(sizeof(double*) * number_of_trajectories);
    if(trajectory_nodes == NULL) {fclose(plik); return 0;}
    for(unsigned int trajectory = 0; trajectory < number_of_trajectories; trajectory++)
        trajectory_nodes[trajectory] = NULL;

    for(unsigned int trajectory = 0; trajectory < number_of_trajectories; trajectory++)
    {
        e = fscanf(plik, "trajectory %i ", &i); if(e <= 0) {fclose(plik); return 0;}
        if((int)trajectory != i) {fclose(plik); return 0;}

        e = fscanf(plik, "particle_volume %e ", &f); if(e <= 0) {fclose(plik); return 0;}
        particle_volume[trajectory] = f;
        e = fscanf(plik, "particle_repetition_time %e ", &f); if(e <= 0) {fclose(plik); return 0;}
        particle_repetition_time[trajectory] = f;
        e = fscanf(plik, "randomization_time %e ", &f); if(e <= 0) {fclose(plik); return 0;}
        randomization_time[trajectory] = f;

        e = fscanf(plik, "trajectory_nodes %i ", &i); if(e <= 0) {fclose(plik); return 0;}
        if(i < 2) {fclose(plik); return 0;}
        nodes_on_trajectory[trajectory] = i;
        trajectory_nodes[trajectory] = (double*) malloc(sizeof(double) * nodes_on_trajectory[trajectory] * 3);

        last_in_tube_particle_index_at_start[trajectory]
        = (int)((nodes_on_trajectory[trajectory]*node_repetition_time - randomization_time[trajectory]) / particle_repetition_time[trajectory]);

        for(i = 0; i < (int)nodes_on_trajectory[trajectory]; i ++)
        {
            float x, y, z;
            e = fscanf(plik, "(%e %e %e) ", &x, &y, &z); if(e < 3) {fclose(plik); return 0;}
            trajectory_nodes[trajectory][i*3] = x;
            trajectory_nodes[trajectory][i*3+1] = y;
            trajectory_nodes[trajectory][i*3+2] = z;
            if(first)
            {
                first = 0;
                minmaxxyz[0][0] = minmaxxyz[0][1] = x;
                minmaxxyz[1][0] = minmaxxyz[1][1] = y;
                minmaxxyz[2][0] = minmaxxyz[2][1] = z;
            }
            else
            {
                if(minmaxxyz[0][0] > x) minmaxxyz[0][0] = x;
                if(minmaxxyz[0][1] < x) minmaxxyz[0][1] = x;
                if(minmaxxyz[1][0] > y) minmaxxyz[1][0] = y;
                if(minmaxxyz[1][1] < y) minmaxxyz[1][1] = y;
                if(minmaxxyz[2][0] > z) minmaxxyz[2][0] = z;
                if(minmaxxyz[2][1] < z) minmaxxyz[2][1] = z;
            }
        }
    }
    return 1;
}

// Writes min and max of coordinates. Returns 1 on success.
void GetMinMaxXYZ(double* xb, double* xt, double* yb, double* yt, double* zb, double* zt)
{
        *xb = minmaxxyz[0][0];
        *xt = minmaxxyz[0][1];
        *yb = minmaxxyz[1][0];
        *yt = minmaxxyz[1][1];
        *zb = minmaxxyz[2][0];
        *zt = minmaxxyz[2][1];
}



// Returns number of trajectories.
unsigned int czNumberOfTrajectories(void)
{
    return number_of_trajectories;
}

// Writes min and max with minimum and maximum index of in-tube particle.
void czGetMinMaxInTube(double current_time, unsigned int trajectory, int* min, int* max)
{
    *max = last_in_tube_particle_index_at_start[trajectory] + (current_time + randomization_time[trajectory]) / particle_repetition_time[trajectory];
    *min = last_in_tube_particle_index_at_start[trajectory] +
            (int)(
                (
                 current_time + randomization_time[trajectory] - (nodes_on_trajectory[trajectory])*node_repetition_time
                ) / particle_repetition_time[trajectory]
               );

   if(*min < 0) *min = 0;
   if(*max < 0) *min = 0;
}

// Gets volume of a particle on a given trajectory. Returns 1 on success.
int czGetVolume(unsigned int trajectory, double* v)
{
    if(particle_volume == NULL || trajectory >= number_of_trajectories) return 0;
    *v = particle_volume[trajectory];
    return 1;
}

// Computes coordinates of a particle on a given trajectory
// after a given time (deltat) from its entrance
int InterpolateCoordinates(double* x, double* y, double* z, unsigned int trajectory, double deltat)
{
    double frac;
    double ffrac;
    if(deltat < 0.0) return 0;
    double findex = (double) deltat / node_repetition_time;
    int trajectory_node = findex;
    if(trajectory_node >= (int) nodes_on_trajectory[trajectory]) return 0;

    if(trajectory_node+1 >= (int) nodes_on_trajectory[trajectory])
    {
        frac = 1.0 + findex - (double)trajectory_node;
        ffrac = (1.0 - frac);

        *x = ffrac*trajectory_nodes[trajectory][3*trajectory_node-3];
        *y = ffrac*trajectory_nodes[trajectory][3*trajectory_node-2];
        *z = ffrac*trajectory_nodes[trajectory][3*trajectory_node-1];

        *x += frac*trajectory_nodes[trajectory][3*trajectory_node];
        *y += frac*trajectory_nodes[trajectory][3*trajectory_node+1];
        *z += frac*trajectory_nodes[trajectory][3*trajectory_node+2];
    }
    else
    {
        frac = findex - (double)trajectory_node;
        ffrac = (1.0 - frac);

        *x = ffrac*trajectory_nodes[trajectory][3*trajectory_node];
        *y = ffrac*trajectory_nodes[trajectory][3*trajectory_node+1];
        *z = ffrac*trajectory_nodes[trajectory][3*trajectory_node+2];

        *x += frac*trajectory_nodes[trajectory][3*trajectory_node+3];
        *y += frac*trajectory_nodes[trajectory][3*trajectory_node+4];
        *z += frac*trajectory_nodes[trajectory][3*trajectory_node+5];
    }
    return 1;
}

// Gets coordinates of a (index) particle located on a given trajectory.
// Returns 1 on success, 0 in case particle outside the tube.
int czGetCoords(double current_time, unsigned int trajectory, unsigned int index, double* x, double* y, double* z)
{
    //Entry time of the index particle
    double deltat = ((int)index - last_in_tube_particle_index_at_start[trajectory])* particle_repetition_time[trajectory] - randomization_time[trajectory];
    //Time in tube after entry
    deltat = current_time - deltat;
    if(deltat < 0.0) return 0;
    if(InterpolateCoordinates(x, y, z, trajectory, deltat)) return 1;
    return 0;
}
