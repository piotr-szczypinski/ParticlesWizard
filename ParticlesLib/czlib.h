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

#ifndef CZLIB_H
#define CZLIB_H
// Loads file with tesselation data and trajectories. Returns true on success.
int czLoadTrajectories(const char* filename);

// Returns number of trajectories.
unsigned int czNumberOfTrajectories(void);

// Writes min and max of coordinates. Returns true on success.
void GetMinMaxXYZ(double* xb, double* xt, double* yb, double* yt, double* zb, double* zt);

// Writes min and max with minimum and maximum index of in-tube particle.
void czGetMinMaxInTube(double current_time, unsigned int trajectory, int* min, int* max);

// Gets volume of a particle on a given trajectory. Returns true on success.
int czGetVolume(unsigned int trajectory, double* v);

// Gets coordinates of a (index) particle located on a given trajectory.
// Returns true on success, false in case particle outside the tube.
int czGetCoords(double current_time, unsigned int trajectory, unsigned int index, double* x, double* y, double* z);
#endif
