/*
 * ParticlesWizard
 *
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

#ifndef TRAJECTORIES_H
#define TRAJECTORIES_H

typedef struct Coordinates { double x[3]; } Coordinates;
#define params_no 6
#define params_area 0
#define params_volume 1
#define params_repetition 2
#define params_random 3
#define params_distance 4
#define params_tilt 5


class Trajectory
{
public:
    Trajectory();
    ~Trajectory();
    void clear();
    void setNumber(unsigned int i);
    void resetNumber(unsigned int i);
    unsigned int getNumber(void);
    bool setCoordinates(unsigned int i, double x[3]);
    bool getCoordinates(unsigned int i, double x[3]);

    double params[params_no];//area, volume; repetition; timeshift, distance;
private:
    unsigned int coordinates_nu;
    Coordinates *coordinates;
    double volume;
    double repetition;
};



class Trajectories
{
public:
    Trajectories();
    Trajectories(unsigned int k);
    ~Trajectories();
    bool setInlet(double c[3], double d);
    void getInlet(double c[3], double* d);
    void getInletPoint(unsigned int i, double c[3]);

    bool getTrajectoryParams(unsigned int k, double params[4]);
    bool setTrajectoryParams(unsigned int k, double params[4]);
    double stepTime;
    unsigned int maxSteps;
    void setNumber(unsigned int k);
    bool setNumber(unsigned int k, unsigned int i);
    bool resetNumber(unsigned int k, unsigned int i);
    unsigned int getNumber(void);
    unsigned int getNumber(unsigned int k);
    unsigned int removeTrajectory(unsigned int k);
    bool setCoordinates(unsigned int k, unsigned int i, double x[3]);
    bool getCoordinates(unsigned int k, unsigned int i, double x[3]);
    double getRange(unsigned int d, unsigned int l);

    bool getProjectedInletPoint(unsigned int i, double c[2]);
    void get2DRanges(double x[2], double y[2]);
    void computeProjection(void);
private:
    unsigned int trajectories_nu;
    Trajectory* trajectories;
    Trajectory inletpoints;
    double ranges[3][2];
    bool ranges_ok;
    double inlet[4];

    double projection[2][5];
};

#endif // TRAJECTORIES_H
