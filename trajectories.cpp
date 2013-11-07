#include "trajectories.h"
#include <QtGlobal>
#include <stddef.h>
//#include <cmath>
#include <math.h>
#include <limits>

Trajectory::Trajectory()
{
    coordinates_nu = 0;
    coordinates = NULL;
    for(int i = 0; i < params_no; i++) params[i] = -1.0;
}
Trajectory::~Trajectory()
{
    clear();
}
void Trajectory::clear()
{
    if(coordinates != NULL) delete[] coordinates;
    coordinates = NULL;
    coordinates_nu = 0;
}
void Trajectory::setNumber(unsigned int i)
{
    coordinates_nu = i;
    if(coordinates != NULL) delete[] coordinates;
    coordinates = new Coordinates[i];
}
void Trajectory::resetNumber(unsigned int i)
{
    if(i < coordinates_nu) coordinates_nu = i;
}
unsigned int Trajectory::getNumber(void)
{
    return coordinates_nu;
}
bool Trajectory::setCoordinates(unsigned int i, double x[3])
{
    if(i < coordinates_nu && coordinates != NULL)
    {
        coordinates[i].x[0] = x[0];
        coordinates[i].x[1] = x[1];
        coordinates[i].x[2] = x[2];
        return true;
    }
    else return false;
}
bool Trajectory::getCoordinates(unsigned int i, double x[3])
{
    if(i < coordinates_nu && coordinates != NULL)
    {
        x[0] = coordinates[i].x[0];
        x[1] = coordinates[i].x[1];
        x[2] = coordinates[i].x[2];
        return true;
    }
    else return false;
}

Trajectories::Trajectories()
{
    trajectories_nu = 0;
    trajectories = NULL;
    stepTime = 0.0;
    ranges_ok = false;
    maxSteps = 0;
    for(int k = 0; k < 4; k++ ) inlet[k] = 0.0;

}
Trajectories::Trajectories(unsigned int k)
{
    trajectories_nu = k;
    if(k > 0) trajectories = new Trajectory[k];
    else trajectories = NULL;
    stepTime = 0.0;
    ranges_ok = false;

    projection[0][3] = 0.0;
    projection[0][4] = -1.0;
    projection[1][3] = 0.0;
    projection[1][4] = -1.0;
}
Trajectories::~Trajectories()
{
    if(trajectories != NULL) delete[] trajectories;
    trajectories_nu = 0;
}

void Trajectories::setNumber(unsigned int k)
{
    ranges_ok = false;
    trajectories_nu = k;
    if(trajectories != NULL) delete[] trajectories;
    if(k > 0) trajectories = new Trajectory[k];
    else trajectories = NULL;
}
unsigned int Trajectories::getNumber(void)
{
    return trajectories_nu;
}
bool Trajectories::setNumber(unsigned int k, unsigned int i)
{
    if(k < trajectories_nu)
    {
        ranges_ok = false;
        trajectories[k].setNumber(i);
        return true;
    }
    else return false;
}
bool Trajectories::resetNumber(unsigned int k, unsigned int i)
{
    if(k < trajectories_nu)
    {
        ranges_ok = false;
        trajectories[k].resetNumber(i);
        return true;
    }
    else return false;
}

unsigned int Trajectories::getNumber(unsigned int k)
{
    if(k < trajectories_nu)
    {
        return trajectories[k].getNumber();
    }
    else return 0;
}

unsigned int Trajectories::removeTrajectory(unsigned int k)
{
    if(k < trajectories_nu)
    {
        ranges_ok = false;
        trajectories_nu--;
        unsigned int i = k;

        trajectories[i].clear();
        Trajectory temp = trajectories[k];
        for(; i < trajectories_nu; i++)
        {
            trajectories[i] = trajectories[i+1];
        }
        trajectories[i] = temp;
        //trajectories[i].~Trajectory();
        return trajectories_nu;
    }
    else return 0;
}

bool Trajectories::setCoordinates(unsigned int k, unsigned int i, double x[3])
{
    ranges_ok = false;
    if(k < trajectories_nu) return trajectories[k].setCoordinates(i, x);
    else return false;
}
bool Trajectories::getCoordinates(unsigned int k, unsigned int i, double x[3])
{
    if(k < trajectories_nu) return trajectories[k].getCoordinates(i, x);
    else return false;

}

bool Trajectories::getTrajectoryParams(unsigned int k, double params[params_no])
{
    if(k < trajectories_nu)
    {
        for(int i = 0; i < params_no; i++) params[i] = trajectories[k].params[i];
        return true;
    }
    else return false;
}
bool Trajectories::setTrajectoryParams(unsigned int k, double params[params_no])
{
    if(k < trajectories_nu)
    {
        for(int i = 0; i < params_no; i++) trajectories[k].params[i] = params[i];
        return true;
    }
    else return false;
}


double Trajectories::getRange(unsigned int d, unsigned int l)
{
    double dd[3];
    if(d >= 3 || l >= 2) return 0.0;
    if(!ranges_ok)
    {
        for(unsigned int k = 0; k < trajectories_nu; k++)
        {
            unsigned int nodes = trajectories[k].getNumber();
            for(unsigned int i = 0; i < nodes; i++)
            {
                trajectories[k].getCoordinates(i, dd);
                if(i == 0 && k == 0)
                {
                    for(unsigned int ii = 0; ii < 3; ii++)
                    {
                        ranges[ii][0] = dd[ii];
                        ranges[ii][1] = dd[ii];
                    }
                }
                else
                {
                    for(unsigned int ii = 0; ii < 3; ii++)
                    {
                        if(ranges[ii][0] > dd[ii]) ranges[ii][0] = dd[ii];
                        if(ranges[ii][1] < dd[ii]) ranges[ii][1] = dd[ii];
                    }
                }
            }
        }
        ranges_ok = true;
    }

    return ranges[d][l];
}

bool Trajectories::setInlet(double c[3], double d)
{
    bool ret = false;

    inlet[0] = c[0];
    inlet[1] = c[1];
    inlet[2] = c[2];
    inlet[3] = d;

    double nanss[3];
    nanss[0] = nanss[1] = nanss[2] = std::numeric_limits<double>::quiet_NaN();//nan("");//qQNaN();
    //nan();

    double n = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
    n = sqrt(n);
    if(trajectories_nu > 0)
    {
        inletpoints.setNumber(trajectories_nu);

        for(unsigned int k = 0; k < trajectories_nu; k++)
        {
            trajectories[k].params[params_distance] = -1;

            unsigned int i = 0;
            unsigned int nodes = trajectories[k].getNumber();
            for(i = 0; i < nodes; i++)
            {
                double dd[3];
                trajectories[k].getCoordinates(i, dd);
                double f = inlet[0]*dd[0] + inlet[1]*dd[1] + inlet[2]*dd[2] - inlet[3];
                if(f > 0)
                {
                    if(i == 0)
                    {
                        inletpoints.setCoordinates(k, nanss);
                    }
                    else
                    {
                        double ddd[3];
                        trajectories[k].getCoordinates(i-1, ddd);

                        double distt;
                        double dist[3];
                        dist[0] = ddd[0] - dd[0];
                        dist[1] = ddd[1] - dd[1];
                        dist[2] = ddd[2] - dd[2];
                        distt = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];
                        distt = sqrt(distt);
                        double div = distt * n;

                        if(div > 0)
                        {
                            double tilt = (c[0] * dist[0] + c[1] * dist[1] + c[2] * dist[2]) / div;
                            trajectories[k].params[params_tilt] = -tilt;
                        }
                        else
                        {
                            trajectories[k].params[params_tilt] = 1.0;
                        }
                        trajectories[k].params[params_distance] = distt;

                        double ff = inlet[0]*ddd[0] + inlet[1]*ddd[1] + inlet[2]*ddd[2] - inlet[3];
                        div = (f - ff);
                        if(div > 0.0000001)
                        {
                            dd[0] *= (-ff); dd[1] *= (-ff); dd[2] *= (-ff);
                            ddd[0] *= (f); ddd[1] *= (f); ddd[2] *= (f);
                            dd[0] += ddd[0]; dd[1] += ddd[1]; dd[2] += ddd[2];
                            dd[0] /= div; dd[1] /= div; dd[2] /= div;
                            ret = true;
                        }
                        inletpoints.setCoordinates(k, dd);
                    }
                    break;
                }
            }
            if(i >= nodes)
            {
                inletpoints.setCoordinates(k, nanss);
            }
        }
    }
    return ret;
}

void Trajectories::getInlet(double c[3], double* d)
{
    c[0] = inlet[0];
    c[1] = inlet[1];
    c[2] = inlet[2];
    *d = inlet[3];
}


void Trajectories::getInletPoint(unsigned int i, double c[3])
{
     if(inletpoints.getNumber() > i)
     {
        inletpoints.getCoordinates(i, c);
     }
}

bool Trajectories::getProjectedInletPoint(unsigned int i, double c[2])
{
    double cc[3];
    if(inletpoints.getNumber() <= i) return false;
    if(!inletpoints.getCoordinates(i, cc)) return false;

    c[0] = projection[0][0] * cc[0] + projection[0][1] * cc[1] + projection[0][2] * cc[2];
    c[1] = projection[1][0] * cc[0] + projection[1][1] * cc[1] + projection[1][2] * cc[2];
    return true;
}

void Trajectories::get2DRanges(double x[2], double y[2])
{
    x[0] = projection[0][3];
    x[1] = projection[0][4];
    y[0] = projection[1][3];
    y[1] = projection[1][4];
}

void Trajectories::computeProjection(void)
{
    double n = sqrt(inlet[0]*inlet[0] + inlet[2]*inlet[2]);
    if(n < 0.0000001)
    {
        n = sqrt(inlet[0]*inlet[0] + inlet[1]*inlet[1]);
        if(n < 0.0000001)
        {
            projection[0][3] = 0.0;
            projection[0][4] = -1.0;
            projection[1][3] = 0.0;
            projection[1][4] = -1.0;
            return;
        }
        projection[0][0] = -inlet[1]/n;
        projection[0][1] = 0.0;
        projection[0][2] = -inlet[0]/n;
    }
    else
    {
        projection[0][0] = inlet[2]/n;
        projection[0][1] = 0.0;
        projection[0][2] = -inlet[0]/n;
    }
    projection[1][0] = projection[0][1]*inlet[2] - projection[0][2]*inlet[1];
    projection[1][1] = projection[0][2]*inlet[0] - projection[0][0]*inlet[2];
    projection[1][2] = projection[0][0]*inlet[1] - projection[0][1]*inlet[0];

    n = sqrt(projection[1][0] * projection[1][0] + projection[1][1] * projection[1][1] + projection[1][2] * projection[1][2]);

    projection[1][0] /= n;
    projection[1][1] /= n;
    projection[1][2] /= n;

    bool first = true;
    int maxi = inletpoints.getNumber();
    double xmin, xmax, ymin, ymax;
    for(unsigned int i = 0; i < maxi; i++)
    {
        double c[2];
        if(getProjectedInletPoint(i, c))
        {
            if(first)
            {
                xmin = xmax = c[0];
                ymin = ymax = c[1];
                first = false;
            }
            else
            {
                if(xmax < c[0]) xmax = c[0];
                if(ymax < c[1]) ymax = c[1];
                if(xmin > c[0]) xmin = c[0];
                if(ymin > c[1]) ymin = c[1];
            }
        }
    }
    if(first)
    {
        projection[0][3] = 0.0;
        projection[0][4] = -1.0;
        projection[1][3] = 0.0;
        projection[1][4] = -1.0;
    }
    else
    {
        projection[0][3] = xmin;
        projection[0][4] = xmax;
        projection[1][3] = ymin;
        projection[1][4] = ymax;
    }
}



