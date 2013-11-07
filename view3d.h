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

#ifndef VIEW3D_H
#define VIEW3D_H

#include <QWidget>
#include "trajectories.h"

class View3D : public QWidget
{
    Q_OBJECT
public:
    explicit View3D(QWidget *parent = 0);
    void setTrajectories(Trajectories *tr);
    void setCurrentStep(unsigned int step);
    void setShowTrajectories(bool b);
    virtual void Render(void);
    double shift[2];
    int mousex;
    int mousey;
    double getZoom(void);
signals:
    
public slots:

private:
    int zoom;
    double W[3][3];
    double WI[3][3];

protected:
    Trajectories* trajectories;
    double resizer;
    unsigned int current_step;
    bool show_trajectories;
    double shift_before[3];
    double shift_after[2];

    void ToScreen(double out[2], double in[2]);
    void FromScreen(double out[2], double in[2]);
    void Project(double *X, double *Y, double x[]);
    void Rotate(int rot, int l1, int l2);

    void wheelEvent(QWheelEvent *e);
    void paintEvent(QPaintEvent *e);
    virtual void mouseMoveEvent(QMouseEvent *e);
    virtual void mousePressEvent(QMouseEvent *e);
};

class View3D_Cross : public View3D
{
public:
    explicit View3D_Cross(QWidget *parent = 0);
    void Render(void);
};

class View3D_Verify : public View3D
{
public:
    explicit View3D_Verify(QWidget *parent = 0);
    void Render(void);
    void setCurrentTime(double t);
private:
    double current_time;
    bool GetCoords(int last_in_tube_particle_index_at_start, int nodes_on_trajectory, unsigned int trajectory, unsigned int index, double params[params_no], double x[3]);
};

extern "C"
{
#include "vvoronoi.h"
}

class View3D_Tessellate : public View3D
{
    Q_OBJECT
public:
    explicit View3D_Tessellate(QWidget *parent = 0);
    void Render(void);
    void setVoronoi(VoronoiData *vd);
signals:
    void radiusChanged(double r);
    void centerXChanged(double x);
    void centerYChanged(double y);

public slots:
    void radiusChangedIn(double r);
    void centerXChangedIn(double x);
    void centerYChangedIn(double y);

private:
    double radius;
    double center[2];
    void mouseMoveEvent(QMouseEvent *e);
    void mousePressEvent(QMouseEvent *e);
    VoronoiData *voronoidata;
};

#endif // VIEW3D_H
