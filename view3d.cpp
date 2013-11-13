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


#include "view3d.h"

#include <math.h>

#include <QRect>
#include <QPainter>
#include <QMouseEvent>
#include <QWheelEvent>


#define ALMOST_ZERO 0.0000000001
#define MAX_ROIS 16

const unsigned int RoiColors[MAX_ROIS]=
{
    0x000000FF, 0x0000FF00, 0x00FF0000,
    0x00FFFF00, 0x00FF00FF, 0x0000FFFF,
    0x000080FF, 0x008000FF, 0x0000FF80,
    0x0080FF00, 0x00FF0080, 0x00FF8000,
    0x0000C4FF, 0x0000FFC4, 0x00FF00C4,
    0x00C4FF00
};
/*------------------------------------------------------------------------*/
double podwyznacznik(double M[2][2], int w, int k)
{
    int w0;
    int k0;

    w0=(w+1)%2;
    k0=(k+1)%2;
    if((w+k)%2) return -M[w0][k0];
    return M[w0][k0];
}

/*------------------------------------------------------------------------*/
double podwyznacznik(double M[3][3], int w, int k)
{
    int w0, w1;
    int k0, k1;

    w0=0; if(w0==w) w0++;
    w1=w0+1; if(w1==w) w1++;
    k0=0; if(k0==k) k0++;
    k1=k0+1; if(k1==k) k1++;

    double r=
            M[w0][k0]*M[w1][k1]-M[w0][k1]*M[w1][k0];

    if((w+k)%2) return -r;
    return r;
}

/*------------------------------------------------------------------------*/
double podwyznacznik(double M[4][4], int w, int k)
{
    int w0, w1, w2;
    int k0, k1, k2;

    w0=0; if(w0==w) w0++;
    w1=w0+1; if(w1==w) w1++;
    w2=w1+1; if(w2==w) w2++;
    k0=0; if(k0==k) k0++;
    k1=k0+1; if(k1==k) k1++;
    k2=k1+1; if(k2==k) k2++;

    double r=
            M[w0][k0]*(M[w1][k1]*M[w2][k2]-M[w1][k2]*M[w2][k1])
            -M[w0][k1]*(M[w1][k0]*M[w2][k2]-M[w1][k2]*M[w2][k0])
            +M[w0][k2]*(M[w1][k0]*M[w2][k1]-M[w1][k1]*M[w2][k0]);

    if((w+k)%2) return -r;
    return r;
}

/*------------------------------------------------------------------------*/
void transpozycja(double MM[4][4], double M[4][4])
{
    for(int k=1; k<4; k++)
        for(int w=0; w<k; w++)
        {
            MM[w][k]=M[k][w];
            MM[k][w]=M[w][k];
        }
    for(int w=0; w<4; w++)
        MM[w][w]=M[w][w];
}

/*------------------------------------------------------------------------*/
void transpozycja(double MM[4][4])
{
    for(int k=1; k<4; k++)
        for(int w=0; w<k; w++)
        {
            double d =MM[w][k];
            MM[w][k]=MM[k][w];
            MM[k][w]=d;
        }
}
/*------------------------------------------------------------------------*/
void transpozycja(double MM[3][3], double M[3][3])
{
    for(int k=1; k<3; k++)
        for(int w=0; w<k; w++)
        {
            MM[w][k]=M[k][w];
            MM[k][w]=M[w][k];
        }
    for(int w=0; w<3; w++)
        MM[w][w]=M[w][w];
}

/*------------------------------------------------------------------------*/
void transpozycja(double MM[3][3])
{
    for(int k=1; k<3; k++)
        for(int w=0; w<k; w++)
        {
            double d =MM[w][k];
            MM[w][k]=MM[k][w];
            MM[k][w]=d;
        }
}

/*------------------------------------------------------------------------*/
void transpozycja(double MM[2][2], double M[2][2])
{
    MM[0][1]=M[1][0];
    MM[1][0]=M[0][1];
    MM[0][0]=M[0][0];
    MM[1][1]=M[1][1];
}

/*------------------------------------------------------------------------*/
void transpozycja(double MM[2][2])
{
    double d =MM[0][1];
    MM[0][1]=MM[1][0];
    MM[1][0]=d;
}


/*------------------------------------------------------------------------*/
double wyznacznik(double M[4][4])
{
    return
            M[0][0]*
            (
            M[1][1]*(M[2][2]*M[3][3]-M[2][3]*M[3][2])
            -M[1][2]*(M[2][1]*M[3][3]-M[2][3]*M[3][1])
            +M[1][3]*(M[2][1]*M[3][2]-M[2][2]*M[3][1])
            )
            -M[0][1]*
            (
            M[1][0]*(M[2][2]*M[3][3]-M[2][3]*M[3][2])
            -M[1][2]*(M[2][0]*M[3][3]-M[2][3]*M[3][0])
            +M[1][3]*(M[2][0]*M[3][2]-M[2][2]*M[3][0])
            )
            +M[0][2]*
            (
            M[1][0]*(M[2][1]*M[3][3]-M[2][3]*M[3][1])
            -M[1][1]*(M[2][0]*M[3][3]-M[2][3]*M[3][0])
            +M[1][3]*(M[2][0]*M[3][1]-M[2][1]*M[3][0])
            )
            -M[0][3]*
            (
            M[1][0]*(M[2][1]*M[3][2]-M[2][2]*M[3][1])
            -M[1][1]*(M[2][0]*M[3][2]-M[2][2]*M[3][0])
            +M[1][2]*(M[2][0]*M[3][1]-M[2][1]*M[3][0])
            );
}

/*------------------------------------------------------------------------*/
double wyznacznik(double M[3][3])
{
    return
            M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
            -M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
            +M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
}

/*------------------------------------------------------------------------*/
double wyznacznik(double M[2][2])
{
    return
            M[0][0]*M[1][1]-M[0][1]*M[1][0];
}

/*------------------------------------------------------------------------*/
double inwersja(double MM[4][4], double M[4][4])
{
    double det = wyznacznik(M);
    if(det==0.0) return det;

    for(int k=0; k<4; k++)
        for(int w=0; w<4; w++)
        {
            MM[w][k] = podwyznacznik(M, w, k)/det;
        }
    transpozycja(MM);
    return det;
}

/*------------------------------------------------------------------------*/
double inwersja(double MM[3][3], double M[3][3])
{
    double det = wyznacznik(M);
    if(det==0.0) return det;

    for(int k=0; k<3; k++)
        for(int w=0; w<3; w++)
        {
            MM[w][k] = podwyznacznik(M, w, k)/det;
        }
    transpozycja(MM);
    return det;
}

/*------------------------------------------------------------------------*/
bool inwersja(double MM[2][2], double M[2][2])
{
    double det = wyznacznik(M);
    if(fabs(det)<ALMOST_ZERO) return false;
    for(int k=0; k<2; k++)
        for(int w=0; w<2; w++)
        {
            MM[w][k] = podwyznacznik(M, w, k)/det;
        }
    transpozycja(MM);
    return true;
}






//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
View3D::View3D(QWidget *parent) :
    QWidget(parent)
{
    zoom = 0;
    shift[0] = 0.0;
    shift[1] = 0.0;

    show_trajectories = true;
    trajectories = NULL;
    current_step = 0;

    W[0][0] = 1; W[0][1] = 0; W[0][2] = 0;
    W[1][0] = 0; W[1][1] = 1; W[1][2] = 0;
    W[2][0] = 0; W[2][1] = 0; W[2][2] = 1;

    WI[0][0] = 1; WI[0][1] = 0; WI[0][2] = 0;
    WI[1][0] = 0; WI[1][1] = 1; WI[1][2] = 0;
    WI[2][0] = 0; WI[2][1] = 0; WI[2][2] = 1;
}

//---------------------------------------------------------------------------
void View3D::setTrajectories(Trajectories* tr)
{
    trajectories = tr;
    update();
}

//---------------------------------------------------------------------------
void View3D::setCurrentStep(unsigned int step)
{
    current_step = step;
    update();
}

//---------------------------------------------------------------------------
void View3D::setShowTrajectories(bool b)
{
    show_trajectories = b;
    update();
}

//---------------------------------------------------------------------------
void View3D::paintEvent(QPaintEvent *e)
{
    Render();
    return;

}

//---------------------------------------------------------------------------
void View3D::ToScreen(double out[2], double in[2])
{
    out[0] = in[0] + shift_before[0];
    out[1] = in[1] + shift_before[1];
    out[0] *= resizer;
    out[1] *= resizer;
    out[0] *= getZoom();
    out[1] *= getZoom();
    out[0] += shift_after[0];
    out[1] += shift_after[1];
    out[0] += shift[0];
    out[1] += shift[1];
}

//---------------------------------------------------------------------------
void View3D::FromScreen(double out[2], double in[2])
{
    out[0] = in[0] - shift[0];
    out[1] = in[1] - shift[1];
    out[0] -= shift_after[0];
    out[1] -= shift_after[1];
    out[0] /= getZoom();
    out[1] /= getZoom();
    out[0] /= resizer;
    out[1] /= resizer;
    out[0] -= shift_before[0];
    out[1] -= shift_before[1];
}

//---------------------------------------------------------------------------
void View3D::Project(double *X, double *Y, double x[3])
{
    x[0] += shift_before[0];
    x[1] += shift_before[1];
    x[2] += shift_before[2];
    *X = WI[0][0]*x[0] + WI[0][1]*x[1] + WI[0][2]*x[2];
    *Y = WI[1][0]*x[0] + WI[1][1]*x[1] + WI[1][2]*x[2];
    *X *= resizer;
    *Y *= resizer;
    *X *= getZoom();
    *Y *= getZoom();
    *X += shift_after[0];
    *Y += shift_after[1];
    *X += shift[0];
    *Y += shift[1];
}

//---------------------------------------------------------------------------
void View3D::Rotate(int rot, int l1, int l2)
{
    double a, b;
    double sina = sin((double)rot/100);
    double cosa = cos((double)rot/100);

    a = cosa * W[0][l1] + sina * W[0][l2];
    b = -sina * W[0][l1] + cosa * W[0][l2];
    W[0][l1] = a;
    W[0][l2] = b;
    a = cosa * W[1][l1] + sina * W[1][l2];
    b = -sina * W[1][l1] + cosa * W[1][l2];
    W[1][l1] = a;
    W[1][l2] = b;
    a = cosa * W[2][l1] + sina * W[2][l2];
    b = -sina * W[2][l1] + cosa * W[2][l2];
    W[2][l1] = a;
    W[2][l2] = b;
    inwersja(WI, W);
}

//---------------------------------------------------------------------------
double View3D::getZoom(void)
{
    return pow(12.0, (double)zoom/2400);
}

//---------------------------------------------------------------------------
void View3D::wheelEvent(QWheelEvent *e)
{
    QRect viewport = geometry();

    double rx = abs((viewport.right()-viewport.left()) / 2);
    double ry = abs((viewport.top()-viewport.bottom()) / 2);
    int cx = e->pos().x() - rx;
    int cy = e->pos().y() - ry;

    double zoomb = getZoom();
    e->accept();
    zoom += e->delta();
    if(zoom > 2400) zoom = 2400;
    if(zoom < -2400) zoom = -2400;
    double zooma = getZoom();
    rx *= zooma;
    ry *= zooma;
    zooma /= zoomb;
    shift[0] = (shift[0] - cx) * zooma + cx;
    shift[1] = (shift[1] - cy) * zooma + cy;
    if(shift[0] > rx) shift[0] = rx;
    if(shift[1] > ry) shift[1] = ry;
    if(shift[0] < -rx) shift[0] = -rx;
    if(shift[1] < -ry) shift[1] = -ry;

    update();
}

//---------------------------------------------------------------------------
void View3D::mouseMoveEvent(QMouseEvent *e)
{

    if(e->buttons() & Qt::LeftButton)
    {
        e->accept();
        Rotate(e->x()-mousex, 0, 2);
        Rotate(e->y()-mousey, 1, 2);
        mousex = e->x();
        mousey = e->y();
        update();
    }
    else if(e->buttons() & Qt::RightButton)
    {
        e->accept();
        QRect viewport = geometry();

        int cx = e->pos().x() - mousex;
        int cy = e->pos().y() - mousey;

        shift[0] = shift[0] + cx;
        shift[1] = shift[1] + cy;

        double zooma = getZoom();
        double rx = abs((viewport.right()-viewport.left()) / 2);
        double ry = abs((viewport.top()-viewport.bottom()) / 2);
        rx *= zooma;
        ry *= zooma;
        if(shift[0] > rx) shift[0] = rx;
        if(shift[1] > ry) shift[1] = ry;
        if(shift[0] < -rx) shift[0] = -rx;
        if(shift[1] < -ry) shift[1] = -ry;

        mousex = e->x();
        mousey = e->y();
        update();
    }
}

//---------------------------------------------------------------------------
void View3D::mousePressEvent(QMouseEvent *e)
{
    e->accept();
    mousex = e->x();
    mousey = e->y();
}

//---------------------------------------------------------------------------
void View3D::Render(void)
{
}















//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

View3D_Cross::View3D_Cross(QWidget *parent) :
    View3D(parent)
{
}

//---------------------------------------------------------------------------
void View3D_Cross::Render(void)
{
    double x[3];
    double X, Y, Xo, Yo;
    int step;

    double mn;
    QPainter painter(this);
    if(trajectories == NULL) return;
    QRect viewport = geometry();

    shift_after[0] = ((viewport.right()+viewport.left())/2);
    shift_after[1] = ((viewport.top()+viewport.bottom())/2);

    shift_before[0] = -(trajectories->getRange(0,1) + trajectories->getRange(0,0))/2.0;
    shift_before[1] = -(trajectories->getRange(1,1) + trajectories->getRange(1,0))/2.0;
    shift_before[2] = -(trajectories->getRange(2,1) + trajectories->getRange(2,0))/2.0;

    mn = trajectories->getRange(0,1) - trajectories->getRange(0, 0); resizer = mn;
    mn = trajectories->getRange(1,1) - trajectories->getRange(1, 0); if(resizer<mn) resizer = mn;
    mn = trajectories->getRange(2,1) - trajectories->getRange(2, 0); if(resizer<mn) resizer = mn;
    if(resizer<ALMOST_ZERO) resizer = ALMOST_ZERO;
    resizer = (double)(viewport.width()) / resizer;

    int liczba_czastek = trajectories->getNumber();
    int liczba_krokow = trajectories->maxSteps;

    if(show_trajectories)
    {
        QPen linepen(Qt::green);
        linepen.setWidth(1);
        painter.setPen(linepen);


        for(int part = 0; part < liczba_czastek; part++)
        {
            if(trajectories->getCoordinates(part, 0, x))
            {
                Project(&Xo, &Yo, x);
                for(step = 1; step < liczba_krokow; step++)
                {
                    if(!trajectories->getCoordinates(part, step, x)) break;
                    Project(&X, &Y, x);
                    painter.drawLine(Xo, Yo, X, Y);
                    Xo = X; Yo = Y;
                }
            }
        }
    }

    QPen pointpen(Qt::black);
    pointpen.setWidth(2);
    painter.setPen(pointpen);
    //liczba_krokow = ui->horizontalSlider->maximum();
    step = current_step;
    for(int part = 0; part < liczba_czastek; part++)
    {
        if(trajectories->getCoordinates(part, step, x))
        {
            Project(&X, &Y, x);
            painter.drawPoint(X, Y);
        }
    }

    pointpen.setColor(Qt::red);
    painter.setPen(pointpen);
    for(int part = 0; part < liczba_czastek; part++)
    {
        trajectories->getInletPoint(part, x);
        if(x[0] == x[0])
        {
            Project(&X, &Y, x);
            painter.drawPoint(X, Y);
        }
    }
}



//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

View3D_Tessellate::View3D_Tessellate(QWidget *parent) :
    View3D(parent)
{
    radius = -1;
    voronoidata = NULL;
}

void View3D_Tessellate::setVoronoi(VoronoiData *vd)
{
    voronoidata = vd;

}

//---------------------------------------------------------------------------
void View3D_Tessellate::mousePressEvent(QMouseEvent *e)
{
    double in[2];
    double out[2];
    e->accept();
    in[0] = mousex = e->x();
    in[1] = mousey = e->y();

    if(e->buttons() & Qt::LeftButton)
    {
        FromScreen(out, in);
        emit centerXChanged(out[0]);
        emit centerYChanged(out[1]);
        emit radiusChanged(0.0);
    }
}

//---------------------------------------------------------------------------
void View3D_Tessellate::mouseMoveEvent(QMouseEvent *e)
{
    if(e->buttons() & Qt::LeftButton)
    {
        double x, y;
        e->accept();
        x = e->x() - mousex;
        y = e->y() - mousey;
        double radius = sqrt(x*x + y*y);
        radius /= getZoom();
        radius /= resizer;
        emit radiusChanged(radius);
        update();
    }
    else if(e->buttons() & Qt::RightButton)
    {
        e->accept();
        QRect viewport = geometry();

        int cx = e->pos().x() - mousex;
        int cy = e->pos().y() - mousey;

        shift[0] = shift[0] + cx;
        shift[1] = shift[1] + cy;

        double zooma = getZoom();
        double rx = abs((viewport.right()-viewport.left()) / 2);
        double ry = abs((viewport.top()-viewport.bottom()) / 2);
        rx *= zooma;
        ry *= zooma;
        if(shift[0] > rx) shift[0] = rx;
        if(shift[1] > ry) shift[1] = ry;
        if(shift[0] < -rx) shift[0] = -rx;
        if(shift[1] < -ry) shift[1] = -ry;

        mousex = e->x();
        mousey = e->y();
        update();
    }
}
//---------------------------------------------------------------------------
void View3D_Tessellate::Render(void)
{
    double x[2];
    double out[2];
    double out2[2];
    double mn;

    //bool voronoi_valid = false;
    QPainter painter(this);
    if(trajectories == NULL) return;
    QRect viewport = geometry();

    shift_after[0] = (viewport.right()+viewport.left())/2;
    shift_after[1] = (viewport.top()+viewport.bottom())/2;

    double rx[2];
    double ry[2];

    trajectories->get2DRanges(rx, ry);
    if(rx[1] <= rx[0] || ry[1] <= ry[0]) return;


    shift_before[0] = -(rx[1] + rx[0])/2.0;
    shift_before[1] = -(ry[1] + ry[0])/2.0;

    mn = rx[1] - rx[0]; resizer = mn;
    mn = ry[1] - ry[0]; if(resizer<mn) resizer = mn;

    if(resizer<ALMOST_ZERO) resizer = ALMOST_ZERO;
    resizer = (double)(viewport.width()) / resizer;


    if(voronoidata != NULL)
    {
        if(voronoidata->validcellsno > 0 && voronoidata->validcellsno <= voronoidata->cellsno)
        {
            QPen pointpen(Qt::cyan);
            pointpen.setWidth(1);
            painter.setPen(pointpen);

            for(int i = 0; i < voronoidata->validcellsno; i++)
            {
                int maxk = voronoidata->cells[i][0];
                int il, ir;
                il = voronoidata->cells[i][maxk];
                for(int k = 1; k <= maxk; k++)
                {
                    ir = voronoidata->cells[i][k];
                    if(ir > 0 && il > 0)
                    {
                        x[0] = voronoidata->vertices[ir*2];
                        x[1] = voronoidata->vertices[ir*2+1];
                        ToScreen(out, x);
                        x[0] = voronoidata->vertices[il*2];
                        x[1] = voronoidata->vertices[il*2+1];
                        ToScreen(out2, x);
                        painter.drawLine(QPoint(out[0], out[1]), QPoint(out2[0], out2[1]));
                    }
                    il = ir;
                }
            }
        }
    }


    int liczba_czastek = trajectories->getNumber();


    double params[params_no];


    //QPen pointpengreen(Qt::darkGreen);
    QPen pointpengreen(Qt::green);
    pointpengreen.setWidth(2);
    QPen pointpen(Qt::red);
    pointpen.setWidth(2);
    for(int part = 0; part < liczba_czastek; part++)
    {
        if(trajectories->getProjectedInletPoint(part, x))
        {
            trajectories->getTrajectoryParams(part, params);
            if(params[params_area] >= 0) painter.setPen(pointpengreen);
            else  painter.setPen(pointpen);

            ToScreen(out, x);
            painter.drawPoint(out[0], out[1]);
        }
    }

    if(radius > 0)
    {
        pointpen.setColor(Qt::blue);
        pointpen.setWidth(1);
        painter.setPen(pointpen);
        ToScreen(out, center);

        double r = radius * getZoom();
        r *= resizer;
        painter.drawEllipse(QPoint(out[0], out[1]), (int)r, (int)r);
    }







    /*
    pointpen.setColor(Qt::red);
    painter.setPen(pointpen);
    for(int part = 0; part < liczba_czastek; part++)
    {
        trajectories->getInletPoint(part, x);
        if(x[0] == x[0])
        {
            Project(&X, &Y, x);
            painter.drawPoint(X, Y);
        }
    }
*/
}


void View3D_Tessellate::radiusChangedIn(double r)
{
    radius = r;
    update();
}
void View3D_Tessellate::centerXChangedIn(double x)
{
    center[0] = x;
    update();
}
void View3D_Tessellate::centerYChangedIn(double y)
{
    center[1] = y;
    update();
}




//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

View3D_Verify::View3D_Verify(QWidget *parent) :
    View3D(parent)
{
}

void View3D_Verify::setCurrentTime(double t)
{
    current_time = t;
}


// Gets coordinates of a (index) particle located on a given trajectory.
// Returns 1 on success, 0 in case particle outside the tube.
bool View3D_Verify::GetCoords(int last_in_tube_particle_index_at_start, int nodes_on_trajectory, unsigned int trajectory, unsigned int index, double params[params_no], double x[3])
{
    //Entry time of the index particle
    double deltat = ((int)index - last_in_tube_particle_index_at_start)* params[params_repetition] - params[params_random];
    //Time in tube after entry
    deltat = current_time - deltat;
    if(deltat < 0.0) return false;

    double frac;
    double ffrac;
    //if(deltat < 0.0) return false;
    double findex = (double) deltat / trajectories->stepTime;
    int trajectory_node = findex;

    if(trajectory_node >= (int) nodes_on_trajectory) return false;

    if(trajectory_node+1 >= (int) nodes_on_trajectory)
    {
        frac = 1.0 + findex - (double)trajectory_node;
        ffrac = (1.0 - frac);

        double xx[3];
        trajectories->getCoordinates(trajectory, trajectory_node-1, xx);
        x[0] = ffrac*xx[0]; x[1] = ffrac*xx[1]; x[2] = ffrac*xx[2];
        trajectories->getCoordinates(trajectory, trajectory_node, xx);
        x[0] += frac*xx[0]; x[1] += frac*xx[1]; x[2] += frac*xx[2];
    }
    else
    {
        frac = findex - (double)trajectory_node;
        ffrac = (1.0 - frac);

        double xx[3];
        trajectories->getCoordinates(trajectory, trajectory_node, xx);
        x[0] = ffrac*xx[0]; x[1] = ffrac*xx[1]; x[2] = ffrac*xx[2];
        trajectories->getCoordinates(trajectory, trajectory_node+1, xx);
        x[0] += frac*xx[0]; x[1] += frac*xx[1]; x[2] += frac*xx[2];
    }
    return true;
}

//---------------------------------------------------------------------------
void View3D_Verify::Render(void)
{
    double X, Y;
    double mn;

    QColor color;
    QPen pointpen;
    pointpen.setWidth(3);
    QPainter painter(this);
    if(trajectories == NULL) return;
    QRect viewport = geometry();

    shift_after[0] = ((viewport.right()+viewport.left())/2);
    shift_after[1] = ((viewport.top()+viewport.bottom())/2);

    shift_before[0] = -(trajectories->getRange(0,1) + trajectories->getRange(0,0))/2.0;
    shift_before[1] = -(trajectories->getRange(1,1) + trajectories->getRange(1,0))/2.0;
    shift_before[2] = -(trajectories->getRange(2,1) + trajectories->getRange(2,0))/2.0;

    mn = trajectories->getRange(0,1) - trajectories->getRange(0, 0); resizer = mn;
    mn = trajectories->getRange(1,1) - trajectories->getRange(1, 0); if(resizer<mn) resizer = mn;
    mn = trajectories->getRange(2,1) - trajectories->getRange(2, 0); if(resizer<mn) resizer = mn;
    if(resizer<ALMOST_ZERO) resizer = ALMOST_ZERO;
    resizer = (double)(viewport.width()) / resizer;

    int max_traj = trajectories->getNumber();
    int liczba_krokow = trajectories->maxSteps;


    if(true)
    {
        double x[3];
        double Xo, Yo;
        QPen linepen(Qt::gray);
        linepen.setWidth(1);
        painter.setPen(linepen);


        for(int part = 0; part < max_traj; part++)
        {
            if(trajectories->getCoordinates(part, 0, x))
            {
                Project(&Xo, &Yo, x);
                for(int step = 1; step < liczba_krokow; step++)
                {
                    if(!trajectories->getCoordinates(part, step, x)) break;
                    Project(&X, &Y, x);
                    painter.drawLine(Xo, Yo, X, Y);
                    Xo = X; Yo = Y;
                }
            }
        }
    }

    for(int traj = 0; traj < max_traj; traj++)
    {
        double params[params_no];
        trajectories->getTrajectoryParams(traj, params);
        if(params[params_area] <= 0) continue;

        int nodes_on_trajectory = trajectories->getNumber(traj);
        int last_in_tube_particle_index_at_start = (nodes_on_trajectory * trajectories->stepTime - params[params_random]) / params[params_repetition];

        int index_range[2];
        index_range[1] = last_in_tube_particle_index_at_start + (current_time + params[params_random]) / params[params_repetition];
        index_range[0] = last_in_tube_particle_index_at_start +
                (int)(
                    (
                        current_time + params[params_random] - (trajectories->getNumber(traj))*trajectories->stepTime
                        ) / params[params_repetition]
                    );

        if(index_range[0] < 0) index_range[0] = 0;
        if(index_range[1] < 0) index_range[1] = 0;

        for(int part = index_range[0]-1; part <= index_range[1]+1; part++)
        {
            double x[3];

            if(GetCoords(last_in_tube_particle_index_at_start, nodes_on_trajectory, traj, part, params, x))
            {
                Project(&X, &Y, x);
                color.setRgb((QRgb) RoiColors[part%MAX_ROIS]);
                pointpen.setColor(color);
                painter.setPen(pointpen);
                painter.drawPoint(X, Y);
            }
        }
    }
}
