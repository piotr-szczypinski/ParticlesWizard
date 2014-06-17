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

#include <QString>
#include <QFileDialog>

#include <math.h>

#include "wizard.h"
#include "ui_wizard.h"
#include <stdlib.h>
#include <stdio.h>
#include <limits>
#include <iostream>
using namespace std;




Wizard::Wizard(QWidget *parent) :
    QWizard(parent),
    ui(new Ui::Wizard)
{
    ui->setupUi(this);

    advance = 0;

    voronoidata.cells = NULL;
    voronoidata.vertices = NULL;
    voronoidata.cellsno = 0;
    voronoidata.validcellsno = 0;
    voronoidata.verticesno = 0;
    voronoidata.trajectoryindex = NULL;

    //animationtimer = new QTimer(this);
    connect(&animationtimer, SIGNAL(timeout()), this, SLOT(timer_update()));
}

Wizard::~Wizard()
{
    delete ui;
    //delete animationtimer;
}


//--------------------------------------------------------------
//--------------------- Load Page ------------------------------
//--------------------------------------------------------------

bool Wizard::readTo(QTextStream& stream, QString str)
{
    QChar c;
    int i = 0;
    int max = str.length();
    do
    {
        stream >> c;
        if(c == str.at(i)) i++;
        else i = 0;
        if(i >= max) return true;
    }while(!stream.atEnd());
    return false;
}


bool Wizard::loadComsolData(QTextStream& stream, QString* error)
{
    int nodes = 0;
    int expre = 0;
    if(!readTo(stream, "% Nodes:")) {*error = QString("\"% Nodes:\" tag missing"); return false;}
    stream >> nodes;
    if(nodes <= 0) {*error = QString("Nodes number incorrect"); return false;}
    if(!readTo(stream, "% Expressions:")) {*error = QString("\"% Expressions:\" tag missing"); return false;}
    stream >> expre;
    if(expre < 6 || (expre % 3 != 0)) {*error = QString("Expressions number incorrect"); return false;}
    expre /= 3;

    if(!readTo(stream, "% x")) {*error = QString("\"% x\" tag missing"); return false;}

    double t[3];
    for(int e = 0; e < expre; e++)
    {
        for(int d = 0; d < 3; d++)
        {
            if(!readTo(stream, "@ t=")) {*error = QString("\"@ t=\" tag missing - node %1, coordinate %2").arg(e).arg(d); return false;}
            stream >> t[d];
            if(stream.atEnd()) {*error = QString("Unexpected end of file"); return false;}
        }
        if(t[0] != t[1] || t[1] != t[2]) {*error = QString("Time indicators differ - node %1").arg(e); return false;}
        if(e == 0 && t[0] > 0.0) {*error = QString("Nonzero start time"); return false;}
    }
    double step = t[0] / (expre - 1);
    if(step <= 0.0) {*error = QString("Incorrect step"); return false;}

    trajectories.stepTime = step;
    trajectories.maxSteps = expre;
    trajectories.setNumber(nodes);
    for(int n = 0; n < nodes; n++)
    {
        int nn;
        stream >> nn;
        if(nn != n+1) {*error = QString("Index mismatch at line %1, index %2").arg(n+1).arg(nn); return false;}

        trajectories.setNumber(n, expre);
        for(int e = 0; e < expre; e++)
        {
            for(int d = 0; d < 3; d++)
            {
                stream >> t[d];
                if(stream.atEnd()) {*error = QString("Unexpected end of coordinates").arg(expre); return false;}
                if(stream.status() != 0 || !(t[d] == t[d]))
                {
                    trajectories.resetNumber(n, e);
                    d = 3; e = expre;
                    readTo(stream, "\n");
                }
            }
            trajectories.setCoordinates(n, e, t);
        }
    }
    return true;
}

bool Wizard::findCrossSection(void)
{
    double a[3];
    double b[3];
    double c[3];
    c[0] = 0.0; c[1] = 0.0; c[2] = 0.0;
    unsigned int k;
    unsigned int trn = trajectories.getNumber();
    for(k = 0; k < trn; k++)
    {
        if(trajectories.getCoordinates(k, 1, a) && trajectories.getCoordinates(k, 2, b))
        {
            c[0] += (b[0]-a[0]);
            c[1] += (b[1]-a[1]);
            c[2] += (b[2]-a[2]);
        }
    }
    double cc = c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
    if(cc < 0.000001) return false;
    cc = sqrt(cc);

    c[0] /= cc;
    c[1] /= cc;
    c[2] /= cc;

    bool first = true;
    double d = 0.0;
    for(k = 0; k < trn; k++)
    {
        if(trajectories.getCoordinates(k, 1, a))
        {
            if(first)
            {
                d = a[0]*c[0] + a[1]*c[1] + a[2]*c[2];
                first = false;
            }
            else
            {
                double dt = a[0]*c[0] + a[1]*c[1] + a[2]*c[2];
                if(dt > d) d = dt;
            }
        }
    }

    if(first) return false;

    ui->crossSectionXSpinBox->setValue(c[0]);
    ui->crossSectionYSpinBox->setValue(c[1]);
    ui->crossSectionZSpinBox->setValue(c[2]);
    ui->crossSectionSpinBox->setValue(d);

    if(trajectories.setInlet(c, d))
    {
        advance = 2;
        trajectories.computeProjection();
    }
    else advance = 1;
    return true;
}

QString Wizard::validateTrajectories(void)
{
    QString ret;
    double min_dist = ui->minimumDistanceSpinBox->value();

    advance = 0;
    if(trajectories.getNumber() <= 0) return QString("Empty");
    unsigned int trn = trajectories.getNumber();
    ret.append(QString("Number of trajectories: %1\n").arg(trn));
    unsigned int k, i;

    unsigned int min;
    unsigned int max;
    unsigned int all = 0;
    for(k = 0; k < trn; k++)
    {
        unsigned int non = trajectories.getNumber(k);
        if(k == 0)
        {
            min = non;
            max = non;
        }


        if(non < min) min = non;
        if(max < non) max = non;
        all += non;
    }

    ret.append(QString("Number of nodes: %1\n").arg(all));
    ret.append(QString("Number of nodes per trajectory: %1 : %2\n").arg(min).arg(max));
    ret.append(QString("Time step: %1\n").arg(trajectories.stepTime));
    ret.append(QString("Range X: %1 : %2\n").arg(trajectories.getRange(0, 0)).arg(trajectories.getRange(0, 1)));
    ret.append(QString("Range Y: %1 : %2\n").arg(trajectories.getRange(1, 0)).arg(trajectories.getRange(1, 1)));
    ret.append(QString("Range Z: %1 : %2\n").arg(trajectories.getRange(2, 0)).arg(trajectories.getRange(2, 1)));

    bool *dups = new bool[trn];
    for(k = 0; k < trn; k++) dups[k] = false;

    for(k = 0; k < trn-1; k++)
    {
        for(i = k+1; i < trn; i++)
        {
            double a[3];
            double b[3];
            if(trajectories.getCoordinates(k, 0, a) && trajectories.getCoordinates(i, 0, b))
            {
                double dist = a[0]-b[0];
                dist *= dist;
                double d = a[1]-b[1];
                dist += (d * d);
                d = a[2]-b[2];
                dist += (d * d);
                dist = sqrt(dist);

                if(dist < min_dist)
                    dups[i] = true;
            }
        }
    }
    int dup = 0;
    for(k = 0; k < trn; k++) if(dups[k]) dup++;
    ret.append(QString("\nDuplicates: %1\n").arg(dup));


    if(dup > 0) ui->removeDupsButton->setEnabled(true);
    else ui->removeDupsButton->setEnabled(false);

    for(k = 0; k < trn; k++) dups[k] = false;

    for(k = 0; k < trn; k++)
    {
        unsigned int imax = trajectories.getNumber(k);
        for(i = 1; i < imax; i++)
        {
            double a[3];
            double b[3];
            if(trajectories.getCoordinates(k, i, a) && trajectories.getCoordinates(k, i-1, b))
            {
                double dist = a[0]-b[0];
                dist *= dist;
                double d = a[1]-b[1];
                dist += (d * d);
                d = a[2]-b[2];
                dist += (d * d);
                dist = sqrt(dist);

                if(dist < min_dist)
                    dups[k] = true;
            }
        }
    }

    dup = 0;
    for(k = 0; k < trn; k++) if(dups[k]) dup++;
    ret.append(QString("Stuck nodes: %1\n").arg(dup));

    if(dup > 0) ui->removeStuckButton->setEnabled(true);
    else ui->removeStuckButton->setEnabled(false);

    ui->timeStepSpinBox->setValue(trajectories.stepTime);
    ui->render3DA->setTrajectories(&trajectories);
    ui->renderTes->setTrajectories(&trajectories);
    ui->renderVerify->setTrajectories(&trajectories);

    ui->viewStepSlider->setValue(0);
    ui->viewStepSpinBox->setValue(0);
    ui->viewStepSlider->setMaximum(trajectories.maxSteps);
    ui->viewStepSpinBox->setMaximum(trajectories.maxSteps);

    ui->render3DA->setCurrentStep(0);

    delete[] dups;

    advance = 1;
    return ret;
}

void Wizard::on_removeStuckButton_clicked()
{
    unsigned int k, i;
    if(trajectories.getNumber() <= 0) return;
    unsigned int trn = trajectories.getNumber();
    bool *dups = new bool[trn];
    for(k = 0; k < trn; k++) dups[k] = false;
    double min_dist = ui->minimumDistanceSpinBox->value();

    for(k = 0; k < trn; k++)
    {
        unsigned int imax = trajectories.getNumber(k);
        for(i = 1; i < imax; i++)
        {
            double a[3];
            double b[3];
            if(trajectories.getCoordinates(k, i, a) && trajectories.getCoordinates(k, i-1, b))
            {
                double dist = a[0]-b[0];
                dist *= dist;
                double d = a[1]-b[1];
                dist += (d * d);
                d = a[2]-b[2];
                dist += (d * d);
                dist = sqrt(dist);

                if(dist < min_dist)
                    dups[k] = true;
            }
        }
    }

    for(k = trn-1; k < trn; k--) if(dups[k]) trajectories.removeTrajectory(k);
    delete[] dups;

    ui->plainTextEdit->setPlainText(validateTrajectories());
}

void Wizard::on_removeDupsButton_clicked()
{
    unsigned int k, i;
    if(trajectories.getNumber() <= 0) return;
    unsigned int trn = trajectories.getNumber();
    bool *dups = new bool[trn];
    for(k = 0; k < trn; k++) dups[k] = false;
    double min_dist = ui->minimumDistanceSpinBox->value();

    for(k = 0; k < trn-1; k++)
    {
        for(i = k+1; i < trn; i++)
        {
            double a[3];
            double b[3];
            if(trajectories.getCoordinates(k, 0, a) && trajectories.getCoordinates(i, 0, b))
            {
                double dist = a[0]-b[0];
                dist *= dist;
                double d = a[1]-b[1];
                dist += (d * d);
                d = a[2]-b[2];
                dist += (d * d);

                if(dist < min_dist) dups[i] = true;
            }
        }
    }
    for(k = trn-1; k > 0; k--) if(dups[k]) trajectories.removeTrajectory(k);
    delete[] dups;

    ui->plainTextEdit->setPlainText(validateTrajectories());
}

void Wizard::on_loadButton_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this,
                                                    tr("Load particle flow file"),
                                                    NULL,//domyslna nazwa pliku
                                                    tr("Flow of particles (*.txt) (*.txt);;All Files (*)"));
    if (!fileName.isEmpty())
    {
        QFile file(fileName);
        if(file.open(QIODevice::ReadOnly | QIODevice::Text))
        {

            QTextStream in(&file);
            QString err;
            if(loadComsolData(in, &err))
            {
                QString ret = validateTrajectories();
                double c[3];
                double d;
                if(findCrossSection())
                {
                    trajectories.getInlet(c, &d);
                    ret.append(QString("Inlet: [%1 %2 %3] %4\n").arg(c[0]).arg(c[1]).arg(c[2]).arg(d));
                }
                else ret.append(QString("Inlet estimation error\n"));
                ui->plainTextEdit->setPlainText(ret);
            }
            else
            {
                ui->plainTextEdit->setPlainText(err);
                advance = 0;
            }
        }
    }
}

void Wizard::on_timeStepSpinBox_valueChanged(double arg1)
{
    trajectories.stepTime = arg1;
    ui->plainTextEdit->setPlainText(validateTrajectories());
}


void Wizard::on_minimumDistanceSpinBox_valueChanged(double arg1)
{
    ui->plainTextEdit->setPlainText(validateTrajectories());
}


//--------------------------------------------------------------
//--------------------- Cross Page -----------------------------
//--------------------------------------------------------------


void Wizard::on_showTrajectoriesCheck_clicked(bool checked)
{
    ui->render3DA->setShowTrajectories(checked);
}

void Wizard::on_viewStepSpinBox_valueChanged(int arg1)
{
    ui->render3DA->setCurrentStep(arg1);
}


bool Wizard::validateCurrentPage()
{
    if(currentId() > advance)
    {
        //button(QWizard::NextButton)->setEnabled(false);
        return false;
    }
    else
    {
        //button(QWizard::NextButton)->setEnabled(false);
        return true;
    }
}

void Wizard::on_crossSectionXSpinBox_valueChanged(double arg1)
{
    on_crossSectionSpinBox_valueChanged(arg1);
}

void Wizard::on_crossSectionYSpinBox_valueChanged(double arg1)
{
    on_crossSectionSpinBox_valueChanged(arg1);
}

void Wizard::on_crossSectionZSpinBox_valueChanged(double arg1)
{
    on_crossSectionSpinBox_valueChanged(arg1);
}

void Wizard::on_crossSectionSpinBox_valueChanged(double arg1)
{
    double c[3];
    double d;
    c[0] = ui->crossSectionXSpinBox->value();
    c[1] = ui->crossSectionYSpinBox->value();
    c[2] = ui->crossSectionZSpinBox->value();
    d = ui->crossSectionSpinBox->value();
    if(trajectories.setInlet(c, d)) advance = 2;
    else advance = 1;

    ui->render3DA->update();
}




//--------------------------------------------------------------
//--------------------- Tessellation ---------------------------
//--------------------------------------------------------------

void Wizard::assignPoints(void)
{
    double params[params_no];
    int i;
    double xl[2], xr[2];
    double Area = 0;
    double AreaC = 0;
    double AreaB = 0;
    for(i = 0; i < 4; i++) params[i] = -1.0;

    double Xc = voronoidata.circle[0];
    double Yc = voronoidata.circle[1];
    double Rc = voronoidata.circle[2];

    if(voronoidata.validcellsno > 0 && voronoidata.validcellsno <= voronoidata.cellsno)
    {
        for(i = 0; i < voronoidata.validcellsno; i++)
        {
            Area = 0;
            int brzeg = 0;
            int maxk = voronoidata.cells[i][0];
            int il, ir;
            il = voronoidata.cells[i][maxk];

            for(int k = 1; k <= maxk; k++)
            {
                ir = voronoidata.cells[i][k];
                if(ir > 0 && il > 0)
                {
                    xr[0] = voronoidata.vertices[ir*2];
                    xr[1] = voronoidata.vertices[ir*2+1];
                    xl[0] = voronoidata.vertices[il*2];
                    xl[1] = voronoidata.vertices[il*2+1];
                    Area += (xl[0]*xr[1] - xr[0]*xl[1]);

                    if((xr[0]-Xc)*(xr[0]-Xc) + (xr[1]-Yc)*(xr[1]-Yc) > Rc*Rc) brzeg++;
                }
                il = ir;
            }
            trajectories.getTrajectoryParams(voronoidata.trajectoryindex[i], params);
            if(brzeg) params[params_area] = -fabs(Area);
            else params[params_area] = fabs(Area);

            params[params_random] = 0.0;
            params[params_repetition] = -1;
            params[params_volume] = -1;
            trajectories.setTrajectoryParams(voronoidata.trajectoryindex[i], params);
        }
        // Korekta powierzchni komorek brzegowych
        for(i = 0; i < voronoidata.validcellsno; i++)
        {
            double a = 0;
            if(trajectories.getTrajectoryParams(voronoidata.trajectoryindex[i], params)) a = params[params_area];
            Area += fabs(a);
            if(a < 0) AreaB += fabs(a);
        }
        AreaC = M_PI * Rc * Rc;
        if(Area > AreaC && Area-AreaC < AreaB && AreaB > 0.0)
        {
            AreaB = (AreaB-Area+AreaC)/AreaB;
            for(i = 0; i < voronoidata.validcellsno; i++)
            {
                double a = 0;
                if(trajectories.getTrajectoryParams(voronoidata.trajectoryindex[i], params)) a = params[params_area];
                if(a < 0)
                {
                    trajectories.getTrajectoryParams(voronoidata.trajectoryindex[i], params);
                    params[params_area] = AreaB * fabs(a);
                    trajectories.setTrajectoryParams(voronoidata.trajectoryindex[i], params);
                }
            }
        }
        for(i = 0; i < voronoidata.validcellsno; i++)
        {
            if(trajectories.getTrajectoryParams(voronoidata.trajectoryindex[i], params))
            {
                trajectories.getTrajectoryParams(voronoidata.trajectoryindex[i], params);
                params[params_area] = fabs(params[params_area]);
                trajectories.setTrajectoryParams(voronoidata.trajectoryindex[i], params);
            }
        }
        advance = 3;
    }
}


void Wizard::tesselate_yes(void)
{
    double rx[2];
    double ry[2];

    trajectories.get2DRanges(rx, ry);
    if(rx[1] <= rx[0] || ry[1] <= ry[0]) return;

    //point = vschPoints;
    //exitcode = qh_new_qhull(dim, vm-vs, point, ismalloc, (char*)("qhull n"), stdout, stderr);

    int liczba_czastek = trajectories.getNumber();
    int pointsincircle = 0;
    int pointsincirclo = 0;
    double r, rr;
    r = ui->radiusSpinBox->value();
    rr = r*r;
    double cx = ui->xSpinBox->value();
    double cy = ui->ySpinBox->value();
    for(int part = 0; part < liczba_czastek; part++)
    {
        double x[2];

        if(trajectories.getProjectedInletPoint(part, x))
        {
            double x0 = x[0] - cx;
            double x1 = x[1] - cy;
            if(x0 * x0 + x1 * x1 < rr)
            {
                pointsincircle++;
            }
        }
    }


    for(int part = 0; part < liczba_czastek; part++)
    {
        double x[2];

        if(trajectories.getProjectedInletPoint(part, x))
        {
            double x0 = x[0] - cx;
            double x1 = x[1] - cy;
            double rxy = x0 * x0 + x1 * x1;

            if(rxy < rr)
            {
                //double rxys = sqrt(rxy);
                if(rxy > rr/1000)
                {
                    pointsincirclo++;
                }
            }
        }
    }


    if(pointsincircle <= 0) return;

    coordT *pointstable = new coordT[(pointsincirclo+pointsincircle)*2];
    coordT *ptr = pointstable;
    for(int part = 0; part < liczba_czastek; part++)
    {
        double x[2];
        if(trajectories.getProjectedInletPoint(part, x))
        {
            double x0 = x[0] - cx;
            double x1 = x[1] - cy;
            if(x0 * x0 + x1 * x1 < rr)
            {
                *ptr = x[0];
                ptr++;
                *ptr = x[1];
                ptr++;
            }
        }
    }

    for(int part = 0; part < liczba_czastek; part++)
    {
        double x[2];
        if(trajectories.getProjectedInletPoint(part, x))
        {
            double x0 = x[0] - cx;
            double x1 = x[1] - cy;
            double rxy = x0 * x0 + x1 * x1;

            if(rxy < rr)
            {
                double rxys = sqrt(rxy);
                if(rxy > rr/1000)
                {
                    rxys = (2*r - rxys) / rxys;
                    x[0] = rxys*x0 + cx;
                    x[1] = rxys*x1 + cy;

                    *ptr = x[0];
                    ptr++;
                    *ptr = x[1];
                    ptr++;
                }
            }
        }
    }

    boolT ismalloc = False;
    int curlong, totlong;

    int exitcode;
    exitcode = qh_new_qhull(2, pointsincircle+pointsincirclo, pointstable, ismalloc, (char*)("qhull v Qbb o"), NULL, stderr);
    //exitcode = qh_new_qhull(2, pointsincircle+pointsincirclo, pointstable, ismalloc, (char*)("qhull v Qbb o"), stdout, stderr);

    if (!exitcode)
    {
        voronoidata.validcellsno = pointsincircle;
        collect_voronoi(&voronoidata);
        voronoidata.circle[0] = cx;
        voronoidata.circle[1] = cy;
        voronoidata.circle[2] = r;
        int p;
        for(p = pointsincircle; p < voronoidata.cellsno; p++)
            voronoidata.trajectoryindex[p]  = -1;
        p = 0;

        for(int part = 0; part < liczba_czastek; part++)
        {
            double x[2];
            if(trajectories.getProjectedInletPoint(part, x))
            {
                double x0 = x[0] - cx;
                double x1 = x[1] - cy;
                if(x0 * x0 + x1 * x1 < rr)
                {
                    voronoidata.trajectoryindex[p] = part;
                    p++;
                }
            }
        }
    }
    else
    {
        voronoidata.validcellsno = 0;
    }

    qh_freeqhull(!qh_ALL);
    qh_memfreeshort (&curlong, &totlong);
    if (curlong || totlong)	fprintf(stderr, "qhull warning: did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);

    delete[] pointstable;

    assignPoints();

    ui->renderTes->setVoronoi(&voronoidata);
}

void Wizard::tesselate_no(void)
{
    double params[params_no];
    int liczba_czastek = trajectories.getNumber();
    int pointsincircle = 0;
    double r, rr;
    r = ui->radiusSpinBox->value();
    rr = r*r;
    double cx = ui->xSpinBox->value();
    double cy = ui->ySpinBox->value();
    double average_area;
    for(int part = 0; part < liczba_czastek; part++)
    {
        double x[2];
        if(trajectories.getProjectedInletPoint(part, x))
        {
            double x0 = x[0] - cx;
            double x1 = x[1] - cy;
            if(x0 * x0 + x1 * x1 < rr)
            {
                pointsincircle++;
            }
        }
    }
    if(pointsincircle > 0)
    {
        average_area = M_PI * r * r / pointsincircle;
    }
    for(int part = 0; part < liczba_czastek; part++)
    {
        double x[2];
        if(trajectories.getProjectedInletPoint(part, x))
        {
            double x0 = x[0] - cx;
            double x1 = x[1] - cy;
            if(x0 * x0 + x1 * x1 < rr)
            {
                trajectories.getTrajectoryParams(part, params);
                params[params_area] = average_area;
                trajectories.setTrajectoryParams(part, params);
            }
        }
    }
    ui->renderTes->setVoronoi(NULL);

}


void Wizard::computeParams(void)
{
    double rr;
    double params[params_no];
    double r = ui->radiusSpinBox->value();
    double cx = ui->xSpinBox->value();
    double cy = ui->ySpinBox->value();
    int liczba_czastek = trajectories.getNumber();
    double original_timestep = trajectories.stepTime;//ui->timeStepSpinBox->value();
    rr = r*r;
    for(int part = 0; part < liczba_czastek; part++)
    {
        double x[2];
        if(trajectories.getProjectedInletPoint(part, x))
        {
            double x0 = x[0] - cx;
            double x1 = x[1] - cy;
            if(x0 * x0 + x1 * x1 < rr)
            {
                trajectories.getTrajectoryParams(part, params);
                double length, volume, repetition;

                switch(ui->whatEqualCombo->currentIndex())
                {
                case 0:
                    volume = ui->valueEqualSpinBox->value();
                    length = volume/params[params_tilt]/params[params_area];
                    repetition = original_timestep * length/params[params_distance];
                    break;
                case 1:
                    length = ui->valueEqualSpinBox->value();
                    volume = length * params[params_area] * params[params_tilt];
                    repetition = original_timestep * length/params[params_distance];
                    break;
                case 2:
                    repetition = ui->valueEqualSpinBox->value();
                    length = repetition * params[params_distance] / original_timestep;
                    volume = length * params[params_area] * params[params_tilt];
                    break;
                default:
                    repetition = -1.0;
                    length = -1.0;
                    volume = -1.0;
                }
                params[params_repetition] = repetition;
                params[params_volume] = volume;
                params[params_random] = 0.0;
                trajectories.setTrajectoryParams(part, params);
            }
        }
    }

    on_randomizeCheckBox_clicked(ui->randomizeCheckBox->checkState() == Qt::Checked);
}


void Wizard::on_addInletButton_clicked()
{
    if(ui->tesselateCheckBox->checkState() == Qt::Checked) tesselate_yes();
    else tesselate_no();

    computeParams();

    ui->renderTes->update();
}

void Wizard::on_clearButton_clicked()
{
    int i;
    double params[params_no];

    int liczba_czastek = trajectories.getNumber();
    for(int i = 0; i < liczba_czastek; i++)
    {
        trajectories.getTrajectoryParams(i, params);
        params[params_area] = -1.0;
        trajectories.setTrajectoryParams(i, params);
    }

    if(voronoidata.trajectoryindex != NULL) free(voronoidata.trajectoryindex);
    if(voronoidata.vertices != NULL) free(voronoidata.vertices);
    if(voronoidata.cells != NULL)
    {
        for(i = 0; i < voronoidata.cellsno; i++) free(voronoidata.cells[i]);
        free(voronoidata.cells);
    }

    voronoidata.cells = NULL;
    voronoidata.vertices = NULL;
    voronoidata.cellsno = 0;
    voronoidata.validcellsno = 0;
    voronoidata.verticesno = 0;
    voronoidata.trajectoryindex = NULL;

    advance = 2;

    ui->renderTes->update();
}

void Wizard::on_radiusSpinBox_valueChanged(double arg1)
{
    double average_distance =0.0;
    double average_tilt =0.0;
    double params[params_no];
    int liczba_czastek = trajectories.getNumber();
    int pointsincircle = 0;
    double r, rr;
    r = ui->radiusSpinBox->value();
    rr = r*r;
    double cx = ui->xSpinBox->value();
    double cy = ui->ySpinBox->value();
    for(int part = 0; part < liczba_czastek; part++)
    {
        double x[2];
        if(trajectories.getProjectedInletPoint(part, x))
        {
            double x0 = x[0] - cx;
            double x1 = x[1] - cy;
            if(x0 * x0 + x1 * x1 < rr)
            {
                pointsincircle++;
                params[params_distance] = 0.0;
                trajectories.getTrajectoryParams(part, params);
                average_distance += params[params_distance];
                average_tilt += params[params_tilt];
            }
        }
    }
    if(pointsincircle > 0)
    {
        average_tilt /= pointsincircle;
        average_distance /= pointsincircle;
        //double circle_area = M_PI * r * r;
        double average_area = M_PI * r * r / pointsincircle;
        double average_length = sqrt(average_area)/average_tilt;
        double average_volume = average_length*average_area;
        double average_repetition = average_length / average_distance;
        average_repetition *= ui->timeStepSpinBox->value();

        switch(ui->whatEqualCombo->currentIndex())
        {
            case 0: ui->valueEqualSpinBox->setValue(average_volume); break;
            case 1: ui->valueEqualSpinBox->setValue(average_length); break;
            case 2: ui->valueEqualSpinBox->setValue(average_repetition); break;
            default: ui->valueEqualSpinBox->setValue(std::numeric_limits<double>::quiet_NaN());
        }
    }
    else ui->valueEqualSpinBox->setValue(std::numeric_limits<double>::quiet_NaN());
}

void Wizard::on_xSpinBox_valueChanged(double arg1)
{
    on_radiusSpinBox_valueChanged(arg1);
}

void Wizard::on_ySpinBox_valueChanged(double arg1)
{
    on_radiusSpinBox_valueChanged(arg1);
}

void Wizard::on_whatEqualCombo_currentIndexChanged(int index)
{
    on_radiusSpinBox_valueChanged(0.0);
}



//--------------------------------------------------------------
//--------------------- Tessellation ---------------------------
//--------------------------------------------------------------

void Wizard::timer_update()
{
    int value = ui->animateSlider->value();
    if(!ui->animateSlider->isSliderDown())
    {
        if(value >= 0) ui->animateSlider->setValue(value - 1);
        else ui->animateSlider->setValue(value + 1);
    }

    double atime = ui->animationTimeSpinBox->value();
    double stime = ui->timeStepSpinBox->value();
    double dtime = value*value*value;
    dtime /= 2000.0;
    dtime *= stime;
    ui->animationTimeSpinBox->setValue(atime+dtime);
}

void Wizard::on_animateSlider_valueChanged(int value)
{
    if(value != 0)
    {
        if(!animationtimer.isActive()) animationtimer.start(66);
    }
    else animationtimer.stop();
}

void Wizard::on_animationTimeSpinBox_valueChanged(double arg1)
{
    ui->renderVerify->setCurrentTime(arg1);
    ui->renderVerify->update();
}

void Wizard::on_randomizeCheckBox_clicked(bool checked)
{
    srand(20131029);
    int liczba_czastek = trajectories.getNumber();

    for(int part = 0; part < liczba_czastek; part++)
    {
        double params[params_no];
        trajectories.getTrajectoryParams(part, params);
        if(params[params_area] > 0.0 && params[params_repetition] > 0)
        {
            if(checked) params[params_random] = (double)rand()*params[params_repetition]/(double)RAND_MAX;
            else params[params_random] = 0.0;
            trajectories.setTrajectoryParams(part, params);
        }
    }
    ui->renderVerify->update();
}

bool Wizard::saveTrajectories(QString fileName)
{
    QFile file(fileName);
    if(!file.open(QIODevice::WriteOnly| QIODevice::Text)) return false;
    {
        QTextStream out(&file);
        out.setRealNumberPrecision(12);

        int notajects = trajectories.getNumber();

        int count = 0;
        for(int part = 0; part < notajects; part++)
        {
            double params[params_no];
            trajectories.getTrajectoryParams(part, params);
            //if(params[params_area] > 0.0 && params[params_repetition] > 0)
            {
                count++;
            }
        }


        out << "number_of_trajectories " << count <<"\n";
        out << "node_repetition_time " << trajectories.stepTime <<"\n";
        count = 0;

        for(int part = 0; part < notajects; part++)
        {
            double params[params_no];
            trajectories.getTrajectoryParams(part, params);
            //if(params[params_area] > 0.0 && params[params_repetition] > 0)
            {
                int nonodes = trajectories.getNumber(part);
                out << "trajectory " << count << "\n";
                out << "\tparticle_volume " << params[params_volume] << "\n";
                out << "\tparticle_repetition_time " << params[params_repetition] << "\n";
                out << "\trandomization_time " << params[params_random] << "\n";
                out << "\ttrajectory_nodes " << nonodes << "\n";

                out << "\t";
                for(int k = 0; k < nonodes; k++)
                {
                    double x[3];
                    char str[256];
                    trajectories.getCoordinates(part, k, x);
                    /*
                    sprintf(str, "(%.12f %.12f %.12f) ", x[0], x[1], x[2]);
                    out << str;
                    */
                    out << "(" << double(x[0]);
                    out << " " << double(x[1]);
                    out << " " << double(x[2]);
                    out << ") ";
                }
                out << "\n";
                count++;
            }
        }
        out.flush();
    }
    file.flush();
    file.close();

    return true;
}



void Wizard::on_saveButton_clicked()
{
    QString fileName = QFileDialog::getSaveFileName(this,
                                    tr("Save trajectories"),
                                    "trajektorie.txt",
                                    tr("Trajektorie (*.txt) (*.txt);;All Files (*)"));
    if (!fileName.isEmpty())
    {
        if(saveTrajectories(fileName)) advance = 4;
    }
}
