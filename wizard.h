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

#ifndef WIZARD_H
#define WIZARD_H

#include <QtGlobal>
#include <QWizard>
#include <QFile>
#include <QFileDialog>
#include <QTextStream>
#include <QMessageBox>
#include <QTimer>

#include "trajectories.h"

extern "C"
{
#include "vvoronoi.h"
}

namespace Ui {
class Wizard;
}

class Wizard : public QWizard
{
    Q_OBJECT
    
public:
    explicit Wizard(QWidget *parent = 0);
    ~Wizard();
    
    bool validateCurrentPage();

private slots:
    void timer_update();

    void on_loadButton_clicked();

    void on_removeDupsButton_clicked();

    void on_showTrajectoriesCheck_clicked(bool checked);

    void on_viewStepSpinBox_valueChanged(int arg1);

    void on_timeStepSpinBox_valueChanged(double arg1);

    void on_removeStuckButton_clicked();

    void on_minimumDistanceSpinBox_valueChanged(double arg1);

    void on_crossSectionXSpinBox_valueChanged(double arg1);

    void on_crossSectionYSpinBox_valueChanged(double arg1);

    void on_crossSectionZSpinBox_valueChanged(double arg1);

    void on_crossSectionSpinBox_valueChanged(double arg1);

    void on_addInletButton_clicked();

    void on_clearButton_clicked();

    void on_radiusSpinBox_valueChanged(double arg1);

    void on_xSpinBox_valueChanged(double arg1);

    void on_ySpinBox_valueChanged(double arg1);

    void on_whatEqualCombo_currentIndexChanged(int index);

    void on_animateSlider_valueChanged(int value);

    void on_animationTimeSpinBox_valueChanged(double arg1);

    void on_randomizeCheckBox_clicked(bool checked);

    void on_saveButton_clicked();

private:
    int advance;

    Ui::Wizard *ui;
    Trajectories trajectories;

    bool readTo(QTextStream& stream, QString str);
    bool loadComsolData(QTextStream &stream, QString *error);
    QString validateTrajectories(void);
    bool findCrossSection(void);
    void assignPoints(void);
    bool saveTrajectories(QString fileName);

    void tesselate_no(void);
    void tesselate_yes(void);
    void computeParams(void);

    QTimer animationtimer;

    VoronoiData voronoidata;
};

#endif // WIZARD_H
