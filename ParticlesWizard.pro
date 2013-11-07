#-------------------------------------------------
#
# Project created by QtCreator 2013-10-24T11:44:46
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = ParticlesWizard
TEMPLATE = app


SOURCES += main.cpp\
    wizard.cpp \
    trajectories.cpp \
    view3d.cpp \
#    libqhull/userprintf_rbox.c \
#    libqhull/usermem.c \
#    libqhull/user.c \
#    libqhull/stat.c \
#    libqhull/rboxlib.c \
#    libqhull/random.c \
#    libqhull/qset.c \
#    libqhull/poly2.c \
#    libqhull/poly.c \
#    libqhull/merge.c \
#    libqhull/mem.c \
#    libqhull/libqhull.c \
#    libqhull/io.c \
#    libqhull/global.c \
#    libqhull/geom2.c \
#    libqhull/geom.c \
#    libqhull/userprintf.c \
    vvoronoi.c


HEADERS  += \
    wizard.h \
    trajectories.h \
    view3d.h \
#    libqhull/user.h \
#    libqhull/stat.h \
#    libqhull/random.h \
#    libqhull/qset.h \
#    libqhull/qhull_a.h \
#    libqhull/poly.h \
#    libqhull/merge.h \
#    libqhull/mem.h \
#    libqhull/libqhull.h \
#    libqhull/io.h \
#    libqhull/geom.h \
    vvoronoi.h


FORMS    += \
    wizard.ui

RESOURCES += \
    particleswizard.qrc

unix|win32: LIBS += -lqhull

