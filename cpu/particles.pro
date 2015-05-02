QT += core gui opengl

TARGET = particles
TEMPLATE = app

CONFIG += c++0x
QMAKE_CXXFLAGS += -std=c++0x

# If you add your own folders, add them to INCLUDEPATH and DEPENDPATH, e.g.
# INCLUDEPATH += folder1 folder2
# DEPENDPATH += folder1 folder2

INCLUDEPATH += src src/solver src/constraint glm
DEPENDPATH += src src/solver src/constraint glm

SOURCES += src/main.cpp \
    src/mainwindow.cpp \
    src/view.cpp \
    src/simulation.cpp \
    src/constraint/distanceconstraint.cpp \
    src/solver/lineareq.cpp \
    src/solver/matrix.cpp \
    src/solver/matrix.inl \
    src/solver/solver.cpp \
    src/solver/particle.cpp \
    src/constraint/totalshapeconstraint.cpp \
    src/constraint/boundaryconstraint.cpp \
    src/constraint/contactconstraint.cpp \
    src/constraint/totalfluidconstraint.cpp \
    src/constraint/rigidcontactconstraint.cpp \
    src/constraint/gasconstraint.cpp

HEADERS += src/mainwindow.h \
    src/view.h \
    src/simulation.h \
    src/particle.h \
    src/includes.h \
    src/constraint/distanceconstraint.h \
    src/solver/lineareq.h \
    src/solver/matrix.h \
    src/solver/solver.h \
    src/constraint/totalshapeconstraint.h \
    src/constraint/boundaryconstraint.h \
    src/constraint/contactconstraint.h \
    src/constraint/totalfluidconstraint.h \
    src/constraint/rigidcontactconstraint.h \
    src/constraint/gasconstraint.h

# UMFPACK
INCLUDEPATH += $$PWD/lib/umfpack/include
LIBS += -L$$PWD/lib/umfpack \
        -lumfpack \
        -lamd \
        #-lblas \
        #-lcerbla

FORMS += src/mainwindow.ui
