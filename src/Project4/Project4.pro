TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += C:/Armadillo/include
DEPENDPATH += C:/Armadillo/include
LIBS += -LC:/Armadillo/newblas/ -llibblas
LIBS += -LC:/Armadillo/newblas/ -lliblapack

SOURCES += main.cpp
