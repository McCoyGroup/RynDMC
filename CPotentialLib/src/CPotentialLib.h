#include "Python.h"
#include "PyExtLib.h"

static PyObject *CPotentialLib_callPot
    ( PyObject *, PyObject * );

static PyObject *CPotentialLib_callPotVec
    ( PyObject *, PyObject * );

static PyObject *CPotentialLib_callPotMPI
    ( PyObject *, PyObject * );

static PyObject *CPotentialLib_getMPIInfo
    ( PyObject *, PyObject * );

static PyObject *CPotentialLib_closeMPI
    ( PyObject *, PyObject * );
