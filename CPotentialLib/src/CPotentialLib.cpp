#include "CPotentialLib.h"
#include <vector>
#include <string>

#ifdef USE_MPI
#include "mpi.h"
#endif

/********************** EXTRACT WALKER ATOM TYPES *******************************/
std::vector<std::string> _getAtomTypes( PyObject* atoms, Py_ssize_t num_atoms ) {

    std::vector<std::string> atom_types(num_atoms);
    for (int i = 0; i<num_atoms; i++) {
        PyObject* atom = PyList_GetItem(atoms, i);
        PyObject* pyStr = NULL;
        const char* atomStr = _GetString(atom, pyStr);
        std::string atomString = atomStr;
        atom_types[i] = atomString;
//        Py_XDECREF(atom);
        Py_XDECREF(pyStr);
    }

    return atom_types;
}

/******************** EXTRACT WALKER COORDINATE VALUES ****************************/
double *_GetDoubleDataBufferArray(Py_buffer *view) {
    return _GetDataBufferArray<double>(view);
}
double *_GetDoubleDataArray(PyObject *data) {
    return _GetDataArray<double>(data);
}

int ind3D(int i, int j, int k, int m, int l) {
    return _idx3D(i, j, k, 0, m, l);
}

std::vector< std::vector<double> > _getWalkerCoords(double* raw_data, int i, Py_ssize_t num_atoms) {
    std::vector< std::vector<double> > walker_coords(num_atoms, std::vector<double>(3));
    for (int j = 0; j<num_atoms; j++) {
        for (int k = 0; k<3; k++) {
            walker_coords[j][k] = raw_data[ind3D(i, j, k, num_atoms, 3)];
        }
    };
    return walker_coords;
}

/******************** CALL POINTERED POTENTIAL ****************************/
// I catch all the errors (I think) and just never let anyone know there was an issue
// this is something to be careful about when using this
double _callPtrPot(
        potFuncPtr potPtr,
        std::vector< std::vector<double> > walker_coords,
        std::vector<std::string> atoms
        ) {
    double pot;
    try {
        pot = (*potPtr)(walker_coords, atoms);
    } catch (int e) {
//        _printOutWalkerStuff(walker_coords);
        pot = 1.0e9;
    } catch (const char* e) {
//        _printOutWalkerStuff(walker_coords);
        pot = 1.0e9;
    } catch (...) {
//        _printOutWalkerStuff(walker_coords);
        pot = 1.0e9;
    }

    return pot;
};

/******************** CALL POTENTIAL PYTHON ***************************/
// First hook into _callPtrPot
// Just does it for one geometry
FUNCWITHARGS(CPotentialLib_callPot) {

    PyObject *ptr, *atoms, *coords;
    PARSEARGS("OOO", &ptr, &atoms, &coords);

    potFuncPtr potPtr = (potFuncPtr) PyCapsule_GetPointer(ptr, "potential");
    // asserts basically that ptr will be a pointer to a potential function
    // but if this doesn't hold you can get segfaults and stuff I bet

    // Assumes we get n atom type names
    Py_ssize_t num_atoms = PyObject_Length(atoms);
    std::vector<std::string> atom_types = _getAtomTypes(atoms, num_atoms);

    // Assumes number of walkers X number of atoms X 3
    double* raw_data = _GetDoubleDataArray(coords);
    CHECKNULL(raw_data);
    std::vector< std::vector<double> > walker_coords = _getWalkerCoords(raw_data, 0, num_atoms);
    double pot = _callPtrPot(potPtr, walker_coords, atom_types);

    PyObject *potVal = Py_BuildValue("f", pot);
    return potVal;

}


/*********************************** MPI Stuff ******************************************/
#ifdef USE_MPI

void _mpiInit(int* world_size, int* world_rank) {
    // Initialize MPI state
    int did_i_do_good_pops = 0;
    MPI_Initialized(&did_i_do_good_pops); // need to check if we called Init once already
    // printf("How'd I do? %d\n", did_i_do_good_pops);
    if (!did_i_do_good_pops){
        MPI_Init(NULL, NULL);
       };
    MPI_Comm_size(MPI_COMM_WORLD, world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, world_rank);
    // printf("This many things %d and I am here %d\n", *world_size, *world_rank);
}

void _mpiFinalize() {
    int did_i_do_bad_pops = 0;
    MPI_Finalized(&did_i_do_bad_pops); // need to check if we called Init once already
    if (!did_i_do_bad_pops){
        MPI_Finalize();
       };
}

FUNCWITHARGS(CPotentialLib_getMPIInfo) {

    PyObject *hello;
    int world_size, world_rank;
    _mpiInit(&world_size, &world_rank);
    hello = Py_BuildValue("(ii)", world_rank, world_size);
    return hello;

}

FUNCWITHARGS(CPotentialLib_closeMPI) {

    _mpiFinalize();
    Py_RETURN_NONE;

}

// what we really want to do here is a bit up in the air...
// to be able to handle walker propagation inside an MPI call requires knowing how to propagate the thing
// ... which requires calling into python for that or similar
std::vector<std::vector<double> > _callMPIPot(
        potFuncPtr potPtr,
        double* raw_data,
        std::vector<std::string> atoms,
        Py_ssize_t num_walkers,
        Py_ssize_t num_atoms,
        int world_size,
        int world_rank
        ) {

    // we're not gonna do this full out for now, but rather save that for later

    // create a buffer for the walkers to be fed into MPI
    double* walker_buf = (double*) malloc(num_atoms*3*sizeof(double));
    // Scatter data buffer to processors
    MPI_Scatter(
                raw_data,  // raw data buffer to chunk up
                3*num_atoms, // three coordinates per atom per num_atoms per walker
                MPI_DOUBLE, // coordinates stored as doubles
                walker_buf, // raw array to write into
                3*num_atoms, // single energy
                MPI_DOUBLE, // energy returned as doubles
                0, // root caller
                MPI_COMM_WORLD // communicator handle
                );

    std::vector< std::vector<double> > walker_coords = _getWalkerCoords(walker_buf, 0, num_atoms);
    double* pot_buf = (double*) malloc(world_size*sizeof(double)); // receive buffer
    double pot = _callPtrPot(potPtr, walker_coords, atoms);

    MPI_Gather(pot&, 1, MPI_DOUBLE, pot_buf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    free(walker_buf);

    // convert double* to std::vector<double>
    std::vector<double> potVals(world_size);
    if( world_rank == 0 ) {
        for (size_t n = 0; n < world_size; n++) {
            potVals[n] = pot_buf[n]
        }
    }
    free(pot_buf);

    return potVals;
}

FUNCWITHARGS(CPotentialLib_callPotMPI) {

    PyObject *ptr, *atoms, *coords;
    PARSEARGS("OOO", &ptr, &atoms, &coords)

    potFuncPtr potPtr = (potFuncPtr) PyCapsule_GetPointer(ptr, "potential");

    // Assumes we get n atom type names
    Py_ssize_t num_atoms = PyObject_Length(atoms);
    std::vector<std::string> atom_types = _getAtomTypes(atoms, num_atoms);

    // Assumes number of walkers X number of atoms X 3
    Py_ssize_t num_walkers = PyObject_Length(coords);
    double* raw_data = _GetDoubleDataArray(coords);
    CHECKNULL(raw_data);

    int world_size, world_rank;
    _mpiInit(&world_size, &world_rank);
    std::vector<double> pot_vals = _mpiCallPot(potPtr, raw_data, atom_types,
                                               num_walkers, num_atoms,
                                               world_size, world_rank
                                               );

    int dims[1] = {num_walkers};
    PyObject* pots = _CreateArray(1, dims, "zeros", "float64");
    double* pot_data = _GetDoubleDataArray(pots);
    memcpy(pot_vals.data(), pot_data, sizeof(double) * world_size);

    return pots;

}

#else

FUNCWITHARGS(CPotentialLib_callPotMPI) {
    PyErr_SetString(PyExc_NotImplementedError, "CPotentialLib was compiled with UseMPI turned off");
    return NULL;
}

FUNCWITHARGS(CPotentialLib_getMPIInfo) {
    PyErr_SetString(PyExc_NotImplementedError, "CPotentialLib was compiled with UseMPI turned off");
    return NULL;
}

FUNCWITHARGS(CPotentialLib_closeMPI) {
    PyErr_SetString(PyExc_NotImplementedError, "CPotentialLib was compiled with UseMPI turned off");
    return NULL;
}

#endif

FUNCWITHARGS(CPotentialLib_callPotVec) {
    // vector version of callPot

    PyObject *ptr, *atoms, *coords;
    PARSEARGS("OOO", &ptr, &atoms, &coords)

    potFuncPtr potPtr = (potFuncPtr) PyCapsule_GetPointer(ptr, "potential");

    // Assumes we get n atom type names
    Py_ssize_t num_atoms = PyObject_Length(atoms);
    std::vector<std::string> atom_types = _getAtomTypes(atoms, num_atoms);

    // Assumes number of walkers X number of atoms X 3
    Py_ssize_t num_walkers = PyObject_Length(coords);
    double* raw_data = _GetDoubleDataArray(coords);
    CHECKNULL(raw_data);

    int dims[1] = {num_walkers};
    PyObject* pots = _CreateArray(1, dims, "zeros", "float64");
    double* pot_data = _GetDoubleDataArray(pots);
    for ( int n = 0; n < num_walkers; n++) {
        pot_data[n] = _callPtrPot(potPtr, _getWalkerCoords(raw_data, n, num_atoms), atom_types);
    };

    return pots;

}

static PyMethodDef CPotentialLibMethods[] = {
    {"callPot", CPotentialLib_callPot, METH_VARARGS, "calls a C-based potential"},
    {"callPotVec", CPotentialLib_callPotVec, METH_VARARGS, "calls a C-based potential on multiple walkers"},
    {"callMPIPot", CPotentialLib_callPotMPI, METH_VARARGS, "calls a C-based potential on multiple walkers using MPI"},
    {"getMPIInfo", CPotentialLib_getMPIInfo, METH_VARARGS, "returns the MPI info for the processor"},
    {"closeMPI", CPotentialLib_closeMPI, METH_VARARGS, "closes out of MPI"},
    {NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION > 2

const char CPotentialLib_doc[] = "CPotentialLib is a little library for hooking C potentials into python";
static struct PyModuleDef CPotentialLibModule = {
    PyModuleDef_HEAD_INIT,
    "CPotentialLib",   /* name of module */
    CPotentialLib_doc, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
   CPotentialLibMethods
};

PyMODINIT_FUNC PyInit_CPotentialLib(void)
{
    return PyModule_Create(&CPotentialLibModule);
}
#else

PyMODINIT_FUNC initCPotentialLib(void)
{
    (void) Py_InitModule("CPotentialLib", CPotentialLibMethods);
}

#endif
