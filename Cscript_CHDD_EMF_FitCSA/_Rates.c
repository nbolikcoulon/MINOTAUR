#include <Python.h>
#include <numpy/arrayobject.h>
#include "Rates.h"


static char module_docstring[] =
    "This module calculates relaxation matrix and relaxation rates using C language.";
    
static char RelaxMat_docstring[] =
    "Calculates the relaxation matrix given a pset of arameters (field and dynamics as an example).";
static char R2calculation_docstring[] =
    "Calculates the carbon R2 given a set of parameters (field and dynamics as an example).";
static char R1calculation_docstring[] =
    "Calculates the carbon R1 given a set of parameters (field and dynamics as an example).";
static char Sigmacalculation_docstring[] =
    "Calculates the dipolar cross-relaxation given a set of parameters (field and dynamics as an example).";
    
    
static PyObject *RelaxMat_RelaxMat(PyObject *self, PyObject *args);
static PyObject *R2calculation_R2calculation(PyObject *self, PyObject *args);
static PyObject *R1calculation_R1calculation(PyObject *self, PyObject *args);
static PyObject *Sigmacalculation_Sigmacalculation(PyObject *self, PyObject *args);


static PyMethodDef module_methods[] = {
    {"RelaxMat", RelaxMat_RelaxMat, METH_VARARGS, RelaxMat_docstring},
    {"R2calculation", R2calculation_R2calculation, METH_VARARGS, R2calculation_docstring},
    {"R1calculation", R1calculation_R1calculation, METH_VARARGS, R1calculation_docstring},
    {"Sigmacalculation", Sigmacalculation_Sigmacalculation, METH_VARARGS, Sigmacalculation_docstring},
    {NULL, NULL, 0, NULL}
};







PyMODINIT_FUNC PyInit__RelaxMat(void)
{
    
    PyObject *module;
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_RelaxMat",
        module_docstring,
        -1,
        module_methods,
        NULL,
        NULL,
        NULL,
        NULL
    };
    module = PyModule_Create(&moduledef);
    if (!module) return NULL;

    /* Load `numpy` functionality. */
    import_array();

    return module;
}

PyMODINIT_FUNC PyInit__R2calculation(void)
{
    
    PyObject *module;
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_R2calculationt",
        module_docstring,
        -1,
        module_methods,
        NULL,
        NULL,
        NULL,
        NULL
    };
    module = PyModule_Create(&moduledef);
    if (!module) return NULL;

    /* Load `numpy` functionality. */
    import_array();

    return module;
}

PyMODINIT_FUNC PyInit__R1calculation(void)
{
    
    PyObject *module;
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_R1calculation",
        module_docstring,
        -1,
        module_methods,
        NULL,
        NULL,
        NULL,
        NULL
    };
    module = PyModule_Create(&moduledef);
    if (!module) return NULL;

    /* Load `numpy` functionality. */
    import_array();

    return module;
}

PyMODINIT_FUNC PyInit__Sigmacalculation(void)
{
    
    PyObject *module;
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_Sigmacalculation",
        module_docstring,
        -1,
        module_methods,
        NULL,
        NULL,
        NULL,
        NULL
    };
    module = PyModule_Create(&moduledef);
    if (!module) return NULL;

    /* Load `numpy` functionality. */
    import_array();

    return module;
}








static PyObject *RelaxMat_RelaxMat(PyObject *self, PyObject *args)
{
    double BVal, TaucVal;
    PyObject *X_obj, *OtherInputs_obj;
    
    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "dOdO", &BVal, &X_obj, &TaucVal, &OtherInputs_obj))
        return NULL;
        
    
    
    /* Interpret the input objects as numpy arrays. */
    PyObject *X_array = PyArray_FROM_OTF(X_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *OtherInputs_array = PyArray_FROM_OTF(OtherInputs_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    
    

    /* If that didn't work, throw an exception. */
    if (X_array == NULL || OtherInputs_array == NULL) {
        Py_XDECREF(X_array);
        Py_XDECREF(OtherInputs_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *Xarr    = (double*)PyArray_DATA(X_array);
    double *OtherInputsArr = (double*)PyArray_DATA(OtherInputs_array);
    
    /* Call the external C function to compute the matrix. */
    struct mat2d MatrixFull;                                   //Is system dependant!!!!!
    MatrixFull = RelaxMatrix( BVal, Xarr, TaucVal, OtherInputsArr );
    

    /* Clean up. */
    Py_DECREF(X_array);
    Py_DECREF(OtherInputs_array);
    
    int Lr = sizeof(MatrixFull.m) / sizeof(MatrixFull.m[0]);

    /* Build the output tuple */
    Py_ssize_t Len = Lr;
    PyObject *MatTransformed = PyTuple_New(Len);
    for (Py_ssize_t i = 0; i < Len; i++) {
        Py_ssize_t Len = Lr;
        PyObject *item = PyTuple_New(Len);
        for (Py_ssize_t j = 0; j < Len; j++) 
            PyTuple_SET_ITEM(item, j, PyFloat_FromDouble(MatrixFull.m[i][j]));
        PyTuple_SET_ITEM(MatTransformed, i, item);
        
    }
    
    
    PyObject *ret = Py_BuildValue("(O)", MatTransformed);
    Py_DECREF(MatTransformed);
    return ret;
}




static PyObject *R2calculation_R2calculation(PyObject *self, PyObject *args)
{
    double BVal, TaucVal;
    PyObject *X_obj, *OtherInputs_obj;
    
    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "dOdO", &BVal, &X_obj, &TaucVal, &OtherInputs_obj))
        return NULL;
        
    
    
    /* Interpret the input objects as numpy arrays. */
    PyObject *X_array = PyArray_FROM_OTF(X_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *OtherInputs_array = PyArray_FROM_OTF(OtherInputs_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    
    

    /* If that didn't work, throw an exception. */
    if (X_array == NULL || OtherInputs_array == NULL) {
        Py_XDECREF(X_array);
        Py_XDECREF(OtherInputs_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *Xarr    = (double*)PyArray_DATA(X_array);
    double *OtherInputsArr = (double*)PyArray_DATA(OtherInputs_array);
    
    /* Call the external C function to compute the rate. */
    double R2val;
    R2val = R2calculation( BVal, Xarr, TaucVal, OtherInputsArr );
    

    /* Clean up. */
    Py_DECREF(X_array);
    Py_DECREF(OtherInputs_array);
    
    PyObject *ret = Py_BuildValue("(d)", R2val);
    return ret;
}

static PyObject *R1calculation_R1calculation(PyObject *self, PyObject *args)
{
    double BVal,TaucVal;
    PyObject *X_obj, *OtherInputs_obj;
    
    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "dOdO", &BVal, &X_obj, &TaucVal, &OtherInputs_obj))
        return NULL;
        
    
    
    /* Interpret the input objects as numpy arrays. */
    PyObject *X_array = PyArray_FROM_OTF(X_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *OtherInputs_array = PyArray_FROM_OTF(OtherInputs_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    
    

    /* If that didn't work, throw an exception. */
    if (X_array == NULL || OtherInputs_array == NULL) {
        Py_XDECREF(X_array);
        Py_XDECREF(OtherInputs_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *Xarr    = (double*)PyArray_DATA(X_array);
    double *OtherInputsArr = (double*)PyArray_DATA(OtherInputs_array);
    
    /* Call the external C function to compute the rate. */
    double R1val;
    R1val = R1calculation( BVal, Xarr, TaucVal, OtherInputsArr );
    

    /* Clean up. */
    Py_DECREF(X_array);
    Py_DECREF(OtherInputs_array);
    
    PyObject *ret = Py_BuildValue("(d)", R1val);
    return ret;
}

static PyObject *Sigmacalculation_Sigmacalculation(PyObject *self, PyObject *args)
{
    double BVal, TaucVal;
    PyObject *X_obj, *OtherInputs_obj;
    
    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "dOdO", &BVal, &X_obj, &TaucVal, &OtherInputs_obj))
        return NULL;
        
    
    
    /* Interpret the input objects as numpy arrays. */
    PyObject *X_array = PyArray_FROM_OTF(X_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *OtherInputs_array = PyArray_FROM_OTF(OtherInputs_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    
    

    /* If that didn't work, throw an exception. */
    if (X_array == NULL || OtherInputs_array == NULL) {
        Py_XDECREF(X_array);
        Py_XDECREF(OtherInputs_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *Xarr    = (double*)PyArray_DATA(X_array);
    double *OtherInputsArr = (double*)PyArray_DATA(OtherInputs_array);
    
    /* Call the external C function to compute the rate. */
    double Sigmaval;
    Sigmaval = Sigmacalculation( BVal, Xarr, TaucVal, OtherInputsArr );
    

    /* Clean up. */
    Py_DECREF(X_array);
    Py_DECREF(OtherInputs_array);
    
    PyObject *ret = Py_BuildValue("(d)", Sigmaval);
    return ret;
}

