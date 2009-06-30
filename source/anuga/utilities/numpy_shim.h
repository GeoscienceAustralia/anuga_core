//******************************************************************************
// 'shim' code to use old Numeric API without deprecation warnings.
//
// This file is #included into every file that uses deprecated Numeric
// functions.  The entry points here are static, so no clash occurs.
//
// The deprecated Numeric functions are renamed from
//     PyArray_****()
// to
//     anuga_****()
// to pick up the definitions here.
//******************************************************************************

#include "numpy/noprefix.h"

#define	anuga_FromDimsAndData(nd, d, type, data)                             \
        anuga_FromDimsAndDataAndDescr(nd, d, PyArray_DescrFromType(type), data)


static PyObject *
anuga_FromDimsAndDataAndDescr(int nd, int *d, PyArray_Descr *descr, char *data)
{
    int i;
    npy_intp newd[MAX_DIMS];

    if (!PyArray_ISNBO(descr->byteorder))
        descr->byteorder = '=';

    for (i = 0; i < nd; i++)
    {
        newd[i] = (npy_intp) d[i];
    }

    return PyArray_NewFromDescr(&PyArray_Type, descr, nd, newd, NULL, data,
                                (data ? CARRAY : 0), NULL);
}

static PyObject *
anuga_FromDims(int nd, int *d, int type)
{
    PyObject *result;

    result = anuga_FromDimsAndDataAndDescr(nd, d,
                                           PyArray_DescrFromType(type), NULL);

    /*
     * Old FromDims set memory to zero --- some algorithms
     * relied on that.  Better keep it the same. If
     * Object type, then it's already been set to zero, though.
     */

    if (result && (PyArray_DESCR(result)->type_num != PyArray_OBJECT))
    {
        memset(PyArray_DATA(result), 0, PyArray_NBYTES(result));
    }

    return result;
}

