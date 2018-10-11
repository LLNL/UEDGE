#include "pybasis.h"

static struct PyMethodDef uedge_methods[] = {
  {NULL,NULL}};

extern void initpybasisC();
extern PyBasisObject * initcom();
extern PyBasisObject * initbbb();
extern PyBasisObject * initaph();
extern PyBasisObject * initapi();
extern PyBasisObject * initflx();
extern PyBasisObject * initgrd();

void inituedge(void) {
  PyObject *m, *d;

  /* Initialize the PyBasis submodule */
  initpybasisC();

  /* Add uedge to Python's list of built-in modules. */
  m = Py_InitModule("uedge",uedge_methods);

  /* Set the uedge packages as module attributes */
  d = PyModule_GetDict(m);
  PyDict_SetItemString(d,"bbb",(PyObject *)initbbb());
  PyDict_SetItemString(d,"aph",(PyObject *)initaph());
  PyDict_SetItemString(d,"com",(PyObject *)initcom());
  PyDict_SetItemString(d,"api",(PyObject *)initapi());
  PyDict_SetItemString(d,"flx",(PyObject *)initflx());
  PyDict_SetItemString(d,"grd",(PyObject *)initgrd());
  if (PyErr_Occurred())
    Py_FatalError("can not initialize module uedge");
}
