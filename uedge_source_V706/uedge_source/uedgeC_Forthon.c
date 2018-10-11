/* Created by Ron Cohen */
/* $Id: uedgeC_Forthon.c,v 7.0 2018/02/28 18:32:42 meyer8 Exp $ */
/* This is a stub module which calls the init functions of all of            */
/* the modules that are part of UEDGE. This is needed since the modules      */
/* depend on each other and so must be incorporated into one shared          */
/* object file.                                                              */
#include "Forthon.h"

/* ######################################################################### */
/* # Method list                                                             */
static struct PyMethodDef uedgeC_methods[] = {
  {NULL,NULL}};

/* ######################################################################### */
/* # The initialization function                                             */
void inituedgeC(void)
{
  PyObject *m, *d;
  PyObject *pystdout;
  m = Py_InitModule("uedgeC", uedgeC_methods);
  d = PyModule_GetDict(m);
  ErrorObject = PyString_FromString("uedgeC.error");
  PyDict_SetItemString(d, "error", ErrorObject);
  if (PyErr_Occurred())
    Py_FatalError("can not initialize module uedgeC");

  pystdout = PySys_GetObject("stdout");
  PyFile_WriteString("Forthon edition\n",pystdout);

  import_array();

  initcompy();
  initaphpy();
  initapipy();
  initbbbpy();
  initflxpy();
  initgrdpy();
  initsvrpy();
  initwdfpy();
}


