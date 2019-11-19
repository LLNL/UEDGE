/* Created by Ron Cohen */
/* $Id: uedgeC_Forthon.c,v 7.1 2019/11/01 22:20:49 meyer8 Exp $ */
/* This is a stub module which calls the init functions of all of            */
/* the modules that are part of UEDGE. This is needed since the modules      */
/* depend on each other and so must be incorporated into one shared          */
/* object file.                                                              */
#include <Python.h>
#define NPY_NO_DEPRECATED_API 8
#include <numpy/arrayobject.h>


#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif

static PyObject *ErrorObject;

#if PY_MAJOR_VERSION < 3
  extern PyMODINIT_FUNC initcompy(void);
  extern PyMODINIT_FUNC initaphpy(void);
  extern PyMODINIT_FUNC initapipy(void);
  extern PyMODINIT_FUNC initbbbpy(void);
  extern PyMODINIT_FUNC initflxpy(void);
  extern PyMODINIT_FUNC initgrdpy(void);
  extern PyMODINIT_FUNC initsvrpy(void);
  extern PyMODINIT_FUNC initwdfpy(void);
#endif

/* ######################################################################### */
/* # Method list                                                             */
static struct PyMethodDef uedgeC_methods[] = {
  {NULL,NULL}};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "uedgeC", /* m_name */
  "uedgeC", /* m_doc */
  -1,                  /* m_size */
  uedgeC_methods,    /* m_methods */
  NULL,                /* m_reload */
  NULL,                /* m_traverse */
  NULL,                /* m_clear */
  NULL,                /* m_free */
  };
#endif



/* ######################################################################### */
/* # The initialization function                                             */
#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_uedgeC(void)
#else
PyMODINIT_FUNC inituedgeC(void)
#endif
{

  PyObject *m, *d;
  /*PyObject *pystdout;*/
#if PY_MAJOR_VERSION >= 3
  m = PyModule_Create(&moduledef);
#else
  m = Py_InitModule("uedgeC", uedgeC_methods);
#endif

  d = PyModule_GetDict(m);
  ErrorObject = PyErr_NewException("uedgeC.error",NULL,NULL);
  PyDict_SetItemString(d, "error", ErrorObject);
  if (PyErr_Occurred())
    Py_FatalError("can not initialize module uedgeC");

  /*pystdout = PySys_GetObject("stdout");*/
  /*PyFile_WriteString("Forthon edition\n",pystdout);*/

  import_array();

#if PY_MAJOR_VERSION < 3
  initcompy();
  initaphpy();
  initapipy();
  initbbbpy();
  initflxpy();
  initgrdpy();
  initsvrpy();
  initwdfpy();
#endif
#if PY_MAJOR_VERSION >= 3
  return m;
#else
  return;
#endif

}


