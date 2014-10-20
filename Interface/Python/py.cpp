/*****************************************************
 *    Matlab to Python Interface
 *****************************************************/

#include <mex.h>
#include <python/Python.h>
#include <string.h>
#include <algorithm>
#include <map>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <iostream>
//#define DEBUG(msg) std::cout << msg << std::endl; std::cout.flush();
#define DEBUG(msg)

//#pragma GCC diagnostic ignored "-fpermissive"


static PyObject* MPy_main = NULL;
static PyObject* MPy_globals = NULL;

enum MPy_status_type {Uninitialized = 0, Initialized};
static MPy_status_type MPy_status = Uninitialized;


/*********************************************************
 *    Array Conversions
 *********************************************************/

/*************************
 * Matlab -> Python
 * **********************/

static PyObject* MPy_mat_to_python(const mxArray* mx);



static bool MPy_is_string(const mxArray *mx)
{
    return mxIsChar(mx) && mxGetM(mx) == 1;
}


static bool MPy_is_scalar(const mxArray *mx)
{
    int nd;
    const size_t *dims;
    nd = mxGetNumberOfDimensions(mx);
    dims = mxGetDimensions(mx);

    return (nd == 1 && dims[0] == 1)
        || (nd == 2 && dims[0] == 1 && dims[1] == 1);
}


static void MPy_reshape_dims(size_t ndims, size_t* dims) {
   //account for c vs fortran order
   std::reverse(dims, dims + ndims);
}



static PyObject* MPy_char_to_python(const mxArray* mx)
{   
   DEBUG("MPy_char_to_python")
   
   char* str;
   PyObject* result;
 
   str =mxArrayToString(mx);
   result = PyString_FromString(str);

   mxFree(str);   
   return result;
}


static PyObject* MPy_scalar_to_python(const mxArray* mx)
{   
   DEBUG("MPy_scalar_to_python")

   char *data     = (char*) mxGetData(mx);
   char *imagData = (char*) mxGetImagData(mx);
   mxClassID cls = mxGetClassID(mx);
   
   PyObject* result = NULL;
   
   if (imagData == NULL) {
   
      switch (cls) {
         case mxLOGICAL_CLASS:
            result = PyBool_FromLong(*((bool*) data)); break;
         case mxDOUBLE_CLASS:
            result = PyFloat_FromDouble(*((double*) data)); break;
         case mxSINGLE_CLASS:
            result = PyFloat_FromDouble(*((float*) data)); break; 
         case mxINT8_CLASS:
            result = PyInt_FromLong(*((char*) data)); break;  
         case mxUINT8_CLASS:
            result = PyInt_FromLong(*((unsigned char*) data)); break;           
         case mxINT16_CLASS:
            result = PyInt_FromLong(*((short*) data)); break;  
         case mxUINT16_CLASS:
            result = PyInt_FromLong(*((unsigned short*) data)); break;
         case mxINT32_CLASS:
            result = PyInt_FromLong(*((int*) data)); break;  
         case mxUINT32_CLASS:
            result = PyLong_FromLongLong(*((unsigned int*) data)); break;        
         case mxINT64_CLASS:
            result = PyLong_FromLongLong(*((long long*) data)); break;  
         case mxUINT64_CLASS:
            result = PyLong_FromLongLong(*((unsigned long long*) data)); break;        
         default:
            mexPrintf("Mpy: MPy_value_to_python: mxClassID %d not a numeric type!", cls);
      }
   
   } else {
      switch (cls) {
         case mxDOUBLE_CLASS:
            result = PyComplex_FromDoubles(*((double*) data), *((double*) imagData)); break;
         case mxSINGLE_CLASS:
            result = PyComplex_FromDoubles(*((float*) data), *((float*) imagData)); break; 
         case mxINT8_CLASS:
            result = PyComplex_FromDoubles(*((char*) data), *((char*) imagData)); break;  
         case mxUINT8_CLASS:
            result = PyComplex_FromDoubles(*((unsigned char*) data), *((unsigned char*) imagData)); break;           
         case mxINT16_CLASS:
            result = PyComplex_FromDoubles(*((short*) data), *((short*) imagData)); break;  
         case mxUINT16_CLASS:
            result = PyComplex_FromDoubles(*((unsigned short*) data), *((unsigned short*) imagData)); break;
         case mxINT32_CLASS:
            result = PyComplex_FromDoubles(*((int*) data), *((int*) imagData)); break;  
         case mxUINT32_CLASS:
            result = PyComplex_FromDoubles(*((unsigned int*) data), *((unsigned int*) imagData)); break;        
         case mxINT64_CLASS:
            result = PyComplex_FromDoubles(*((long long*) data), *((long long*) imagData)); break;  
         case mxUINT64_CLASS:
            result = PyComplex_FromDoubles(*((unsigned long long*) data), *((unsigned long long*) imagData)); break;        
         default:
            mexPrintf("Mpy: MPy_value_to_python: mxClassID %d not a numeric complex type!", cls);
      }
   }
   
   return result;
}



static PyObject* MPy_array_to_python(const mxArray* mx)
{
   DEBUG("MPy_array_to_python")
   
   const void *real, *imag;
   size_t n, nd, *dims;
   
   PyObject* result;
   npy_intp npytype = 0;

   real = mxGetData(mx);
   imag = mxGetImagData(mx);
   n    = mxGetNumberOfElements(mx);
   nd   = mxGetNumberOfDimensions(mx);
   dims = (size_t*) mxGetDimensions(mx);

   //check if vector 
   if (nd == 2 && dims[1] == 1)
      nd = 1;

   MPy_reshape_dims(nd, dims);

   switch (mxGetClassID(mx)) {
      case mxDOUBLE_CLASS:
            npytype = imag ? NPY_CDOUBLE : NPY_DOUBLE;
            if (imag) {
                  result = PyArray_SimpleNew(nd, (npy_intp*) dims, npytype);
                  result = (PyObject*) PyArray_GETCONTIGUOUS((PyArrayObject*) result);
                  Py_DECREF(result);
                  for (size_t i = 0; i < n; i++) {
                     ((double *)PyArray_DATA((PyArrayObject *)result))[2 * i] = ((double *)real)[i];
                     ((double *)PyArray_DATA((PyArrayObject *)result))[2 * i + 1] = ((double *)imag)[i];
                  }
                  return result;
            }
            break;
            
      case mxSINGLE_CLASS:
            npytype = imag ? NPY_CFLOAT : NPY_FLOAT;
            if (imag) {
                  result = PyArray_SimpleNew(nd, (npy_intp*) dims, npytype);
                  result = (PyObject*) PyArray_GETCONTIGUOUS((PyArrayObject*) result);
                  Py_DECREF(result);
                  for (size_t i = 0; i < n; i++) {
                     ((float *)PyArray_DATA((PyArrayObject *)result))[2 * i] = ((float *)real)[i];
                     ((float *)PyArray_DATA((PyArrayObject *)result))[2 * i + 1] = ((float *)imag)[i];
                  }
                  return result;
            }
            break;
            
      case mxINT8_CLASS:   npytype = NPY_BYTE;     break;
      case mxUINT8_CLASS:  npytype = NPY_UBYTE;    break;
      case mxINT16_CLASS:  npytype = NPY_SHORT;    break;
      case mxUINT16_CLASS: npytype = NPY_USHORT;   break;
      case mxINT32_CLASS:  npytype = NPY_INT;      break;
      case mxUINT32_CLASS: npytype = NPY_UINT;     break;
      case mxINT64_CLASS:  npytype = NPY_LONGLONG; break;
      case mxUINT64_CLASS: npytype = NPY_ULONGLONG;break;
      case mxLOGICAL_CLASS:npytype = NPY_BOOL;     break;
      case mxCHAR_CLASS:   npytype = NPY_CHAR;     break;
      default:
         mexPrintf("MPy: MPy_numeric_to_python, Unexpected class ID %d !", mxGetClassID(mx));
         return NULL;
   }
   
   DEBUG("MPy_numeric_to_python npytype=" << npytype << " mxclass=" << mxGetClassID(mx))

   result = PyArray_SimpleNew(nd, (npy_intp*) dims, npytype);
   if (result == NULL) {
      mexPrintf("MPy: MPy_numeric_to_python, cannot create array!");
      return NULL;
   }
   
   result = (PyObject*) PyArray_GETCONTIGUOUS((PyArrayObject*) result);
   Py_DECREF(result);
   //DEBUG(result->ob_refcnt)

   memcpy(PyArray_DATA((PyArrayObject*) result), real, n * PyArray_ITEMSIZE((PyArrayObject*)result));
  
   return result;
}


static PyObject* MPy_value_to_python(const mxArray* mx)
{   
   DEBUG("MPy_value_to_python")
   
   if (MPy_is_string(mx))
      return MPy_char_to_python(mx);
   else if (MPy_is_scalar(mx))
      return MPy_scalar_to_python(mx);
   else
      return MPy_array_to_python(mx); 
}




static PyArray_Descr* MPy_struct_dtype(const mxArray* mx)
{
   DEBUG("MPy_struct_dtype")
   
   int nfields;
   
   PyArray_Descr *descr = NULL;
   PyObject *list;

   nfields = mxGetNumberOfFields(mx);
   list = PyList_New(nfields);
   for (int i = 0; i < nfields; i++) {
      DEBUG("MPy_struct_dtype: fieldname:" << mxGetFieldNameByNumber(mx, i));  
      PyList_SetItem(list, i, Py_BuildValue("(s,s)", mxGetFieldNameByNumber(mx, i), "object"));
   }
   
   PyDict_SetItemString(MPy_globals, "test", list);
   PyRun_String("print test" , Py_file_input, MPy_globals, MPy_globals);
   
   
   if (PyArray_DescrConverter(list, &descr) != NPY_SUCCEED) {
      mexPrintf("MPy: MPy_struct_dtype: Failed to determine dtype of struct array.");
   }
   
   Py_DECREF(list);
   

   return descr;
}

static PyObject* MPy_struct_to_python(const mxArray* mx)
{
   DEBUG("MPy_struct_to_python")

   size_t n, nd, *dims;
   int nfields;

   PyArray_Descr *dtype;
   PyArrayIterObject *iter;
   PyObject *result;
   PyObject ** data;
      
   n       = mxGetNumberOfElements(mx);
   nfields = mxGetNumberOfFields(mx);
   nd      = mxGetNumberOfDimensions(mx);
   dims    = (size_t*) mxGetDimensions(mx);
   
   if (nd == 2 && dims[1] == 1)
      nd = 1;
   
   
   
 
    PyObject* op = Py_BuildValue("[(s, s), (s, s)]", "aaaa", "i4", "bbbb", "f4");
    //  PyObject* op = Py_BuildValue("{s:(si),s:(si)}", "aaaa", "i4",0, "bbbb", "f4",4 );
      PyDict_SetItemString(MPy_globals, "testa", op);
      PyRun_String("print 'doll'; print testa" , Py_file_input, MPy_globals, MPy_globals);
   

       int x = PyArray_DescrConverter(op, &dtype);
       DEBUG(x << "," << NPY_SUCCEED);
      //PyArray_DescrAlignConverter(op, &dtype );
      //Py_DECREF(op);
 
     PyDict_SetItemString(MPy_globals, "testb", dtype->names);
      PyRun_String("print 'doll'; print testb" , Py_file_input, MPy_globals, MPy_globals);
       
       DEBUG("JOJO")
      
     // PyObject_Print( (PyObject*)dtype , stdout, 0);
       
       size_t dims2[1] = {3};
       
      result = PyArray_SimpleNewFromDescr(1, (npy_intp*) dims2, dtype);
      Py_DECREF(dtype); 
      DEBUG("result ptr = " << result);
   
      PyDict_SetItemString(MPy_globals, "testc", result);
      PyRun_String("print testc" , Py_file_input, MPy_globals, MPy_globals);
      DEBUG("TT")
   
       
   

   
   
   MPy_reshape_dims(nd, dims);

   dtype = MPy_struct_dtype(mx);
   DEBUG("element size: " << dtype->elsize)
   
   result = PyArray_SimpleNewFromDescr(nd, (npy_intp*) dims, dtype);
   Py_DECREF(dtype); 
   DEBUG("result ptr = " << result);
   
   PyDict_SetItemString(MPy_globals, "test2", result);
   PyRun_String("print test2" , Py_file_input, MPy_globals, MPy_globals);
   DEBUG("TT")
   
   
   iter = (PyArrayIterObject *)PyArray_IterNew(result);
   if (iter == NULL) {
      mexPrintf("MPy: MPy_struct_to_python: Failed to convert struct to python\n");
      Py_CLEAR(result);
      return result;
   }

   PyDict_SetItemString(MPy_globals, "test", result);
   PyRun_String("print test" , Py_file_input, MPy_globals, MPy_globals);
   DEBUG("TEST")
   
   for (size_t i = 0; i < n; i++) {
      for (int k = 0; k < nfields; k++) {
         
         data = (PyObject**) PyArray_ITER_DATA(iter);
         //data = (PyArrayObject*) Py_None; Py_INCREF(Py_None);
         *data = MPy_mat_to_python(mxGetFieldByNumber(mx, i, k));
         PyArray_ITER_NEXT(iter);
      }
   }
 
   PyDict_SetItemString(MPy_globals, "test", result);
   PyRun_String("print 'hellow'; print test" , Py_file_input, MPy_globals, MPy_globals);
   DEBUG("TEST")
   
 
   return result;
}



/*
static PyObject* MPy_struct_to_python_dict(const mxArray* mx)
{
    int nfields;
    int field_number;
    
    PyObject *result;

    nfields = mxGetNumberOfFields(mx);

    result = PyDict_New();

    for (field_number = 0; field_number < nfields; ++field_number) {
        const char *name;  
        name = mxGetFieldNameByNumber(mx, field_number);

        if (MPy_is_scalar(mx)) {
            mxArray *a = mxGetFieldByNumber(mx, 0, field_number);
            PyDict_SetItemString(result, name, MPy_mat_to_python(a));
        } else {
            //todo
        }
    }
    return result;
}
*/


static PyObject* MPy_cell_to_python(const mxArray* mx)
{
   DEBUG("MPy_cell_to_python")
   
   const mxArray *cell;
   size_t n, nd, *dims;
   
   PyObject *result, **itemptr;

   n    = mxGetNumberOfElements(mx);
   nd   = mxGetNumberOfDimensions(mx);
   dims = (size_t*) mxGetDimensions(mx);
   
   if (nd == 2 && dims[1] == 1)
      nd = 1;
   
   MPy_reshape_dims(nd, dims); 
   
   
   result = PyArray_SimpleNew(nd, (npy_intp*) dims, NPY_OBJECT);
   result = (PyObject*) PyArray_GETCONTIGUOUS((PyArrayObject*) result);
   Py_DECREF(result);
     
   itemptr = (PyObject**) (PyArrayObject**) PyArray_DATA((PyArrayObject*) result);
   
   for (size_t i = 0; i < n; i++) {
      cell = mxGetCell(mx, i);
      if (cell) {
         *itemptr = MPy_mat_to_python(cell);
      } else {
         *itemptr = Py_None;
         Py_INCREF(*itemptr);
      }
      itemptr++;
   }

   return result;
}

static PyObject* MPy_function_to_python(const mxArray* mx)
{
   DEBUG("MPy_function_to_python")
   
   return PyCObject_FromVoidPtr((void *)mx, NULL);
}



static PyObject* MPy_mat_to_python(const mxArray* mx)
{
   DEBUG("MPy_mat_to_python")

   switch (mxGetClassID(mx)) {
      case mxCHAR_CLASS: 
      case mxDOUBLE_CLASS:
      case mxSINGLE_CLASS:
      case mxINT8_CLASS:
      case mxUINT8_CLASS:
      case mxINT16_CLASS:
      case mxUINT16_CLASS:
      case mxINT32_CLASS:
      case mxUINT32_CLASS:
      case mxINT64_CLASS:
      case mxUINT64_CLASS:
      case mxLOGICAL_CLASS:
         return MPy_value_to_python(mx);
      case mxSTRUCT_CLASS:
         return MPy_struct_to_python(mx);
      case mxCELL_CLASS:
         return MPy_cell_to_python(mx);
      case mxFUNCTION_CLASS:
         return MPy_function_to_python(mx);
            
      case mxUNKNOWN_CLASS:
      default:
         mexPrintf("MPy: MPy_to_python: unknown ClassID %d to convert to Python.", mxGetClassID(mx));
         return NULL;
   }
}







/*************************
 * Python -> Matlab
 * **********************/

static mxArray *MPy_python_to_mat(PyObject *object);

static char * MPy_dimstrcpy(char *buffer, size_t ndims, size_t* dims)
{
   DEBUG("MPy_dimstrcpy")
   
   int idx;
   char *ptr = buffer;
   *ptr++='(';
   for (idx = 0; idx<ndims; idx++)
         ptr += sprintf(ptr,"%d,",dims[idx]);
   ptr[-1] = ')';
   return buffer;
}

static long MPy_get_dtype(PyObject *dtype)
{
   DEBUG("MPy_get_dtype")

   PyObject *dtype_num = PyObject_GetAttrString(dtype, "num");
   
   long npy_type = NPY_NOTYPE;
   if (! (dtype_num && PyInt_Check(dtype_num))) {
         Py_DECREF(dtype_num);
         mexErrMsgTxt("MPy: numpy.ndarray.dtype did not have a proper num!\n");
   }
   
   npy_type = PyInt_AsLong(dtype_num);
   Py_DECREF(dtype_num);
   
   return npy_type;
}



static mxArray *MPy_none_to_mat(PyObject* obj)
{
    size_t dims[1] = { 0 };
    return mxCreateCellArray(1, dims);
}


static mxArray *MPy_bool_to_mat(PyObject* obj)
{
    size_t dims[1] = { 1 };
    mxArray *result;
    
    result = mxCreateNumericArray(1, dims, mxLOGICAL_CLASS, mxREAL);

    if (PyObject_Compare(obj, Py_True) == 0) {
        *((char*)mxGetData(result)) = 1;
    } else {
        *((char*)mxGetData(result)) = 0;
    }
    return result;
}

static mxArray* MPy_string_to_mat(PyObject* obj)
{
    DEBUG("MPy_string_to_mat")

    size_t dims[2];
    
    char* buffer;
    Py_ssize_t length;

    mxArray *result;
    mxChar* ptr;
    
    PyString_AsStringAndSize(obj, &buffer, &length);

    dims[0] = 1;
    dims[1] = length;
    result = mxCreateCharArray(2, dims);
    ptr    = (mxChar*) mxGetData(result);
 
    for (; length > 0; --length) {
        *ptr = *buffer;
        ++ptr; ++buffer;
    }
    //memcpy(mxGetData(result), buffer, ((int) length) * sizeof(char)); 
    
    return result;
}

        
static mxArray* MPy_complex_to_mat(PyObject* obj)
{
    Py_complex complex;
    mxArray *result = NULL;
    size_t dims[1] = { 1 };
    complex = PyComplex_AsCComplex(obj);
    result = mxCreateNumericArray(1, dims, mxDOUBLE_CLASS, mxCOMPLEX);
    *mxGetPr(result) = complex.real;
    *mxGetPi(result) = complex.imag;
    
    return result;
}


static mxArray* MPy_float_to_mat(PyObject* obj)
{
   mxArray* result;
   size_t dims[1] = { 1 };
   result = mxCreateNumericArray(1, dims, mxDOUBLE_CLASS, mxREAL);
   *((double*)mxGetData(result)) = PyFloat_AS_DOUBLE(obj);

   return result;
}

static mxArray* MPy_int_to_mat(PyObject* obj)
{
    mxArray* result;
    size_t dims[1] = { 1 };
    result = mxCreateNumericArray(1, dims, mxINT32_CLASS, mxREAL);
    *((long*)mxGetData(result)) = PyInt_AS_LONG(obj);
    
    return result;
}



static mxArray* MPy_ndstruct_to_mat(PyObject *object, PyObject *dtype, size_t ndims, const size_t *dims, size_t length)
{
   DEBUG("MPy_ndstruct_to_mat");
     
   PyObject *contiguous = NULL;
   PyObject *fields = NULL;
   int nfields, npy_type;
   char **keys = NULL;
   PyObject* stringized = NULL;
   PyObject* names = NULL;
   PyObject** data;
   
   mxArray *result = NULL;
   
   fields = PyObject_GetAttrString(dtype,"fields");
   if (! PyMapping_Check(fields)) {
      mexPrintf("MPy: Failed to get fields from dtype\n");
      goto exit;
   }
   names = PyObject_GetAttrString(dtype,"names");
   if (! PyTuple_Check(names)) {  
      mexPrintf("MPy: Failed to get names from dtype\n");
      goto exit;
   }
   nfields = PyMapping_Size(fields);

   keys = (char**) mxCalloc(nfields, sizeof(char*));
   
   for (int i=0; i<nfields; i++) {
      PyObject *key = PyTuple_GetItem(names, i);
      if (!key) {
         mexPrintf("MPy: Failed to get key from tuple for field # %d\n",i);
         goto exit;
      }
      stringized = PyObject_Str(key);
      if (! PyString_Check(stringized)) {
         mexPrintf("MPy: Key # %d couldn't be stringized\n",i);
         goto exit;
      }
      keys[i] = PyString_AsString(stringized);
      Py_DECREF(stringized);
   }


   for (int i=0; i<nfields; i++) {
      PyObject *subtype;
      PyObject *item;

      item = PyMapping_GetItemString(fields, keys[i]);
      if (! PyTuple_Check(item)) {
         mexPrintf("MPy: Did not get a tuple item for field %s\n", (char*) keys[i]);
         goto exit;
      }
      
      subtype = PyTuple_GetItem(item,0);
      if (! subtype) {
         mexPrintf("MPy: Failed to get sub-type from tuple for field %s\n",keys[i]);
         Py_DECREF(item);
         goto exit;
      }
      if (! PyArray_DescrCheck(subtype)) {
         mexPrintf("MPy: Subtype for field %s was not an array descriptor.\nIt was type %s\n",keys[i],subtype->ob_type->tp_name);
         Py_DECREF(item);
         goto exit;
      }
      
      npy_type = MPy_get_dtype(subtype);
      Py_DECREF(item);
      
      if (npy_type != NPY_OBJECT) {
         mexPrintf("MPy: Subtype for field %s is not object - not supported\n",keys[i]);
         goto exit;
      }
   }
   
   result = mxCreateStructArray(ndims, dims, nfields, (const char**) keys);
   if (! result) {
         mexPrintf("MPy: Failed to allocate mx struct array");
         goto exit;
   }
   
   contiguous = (PyObject*) PyArray_GETCONTIGUOUS((PyArrayObject*) object);
   data = (PyObject **) PyArray_DATA((PyArrayObject*) contiguous);

   for (size_t i = 0; i < length; i++) {
      for (int k = 0; k < nfields; k++) {
         mxArray *sub;
         sub = MPy_python_to_mat(*data++);
         if (! sub) {
            mexPrintf("MPy: Failed to encode object at position %d, field %s",i,keys[k]);
            goto exit;
         }
         mxSetFieldByNumber(result, i, k, sub);
      }
   }
   
exit:
   Py_XDECREF(fields);
   Py_XDECREF(names);
   Py_XDECREF(contiguous);
   Py_XDECREF(stringized);
   if (keys) mxFree(keys);

   return result;
}



static mxArray* MPy_ndarray_to_mat(PyObject *object)
{
   DEBUG("MPy_ndarray_to_mat");

   Py_ssize_t ndims;
   PyObject *shape = NULL;
   PyObject *dtype = NULL;
   PyObject *iter = NULL;
   PyObject *contiguous = NULL;
   PyObject *dtype_num = NULL;
   mxClassID mxClass = (mxClassID) 0;

   mxArray *result = NULL;
   void *data = NULL;
   int length=1;
   size_t* dims;
   
   long npy_type = NPY_NOTYPE;
   char buffer[1024];

   /* Get the array dimensions */

   shape = PyObject_GetAttrString(object, "shape");
   if (! shape) {
      mexPrintf("MPy: numpy.ndarray did not have a shape\n");
      goto exit;
   }
   if (! PyTuple_Check(shape)) {
      mexPrintf("MPy: numpy.ndarray did not have a tuple shape\n");
      goto exit;
   }
   
   ndims = PyTuple_Size(shape);
   dims = (size_t *) mxMalloc(sizeof(size_t)*ndims);
   for(Py_ssize_t i=0; i< ndims; i++) {
      PyObject *obDim = PyTuple_GetItem(shape,i);
      if (! PyInt_Check(obDim)) {
         mexPrintf("MPy: Array shape is not a tuple of integers\n");
         goto exit;
      }
      dims[i] = (size_t) PyInt_AsLong(obDim);
      length *= dims[i];
   }
   
   // correct matlab vs python arrays
   MPy_reshape_dims(ndims, dims);

   
   // Handle dtype
   DEBUG("MPy_ndarray_to_mat: dtype")
   
   dtype = PyObject_GetAttrString(object, "dtype");
   if (! dtype) {
      mexPrintf("MPy: numpy.ndarray did not have a dtype\n");
      goto exit;
   }
   npy_type = MPy_get_dtype(dtype);
   
   
   DEBUG("MPy_ndarray_to_mat: switch")
   
   switch(npy_type) {
      case NPY_VOID: 
         if (PyArray_HASFIELDS((PyArrayObject*) object)) {
            result = MPy_ndstruct_to_mat(object, dtype, ndims, dims, (size_t) length);
         } else {
            mexPrintf("MPy: NPY_VOID array does not have fields\n");
         }
         break;
            
      case NPY_OBJECT: {
         PyObject *iter = PyObject_GetAttrString(object,"flat");

         if (! (iter && PyIter_Check(iter))) {
            mexPrintf("MPy: numpy.ndarray.flat did not return a proper iterator\n");
            goto exit;
         }
         
         result = mxCreateCellArray(ndims, dims);
         if (! result) {
            strcpy(buffer,"MPy: mxCreateCellArray returned NULL for an array of dimension ");
            strcpy(MPy_dimstrcpy(buffer+strlen(buffer), ndims, dims),"\n");
            mexPrintf(buffer);
            goto exit;
         }
         
         for (int idx = 0; idx < length; idx++) {
            PyObject *py_value;
            mxArray *mx_value;
            py_value = (*iter->ob_type->tp_iternext)(iter);
            if (! py_value) {
               mexPrintf("MPy: numpy.ndarray's __iter__'s next function failed\n");
               goto exit;
            }
            
            Py_INCREF(py_value);
            mx_value = MPy_python_to_mat(py_value);
            Py_DECREF(py_value);
            if (! mx_value) {
               mexPrintf("MPy: ndarray_to_mxarray failed while trying to translate a subobject to an mxarray\n");
               goto exit;
            }
            
            mxSetCell(result, idx, mx_value);
         }
         break;
      }
      
      case NPY_DOUBLE:
         mxClass = mxDOUBLE_CLASS;
         goto NPY_ALL;
            
      case NPY_FLOAT:
         mxClass = mxSINGLE_CLASS;
         goto NPY_ALL;
            
      case NPY_BYTE:
         mxClass = mxINT8_CLASS;
         goto NPY_ALL;
            
      case NPY_UBYTE:
         mxClass = mxUINT8_CLASS;
         goto NPY_ALL;
            
      case NPY_SHORT:
         mxClass = mxINT16_CLASS;
         goto NPY_ALL;
            
      case NPY_USHORT:
         mxClass = mxUINT16_CLASS;
         goto NPY_ALL;
            
      case NPY_INT:
      case NPY_LONG:
         mxClass = mxINT64_CLASS;
         goto NPY_ALL;
            
      case NPY_UINT:
      case NPY_ULONG:
         mxClass = mxUINT64_CLASS;
         goto NPY_ALL;
            
      case NPY_LONGLONG:
         mxClass = mxINT64_CLASS;
         goto NPY_ALL;
            
      case NPY_ULONGLONG:
         mxClass = mxUINT64_CLASS;
       
      NPY_ALL:
         DEBUG("MPy_ndarray_to_mat: ALL")

         contiguous = (PyObject*) PyArray_GETCONTIGUOUS((PyArrayObject*)object); 
         result = mxCreateNumericArray(ndims, (size_t*) dims, mxClass, mxREAL);
         if (! result) {               
            mexPrintf("MPy: cannot create maltab array.\n");
            goto exit;
         }

         memcpy(mxGetData(result), PyArray_DATA((PyArrayObject*) contiguous), length * PyArray_ITEMSIZE((PyArrayObject*) contiguous));
         
         Py_DECREF(contiguous);
         break;
            
            
      case NPY_CFLOAT: {
            float *real,*imaginary;
            float *npy_data;

            contiguous = (PyObject*) PyArray_GETCONTIGUOUS((PyArrayObject*)object);
            result = mxCreateNumericArray(ndims, (size_t*) dims, mxSINGLE_CLASS, mxCOMPLEX);
            
            real = (float *)mxGetData(result);
            imaginary = (float *)mxGetImagData(result);
            
            npy_data = (float *)PyArray_DATA((PyArrayObject*) contiguous);

            for (int i = 0; i < length; i++) {
               real[i] = npy_data[2*i];
               imaginary[i] = npy_data[2*i+1];
            }
            
            Py_DECREF(contiguous);
            break;
      }

      case NPY_CDOUBLE: {
         double *real,*imaginary;
         double *npy_data;

         contiguous = (PyObject*) PyArray_GETCONTIGUOUS((PyArrayObject*)object);

         result = mxCreateNumericArray(ndims, (size_t*) dims, mxDOUBLE_CLASS, mxCOMPLEX);

         real = (double *)mxGetData(result);
         imaginary = (double *)mxGetImagData(result);

         npy_data = (double *)PyArray_DATA((PyArrayObject*) contiguous);

         for (int i = 0; i < length; i++) {
            real[i] = npy_data[2*i];
            imaginary[i] = npy_data[2*i+1];
         }

         Py_DECREF(contiguous);
         break;
      }
      
      default: {
         PyObject *dtypename=PyObject_GetAttrString(dtype,"str");
         mexPrintf("MPy: Unhandled dtype: %s",PyString_AsString(dtypename));
         Py_DECREF(dtypename);
      }
   }
   
   DEBUG("MPy_ndarray_to_mat: exit")
   
        
exit:
   if (dims) mxFree(dims);
   Py_XDECREF(shape);
   Py_XDECREF(dtype);
   Py_XDECREF(dtype_num);
   Py_XDECREF(iter);

   DEBUG("MPy_ndarray_to_mat: return")
   
   return result;
}


  
mxArray *MPy_dict_to_mat(PyObject* obj)
{

    size_t dims[1] = { 1 };
    PyObject *items;
    int nitems, k;
    char *buf;
    Py_ssize_t len;
    char **fieldnames;
    PyObject *repr = NULL;
    
    mxArray* result = NULL;

    items = PyDict_Items(obj);
    if (!items) goto error;

    nitems = PyList_Size(items);
    fieldnames = (char**) mxCalloc(nitems, sizeof(char*));

    for (k = 0; k < nitems; ++k) {
        PyObject *o;
        
        o = PyList_GetItem(items, k);
        if (o == NULL) {
           mexPrintf("MPy: Failed to get python item from dict.\n");
           goto error;
        }

        o = PyTuple_GetItem(o, 0);
        if (o == NULL) {
           mexPrintf("MPy: Failed to get python item name from dict.\n");
           goto error;
        }

        if (PyString_Check(o)) {
            PyString_AsStringAndSize(o, &buf, &len);
        } else {
            repr = PyObject_Repr(o);
            buf = PyString_AsString(repr);
            len = strlen(buf);
        }

        fieldnames[k] = (char*) mxCalloc(len + 1, sizeof(char));
        memcpy(fieldnames[k], buf, len);
        fieldnames[k][len] = '\0';

        Py_CLEAR(repr);
    }

    result = mxCreateStructArray(1, dims, nitems, (const char**)fieldnames);
    if (!result) {        
       mexPrintf("MPy: Failed to create matalb array from dict.\n");
       goto error;
    }

    for (k = 0; k < nitems; ++k) {
        PyObject *o;
        o = PyList_GetItem(items, k);
        if (o == NULL) {
           mexPrintf("MPy: Failed to get python item from dict.\n");
           goto error;
        }

        o = PyTuple_GetItem(o, 1);
        if (o == NULL) {
           mexPrintf("MPy: Failed to get python item name from dict.\n");
           goto error;
        }

        mxSetFieldByNumber(result, 0, k, MPy_python_to_mat(o));
    }

    Py_DECREF(items);
    
    return result;

 error:
    Py_XDECREF(items);
    Py_CLEAR(result);
    return result;
}




static mxArray *MPy_sequence_to_mat(PyObject *object)
{
     mxArray *result=NULL;
     size_t length = 0;

     length = PySequence_Size(object);
     if (length == -1) {
         mexPrintf("MPy: Failed to get python sequence size.\n");
         return result;
     }

     result = mxCreateCellArray(1,&length);
     
     for (int i = 0; i < length; i++) {
          PyObject *item = PySequence_GetItem(object, i);
          if (! item) {
               mexPrintf("MPy: Failed to get sequence item # %d",i);
          } else {
               mxArray* mx = MPy_python_to_mat(item);
               Py_DECREF(item);
               if (!mx) {
                  mexPrintf("MPy: Failed to encode item # %d",i);
               } else {
                  mxSetCell(result,i,mx);
               }
          }
     }

     return result;
}



static mxArray* MPy_unicode_to_mat(PyObject *object) {
   PyObject *encoded;

   encoded = PyUnicode_AsUTF16String(object);

   size_t dims[2];
   char* buffer;
   Py_ssize_t length;

   mxArray *result;

   PyString_AsStringAndSize(encoded, &buffer, &length);
   /* First character is the encoding indicator */
   length -= sizeof(mxChar);
   buffer += sizeof(mxChar);

   dims[0] = 1;
   dims[1] = length;
   result = mxCreateCharArray(2, dims);

   memcpy(mxGetData(result), buffer, length * sizeof(char));

   Py_XDECREF(encoded);
   
   return result;
}



static mxArray *MPy_python_to_mat(PyObject *object)
{
     void *data;
     mxArray *result = NULL;

     if (strcmp(object->ob_type->tp_name,"numpy.ndarray")==0) {
          result = MPy_ndarray_to_mat(object);
          
     } else if (PyInt_Check(object)) {
          result = MPy_int_to_mat(object);
          
     } else if (PyFloat_Check(object)) {
          result = MPy_float_to_mat(object);
    
     } else if (PyBool_Check(object)) {
          result = MPy_bool_to_mat(object);
          
     } else if (PyString_Check(object)) {
          result = MPy_string_to_mat(object);

     } else if (PyComplex_Check(object)) {
          result = MPy_complex_to_mat(object);
 
     } else if (PyUnicode_Check(object)) {
          result = MPy_unicode_to_mat(object);
  
     } else if (PyMapping_Check(object)) {
          result = MPy_dict_to_mat(object);
          
     } else if (PySequence_Check(object)) {
          result = MPy_sequence_to_mat(object);
      
     } else if (PyObject_Compare(object, Py_None) == 0) {
          result = MPy_none_to_mat(object);
          
     } else {
          char buffer[1024];
          sprintf(buffer,"MPy: Unknown object type: %s\n",object->ob_type->tp_name);
          mexPrintf(buffer);
     }

     return result;
}






  








/*********************************************************
 *  Running Python
 *********************************************************/


/************************
 *  Python Messages
 ************************/

static PyObject* MPy_message(PyObject *self, PyObject *args) {
   const char *what;
   
   if (!PyArg_ParseTuple(args, "s", &what)) return NULL;
   
   printf("%s", what);
   
   return Py_BuildValue("");
}

static PyMethodDef MPy_message_method[] = {
   {"write", MPy_message, METH_VARARGS, "message"},
   {NULL, NULL, 0, NULL}
};

static void MPy_message_inititialize() {
   PyObject *m = Py_InitModule("MPyMessage", MPy_message_method);
   if (m == NULL) return;
   PySys_SetObject((char*) "stdout", m);
   PySys_SetObject((char*) "stderr", m);
}




/***************************
 *  Python Initialization
 ***************************/

static void MPy_initialize_builtins() {
   if (PyDict_GetItemString(MPy_globals, "__builtins__") == NULL) {      
      DEBUG("MPy: Initializing loading builtins...")
      PyObject* builtinMod = PyImport_ImportModule("__builtin__");
      if (builtinMod == NULL || PyDict_SetItemString(MPy_globals, "__builtins__", builtinMod) != 0) {
         Py_XDECREF(MPy_globals);
         mexErrMsgTxt("MPy: failed to load buildins!");
         return;
      }
   }
}


static void MPy_initialize() { 
   
   DEBUG("MPy: Initializing Python...")
   Py_Initialize();
   
   DEBUG("MPy: Initializing Python Messages...")
   MPy_message_inititialize();
   
   DEBUG("MPy: Initializing Python Main...")
   if (!MPy_main) 
      MPy_main = PyImport_AddModule("__main__");
   
   if (!MPy_globals)
      MPy_globals = PyModule_GetDict(MPy_main);
   
   DEBUG("MPy: Initializing import_array...")
   import_array();
   
   DEBUG("MPy: Initializing builtins...")
   MPy_initialize_builtins();

   MPy_status = Initialized;
}


static void MPy_finalize() {  
   DEBUG("MPy: Finalizing Python")
           
   DEBUG("MPy: Finalizing Python, clearing dict...")
   
   if (MPy_globals) PyDict_Clear(MPy_globals);
   Py_CLEAR(MPy_globals);
   
   
   DEBUG("MPy: Finalizing Python, clearing main...")
   Py_CLEAR(MPy_main);
    
   //DEBUG("MPy: Finalizing Python, finalize...")
   //Py_Finalize(); //this causes trouble when restarting using numpy as Py_Finalize is buggy

   MPy_status = Uninitialized;
}



static void MPy_clear() {  
   DEBUG("MPy: Clear Python Dict...")
   
   if (MPy_globals) PyDict_Clear(MPy_globals);
   Py_CLEAR(MPy_globals);
   
   DEBUG("MPy: Clear Python reinitialize...") 
   
   MPy_globals = PyModule_GetDict(MPy_main);
   
   MPy_initialize_builtins();
}



static void MPy_evaluate(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

   if (nrhs != 2 || !mxIsChar(prhs[1])) {
      mexErrMsgTxt("MPy: syntax error! Usage: py('eval', stmt)");
   }
   
   char *stmt = mxArrayToString(prhs[1]);   
   DEBUG("MPy: evaluate:") DEBUG(stmt)
   
   PyObject *o = PyRun_String(stmt, Py_file_input, MPy_globals, MPy_globals);
   mxFree(stmt);
  
   if (o == NULL) {
      PyErr_Print();
      mexErrMsgTxt("MPy: Error while evaluating Python statement");
   }

   if (nlhs > 0) {
      plhs[0] = MPy_python_to_mat(o);
      Py_DECREF(o);

      if (plhs[0] == NULL) {

         mexErrMsgTxt("MPy: Error converting output to matlab");
      }
   }
}


static void MPy_set(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   if (nrhs != 3 || !mxIsChar(prhs[1])) {
      mexErrMsgTxt("MPy: syntax error! Usage: py('set', var_name, var)");
   }
   
   char *var_name = mxArrayToString(prhs[1]);
      
   DEBUG("MPy: set, var_name:") DEBUG(var_name)

   PyObject *var = MPy_mat_to_python(prhs[2]);
   
   if (var) {
      PyDict_SetItemString(MPy_globals, var_name, var);
      Py_DECREF(var);
      mxFree(var_name);
   } else {
      mxFree(var_name);
      mexErrMsgTxt("MPy: error converting expression to python!");
   }
}


static void MPy_get(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   if (nlhs != 1 || nrhs != 2 || !mxIsChar(prhs[1])) {
      mexErrMsgTxt("MPy: syntax error! Usage: var = py('get', expr)");
   }
   
   char *expr = mxArrayToString(prhs[1]);
   PyCodeObject *code = (PyCodeObject*) Py_CompileString(expr, "no_source", Py_eval_input);
   mxFree(expr);
   
   if (code == NULL) {
      PyErr_Print();
      mexErrMsgTxt("MPy: Error compiling expression");
   }
   
   PyObject *o = PyEval_EvalCode(code, MPy_globals, MPy_globals);
   Py_DECREF(code);
   
   if (o == NULL) {
      mexErrMsgTxt("MPy: Error evaluating Python expression");
   }
   
   if (nlhs > 0) {
      plhs[0] = MPy_python_to_mat(o);
      Py_DECREF(o);

      if (plhs[0] == NULL) {
         mexErrMsgTxt("MPy: Error converting output to matlab");
      }
   } else {
      Py_DECREF(o);
   }
}







/************************
 *   Mex Function
 ************************/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

   if (!MPy_status) MPy_initialize();

   //mexAtExit(MPy_finalize);

   if (nrhs == 0 || !mxIsChar(prhs[0])) {
      mexErrMsgTxt("MPy: syntax error. Usage: py(cmd, varargin)");
   }
   
   char *cmd = mxArrayToString(prhs[0]);
   
   if (!strcmp(cmd, "eval")) {
      MPy_evaluate(nlhs, plhs, nrhs, prhs);
      
   } else if (!strcmp(cmd, "set")) {
      MPy_set(nlhs, plhs, nrhs, prhs);
      
   } else if (!strcmp(cmd, "get")) {
      MPy_get(nlhs, plhs, nrhs, prhs);

   } else if (!strcmp(cmd, "clear")) {
      MPy_clear();
      
   } else if (!strcmp(cmd, "exit")) {
      MPy_finalize();
            
   } else {
      mexErrMsgTxt("MPy: unknown command, use eval, set, get, clear or exit");
   }
   
   mxFree(cmd);
}
