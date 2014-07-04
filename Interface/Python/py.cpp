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

static PyObject* MPy_globals;

enum MPy_status_type {Uninitialized = 0, Initialized};
static MPy_status_type MPy_status = Uninitialized;

//#pragma GCC diagnostic ignored "-fpermissive"

#include <iostream>
#define DEBUG(msg) std::cout << msg << std::endl; std::cout.flush();




/*********************************************************
 *    Array Conversions
 *********************************************************/

/*************************
 * Matlab -> Python
 * **********************/

static PyObject* MPy_mat_to_python(const mxArray* mx);


static void MPy_swap_dims(size_t ndims, size_t* dims) {
   if (ndims > 1)
      std::swap(dims[0], dims[1]);
}

static PyObject* MPy_char_to_python(const mxArray* mx)
{   
   DEBUG("MPy_char_to_python")
   
   char* buf; size_t buflen;
   PyObject* result;
   
   if (mxGetM(mx) != 1) {
      mexErrMsgTxt("MPy: String expected to be 1d row of chars in Matlab!\n");
   }

   buflen = mxGetN(mx) + 1;
   buf = (char*) mxCalloc(buflen, sizeof(char)); 
   if (mxGetString(mx, buf, buflen) != 0) {
      mxFree(buf);
      mexErrMsgTxt("MPy: MPy_char_to_python: error when converting to string!\n");
   }
   
   result = PyString_FromStringAndSize(buf, buflen-1);

   mxFree(buf);   

   return result;
}


static PyObject* MPy_numeric_to_python(const mxArray* mx)
{
   DEBUG("MPy_numeric_to_python")
   
   const void *real, *imag;
   size_t n, nd, *dims;
   
   PyObject* result;
   npy_intp npytype = 0;

   real = mxGetData(mx);
   imag = mxGetImagData(mx);
   n    = mxGetNumberOfElements(mx);
   nd   = mxGetNumberOfDimensions(mx);
   dims = (size_t*) mxGetDimensions(mx);
   MPy_swap_dims(nd, dims);
   
   switch (mxGetClassID(mx)) {
      case mxDOUBLE_CLASS:
            npytype = imag ? NPY_CDOUBLE : NPY_DOUBLE;
            if (imag) {
                  result = PyArray_SimpleNew(nd, (npy_intp*) dims, npytype);
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
      default:
         DEBUG("MPy_numeric_to_python: unknown class id:" << mxGetClassID(mx))
         mexErrMsgTxt("MPy: MPy_numeric_to_python, Unexpected class ID.");
   }
   
   DEBUG("MPy_numeric_to_python npytype=" << npytype << " mxclass=" << mxGetClassID(mx))

   result = PyArray_SimpleNew(nd, (npy_intp*) dims, npytype);
   result = (PyObject*) PyArray_GETCONTIGUOUS((PyArrayObject*) result);
   
   memcpy(PyArray_DATA((PyArrayObject*) result), real, n * PyArray_ITEMSIZE((PyArrayObject*)result));
   
   return result;
}

static PyArray_Descr* MPy_struct_dtype(const mxArray* mx)
{
   DEBUG("MPy_struct_dtype")
   
   int nfields;
   
   PyArray_Descr *descr;
   PyObject *list;

   nfields = mxGetNumberOfFields(mx);
   list = PyList_New(nfields);
   for (int i = 0; i < nfields; i++) {
         PyList_SetItem(list, i, Py_BuildValue("(s, s)", mxGetFieldNameByNumber(mx, i), "object"));
   }
   if (PyArray_DescrConverter(list, &descr) != NPY_SUCCEED) {
         PyErr_Print();
         mexErrMsgTxt("MPy: MPy_struct_dtype: Failed to determine dtype of struct array.");
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
   PyObject *result, **itemptr;

   n       = mxGetNumberOfElements(mx);
   nfields = mxGetNumberOfFields(mx);
   nd      = mxGetNumberOfDimensions(mx);
   dims    = (size_t*) mxGetDimensions(mx);
   MPy_swap_dims(nd, dims);

   
   dtype = MPy_struct_dtype(mx);
   
   result = PyArray_SimpleNewFromDescr(nd, (npy_intp*) dims, dtype);

   itemptr = (PyObject**) (PyArrayObject**) PyArray_DATA((PyArrayObject*) result);
   for (size_t i = 0; i < n; i++)
         for (int k = 0; k < nfields; k++)
            *itemptr++ = MPy_mat_to_python(mxGetFieldByNumber(mx, i, k));
   
   return result;
}

static PyObject* MPy_cell_to_python(const mxArray* mx)
{
   DEBUG("MPy_cell_to_python")
   
   const mxArray *cell;
   size_t n, nd, *dims;
   
   PyObject *result, **itemptr, *value;

   
   n    = mxGetNumberOfElements(mx);
   nd   = mxGetNumberOfDimensions(mx);
   dims = (size_t*) mxGetDimensions(mx);
   MPy_swap_dims(nd, dims); 

   result = PyArray_SimpleNew(nd, (npy_intp*)  dims, NPY_OBJECT);
   itemptr = (PyObject**) (PyArrayObject**) PyArray_DATA((PyArrayObject*) result);
   
   for (size_t i = 0; i < n; i++) {
      cell = mxGetCell(mx, i);
      if (cell)
         value = MPy_mat_to_python(cell);
      else {
         value = Py_None;
         Py_INCREF(value);
      }
      *itemptr++ = value;
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
            return MPy_char_to_python(mx);
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
            return MPy_numeric_to_python(mx);
      case mxSTRUCT_CLASS:
            return MPy_struct_to_python(mx);
      case mxCELL_CLASS:
            return MPy_cell_to_python(mx);
      case mxFUNCTION_CLASS:
            return MPy_function_to_python(mx);
      case mxUNKNOWN_CLASS:
            mexErrMsgTxt("MPy: MPy_to_python: unknown class to convert to Python.");
      default:
            mexErrMsgTxt("MPy: MPy_to_python: unexpected type to convert to Python.");
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


static mxArray *MPy_ndstruct_to_mat(PyObject *object, PyObject *dtype, size_t ndims, const size_t *dims, size_t length)
{
   DEBUG("MPy_ndstruct_to_mat");
     
   PyObject *contiguous = NULL;
   PyObject *fields = NULL;
   int nfields, npy_type;
   char **keys = NULL;
   PyObject **stringized = NULL;
   PyObject *names = NULL;
   PyObject **data;
   
   mxArray *result;
   
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
   stringized = (PyObject**) mxCalloc(nfields, sizeof(PyObject *));
   
   for (int i=0; i<nfields; i++) {
      PyObject *key = PyTuple_GetItem(names, i);
      if (!key) {
         mexPrintf("MPy: Failed to get key from tuple for field # %d\n",i);
         goto exit;
      }
      stringized[i] = PyObject_Str(key);
      if (! PyString_Check(stringized[i])) {
         mexPrintf("MPy: Key # %d couldn't be stringized\n",i);
         goto exit;
      }
      keys[i] = PyString_AsString(stringized[i]);
   }


   for (int i=0; i<nfields; i++) {
      PyObject *subtype;
      PyObject *item;

      item = PyMapping_GetItemString(fields, keys[i]);
      if (! PyTuple_Check(item)) {
         mexPrintf("MPy: Did not get a tuple item for field %s\n", (char*) keys[i]);
         Py_XDECREF(item);
         goto exit;
      }
      
      subtype = PyTuple_GetItem(item,0);
      if (! subtype) {
         mexPrintf("MPy: Failed to get sub-type from tuple for field %s\n",keys[i]);
         Py_XDECREF(subtype); Py_XDECREF(item);
         goto exit;
      }
      if (! PyArray_DescrCheck(subtype)) {
         mexPrintf("MPy: Subtype for field %s was not an array descriptor.\nIt was type %s\n",keys[i],subtype->ob_type->tp_name);
         Py_XDECREF(subtype); Py_XDECREF(item);
         goto exit;
      }
      
      npy_type = MPy_get_dtype(subtype);
      Py_XDECREF(subtype); Py_XDECREF(item);
      
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
   if (keys) mxFree(keys);
   if (stringized) {
      for (int i = 0; i<nfields; i++) {
         Py_XDECREF(stringized[i]);
      }
      mxFree(stringized);
   }
   Py_XDECREF(contiguous);

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
   
   size_t* dims = NULL;
   mxArray *result = NULL;
   void *data = NULL;
   int length=1;
   
   long npy_type = NPY_NOTYPE;
   char buffer[1024];

   /* Get the array dimensions */

   shape = PyObject_GetAttrString(object, "shape");
   if (! shape) {
      mexPrintf("MPy: numpy.ndarray did not have a shape\n");
      goto error;
   }
   if (! PyTuple_Check(shape)) {
      mexPrintf("MPy: numpy.ndarray did not have a tuple shape\n");
      goto error;
   }
   
   ndims = PyTuple_Size(shape);
   dims = (size_t *) mxMalloc(sizeof(size_t)*ndims);
   for(Py_ssize_t i=0; i< ndims; i++) {
      PyObject *obDim = PyTuple_GetItem(shape,i);
      if (! PyInt_Check(obDim)) {
         mexPrintf("MPy: Array shape is not a tuple of integers\n");
         goto error;
      }
      dims[i] = (size_t) PyInt_AsLong(obDim);
      length *= dims[i];
   }
   
   // correct natlab vs python arrays
   MPy_swap_dims(ndims, dims);

   DEBUG("MPy_ndarray_to_mat: dtype")
   
   /*
   * Handle different dtypes
   */
   dtype = PyObject_GetAttrString(object, "dtype");
   if (! dtype) {
      mexPrintf("MPy: numpy.ndarray did not have a dtype\n");
      goto error;
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
            goto error;
         }
         result = mxCreateCellArray(ndims, dims);
         if (! result) {
            strcpy(buffer,"MPy: mxCreateCellArray returned NULL for an array of dimension ");
            strcpy(MPy_dimstrcpy(buffer+strlen(buffer), ndims, dims),"\n");
            mexPrintf(buffer);
            goto error;
         }
         for (int idx = 0; idx < length; idx++) {
            PyObject *py_value;
            mxArray *mx_value;
            py_value = (*iter->ob_type->tp_iternext)(iter);
            if (! py_value) {
               mexPrintf("MPy: numpy.ndarray's __iter__'s next function failed\n");
               goto error;
            }
            Py_INCREF(py_value);
            mx_value = MPy_python_to_mat(py_value);
            Py_DECREF(py_value);
            if (! mx_value) {
               mexPrintf("MPy: ndarray_to_mxarray failed while trying to translate a subobject to an mxarray\n");
               goto error;
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
         goto error;
      }
   }
   
   DEBUG("MPy_ndarray_to_mat: exit")
   
   goto exit;
error:
   if (dims) mxFree(dims);
   result=mxCreateNumericMatrix(1,1,mxLOGICAL_CLASS, mxREAL);
   data = mxGetData(result);
   *((bool *)data) = 0;
     
exit:
   if (dims)
      mxFree(dims);
   if (shape)
      Py_DECREF(shape);
   if (dtype)
      Py_DECREF(dtype);
   if (dtype_num)
      Py_DECREF(dtype_num);
   if (iter)
      Py_DECREF(iter);

   DEBUG("MPy_ndarray_to_mat: return")
   
   return result;
}



static mxArray *MPy_dict_to_mat(PyObject *object)
{
     int nItems,pass;
     
     PyObject *items=NULL;
     mxArray *result=NULL;
     char **keys = NULL;
     PyObject **stringized = NULL;
     mxArray *mxValue;


     nItems = PyMapping_Length(object);

     keys = (char **)mxCalloc(nItems, sizeof(char *));
     stringized = (PyObject **)mxCalloc(nItems, sizeof(PyObject *));
     if (! keys) {
          mexPrintf("MPy: Not enough memory for keys array\n");
          goto exit;
     }

     items = PyMapping_Items(object);

     if (! items) {
          mexPrintf("MPy: PyMapping_Items returned null\n");
          goto exit;
     }
     if (! PyList_Check(items)) {
          mexPrintf("MPy: Failed to get list of items from dictionary.\n");
          goto exit;
     }

     /*
     ** Loop twice - once to get the keys from the tuples and make the array
     **              a second time to fill in the values for the array
     */
     for (pass = 0; pass < 2; pass++) {

          if (pass == 1) {
               /*
               ** Create the array on the second pass (after we have the keys)
               */
               result = mxCreateStructMatrix(1,1,nItems,(const char**) keys);
          }
          for (int i=0; i<nItems; i++) {
               PyObject *kv;
               PyObject *item = PyList_GetItem(items,i);
               if (! PyTuple_Check(item)) {
                    mexPrintf("MPy: Item # %d was not a tuple.\n",i);
                    goto exit;
               } 
               kv = PyTuple_GetItem(item, pass);
               if (! kv) {
                    mexPrintf("MPy: Failed to get %s for item # %d\n",pass?"key":"value",i);
                    goto exit;
               }
               switch(pass) {
               case 0:  /* Key */
                    stringized[i] = PyObject_Str(kv);
                    if (! PyString_Check(stringized[i])) {
                         mexPrintf("MPy: Key # %d did not stringize\n",i);
                         goto exit;
                    }
                    keys[i] = PyString_AsString(stringized[i]);
                    break;
               case 1: /* value */
                    mxValue = MPy_python_to_mat(kv);
                    if (! mxValue) {
                         mexPrintf("MPy: Value for %s was untranslatable\n",keys[i]);
                         goto exit;
                    }
                    mxSetField(result, 0, keys[i], mxValue);
                    break;
               }
          }
     }
  exit:
     if (items)
          Py_DECREF(items);
     if (keys)
          mxFree(keys);
     if (stringized) {
          for (int i=0; i<nItems; i++) {
               if (stringized[i])
                    Py_DECREF(stringized[i]);
          }
          mxFree(stringized);
     }

     return result;
}



static mxArray *MPy_sequence_to_mat(PyObject *object)
{
     mxArray *result=NULL;
     int nItems = 0;

     nItems = PySequence_Size(object);
     if (nItems == -1) {
          mexPrintf("MPy: Failed to get sequence size.\n");
          goto exit;
     }

     result = mxCreateCellMatrix(1,nItems);
     for (int i=0; i<nItems; i++) {
          PyObject *item = PySequence_GetItem(object, i);
          mxArray *mxItem;
          if (! item) {
               mexPrintf("MPy: Failed to get sequence item # %d",i);
               goto exit;
          }
          mxItem = MPy_python_to_mat(item);
          Py_DECREF(item);
          if (! mxItem) {
               mexPrintf("MPy: Failed to encode item # %d",i);
               goto exit;
          }
          mxSetCell(result,i,mxItem);
     }

  exit:
     return result;
}






static mxArray *MPy_python_to_mat(PyObject *object)
{
     void *data;
     mxArray *result=0;

     if (strcmp(object->ob_type->tp_name,"numpy.ndarray")==0) {
          result = MPy_ndarray_to_mat(object);
          
     } else if (PyInt_Check(object)) {
          result=mxCreateNumericMatrix(1,1,mxINT32_CLASS, mxREAL);
          data = mxGetData(result);
          *((long *)data) = PyInt_AsLong(object);
          
     } else if (PyFloat_Check(object)) {
          result=mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS, mxREAL);
          data = mxGetData(result);
          *((double *)data) = PyFloat_AS_DOUBLE(object);
          
     } else if (!PyObject_IsTrue(object)) {
          result=mxCreateNumericMatrix(1,1,mxLOGICAL_CLASS, mxREAL);
          data = mxGetData(result);
          *((bool *)data) = 0;
          
     } else if (PyString_Check(object)) {
          Py_ssize_t length;
          PyObject *encoded;
          char *ptr;

          encoded = PyObject_CallMethod(object, "encode", "s,s", "utf_16", "replace");
          if (! encoded) {
               mexPrintf("MPy: Failed to encode string in python_to_mxarray\n");
               return NULL;
          }
          PyString_AsStringAndSize(encoded,&ptr, &length);
          /* First character is the encoding indicator */
          length-=sizeof(mxChar);
          ptr += sizeof(mxChar);

          result=mxCreateNumericMatrix(1,length/sizeof(mxChar),mxCHAR_CLASS, mxREAL);
          data = mxGetData(result);
          memcpy(data,ptr,length);
          
          Py_XDECREF(encoded);
          
     } else if (PyUnicode_Check(object)) {
          Py_ssize_t length;
          PyObject *encoded;
          char *ptr;

          encoded = PyUnicode_AsUTF16String(object);
          PyString_AsStringAndSize(encoded,&ptr, &length);
          /* First character is the encoding indicator */
          length-=sizeof(mxChar);
          ptr += sizeof(mxChar);

          result=mxCreateNumericMatrix(1,length/sizeof(mxChar),mxCHAR_CLASS, mxREAL);
          data = mxGetData(result);
          memcpy(data,ptr,length);
          
          Py_XDECREF(encoded);
          
          
     } else if (PyMapping_Check(object)) {
          result = MPy_dict_to_mat(object);
          
     } else if (PySequence_Check(object)) {
          result = MPy_sequence_to_mat(object);
          
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
    PySys_SetObject("stdout", m);
    PySys_SetObject("stderr", m);
}


/***************************
 *  Python Initialization
 ***************************/


void MPy_initialize() { 
   
   DEBUG("MPy: Initializing Python...")

   Py_Initialize();
   
   MPy_message_inititialize();
   
   MPy_globals = PyModule_GetDict(PyImport_AddModule("__main__"));
   
   import_array();

   MPy_status = Initialized;
}


void MPy_finalize() {  
   DEBUG("MPy: Finalizing Python...")
   
   if (MPy_globals) PyDict_Clear(MPy_globals);
   
   DEBUG("MPy: Finalizing Python...")
   
   //Py_XDECREF(MPy_globals);
   
   DEBUG("MPy: Finalizing Python...")
   
   Py_Finalize();
   
   DEBUG("MPy: Finalizing Python...")

   MPy_status = Uninitialized;
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

   if (nlhs > 0) plhs[0] = MPy_python_to_mat(o);
   
   Py_DECREF(o);
}


static void MPy_set(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   if (nrhs != 3 || !mxIsChar(prhs[1])) {
      mexErrMsgTxt("MPy: syntax error! Usage: py('set', var_name, var)");
   }
   
   char *var_name = mxArrayToString(prhs[1]);
      
   DEBUG("MPy: set, var_name:") DEBUG(var_name)

   PyObject *var = MPy_mat_to_python(prhs[2]);
   
   PyDict_SetItemString(MPy_globals, var_name, var);
   
   //Py_DECREF(var);
   mxFree(var_name);
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
   
   // PyObject *locals = PyDict_New();
   PyObject *o = PyEval_EvalCode(code, MPy_globals, MPy_globals);
   // Py_DECREF(locals);
   if (o == NULL) {
      mexErrMsgTxt("MPy: Error evaluating Python expression");
   }
   Py_DECREF(code);
   
   if (nlhs > 0) plhs[0] = MPy_python_to_mat(o);
   //else //display ?
      
   
   if (plhs[0] == NULL) {
      mexErrMsgTxt("MPy: Error converting to MATLAB variable");
   }
}









/************************
 *   Mex Function
 ************************/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

   if (!MPy_status) MPy_initialize();

   mexAtExit(MPy_finalize);

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

   } else if (!strcmp(cmd, "exit")) {
      MPy_finalize();
            
   } else {
      mexErrMsgTxt("MPy: unknown command");
   }
   
   mxFree(cmd);
}
