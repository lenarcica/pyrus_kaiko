#![allow(unused_imports)] 

use pyo3::{Python, Bound, types::{IntoPyDict, PyAnyMethods, PyDict}};
use numpy::{Element};
use numpy::{PyArrayDescr, get_array_module, PyArrayDescrMethods, PyUntypedArrayMethods};
//use std::ffi::{CString};
//yuse pyo3::ffi::{c_str};
use pyo3::prelude::*;
use std::str;
use std::convert::Into;

// Why numpy requires CStr seems like strangeness?
use std::ffi::CStr;
use std::ffi::CString;

//use pyo3_polars::{PyDataFrame}
// get_array_module
//use numpy::{ dtype_bound
//use numpy::ndarray::{Array};

const NSIDEL:usize = 5;

pub type TTime = i64;
pub type TPrice = f64;
pub type TSide = [ u8; NSIDEL ];

pub struct TTimePrice{ pub ttime: TTime, pub tprice: TPrice }

pub struct TSideS{ tside:  TSide }


pub struct TTmPrSd{pub ttime: TTime, pub tprice: TPrice, pub tside_s: TSideS}

impl Clone for TTimePrice {
  fn clone(&self) -> Self {
    TTimePrice{ttime:(*self).ttime as TTime, tprice: (*self).tprice as TPrice}
  }
}

impl Clone for TTmPrSd {
  fn clone(&self) -> Self {
    TTmPrSd{ ttime: (*self).ttime as TTime, tprice: (*self).tprice as TPrice, tside_s: (*self).tside_s.clone() as TSideS }
  }
}
impl Clone for TSideS {
  fn clone(&self) -> Self {
    TSideS{tside: (*self).tside}
  }	  
}

/***
impl From<&str> for TSide {
  fn from(x:&str) -> TSide {
    let mut aout:TSide =  [b'\0',b'\0',b'\0',b'\0',b'\0']; //[0u8;NSIDEL] as TSide; //['\0','\0','\0','\0','\0'];
    for itt in 0..4 {
      //aout[itt as usize] = (x.chars().nth(itt as usize)) as u8;
	  aout[itt as usize] = x.as_bytes()[itt as usize];
    }
    //aout[<usize as Into<i64>>::into(NSIDEL-1)] = '\0';
	aout[4] = b'\0'; 
    aout
  }
}
***/

impl From<&str> for TSideS {
  fn from(x:&str) -> TSideS {
    let mut aout:TSide =  [b'\0',b'\0',b'\0',b'\0',b'\0']; //[0u8;NSIDEL]; //
	const NSM1:usize = NSIDEL-1;
    for itt in 0..4 {
     // aout[itt as usize] = (x.chars().nth(itt as usize)).unwrap();
	  aout[itt as usize] = x.as_bytes()[itt as usize];
    }
    aout[NSM1] = b'\0';
    TSideS{tside:aout}
  }
}
// "Element" Trait implementation for the TTimePrice
//
// To automate trnaslation of TTimePrice to a Numpy RecArray
//  we need Element trait impelmented.
// Note this challenge requires knowledge of how to construct 'dtype'
//  properly inside running python.  This is a challenge
//  because python "np.dtype" can be finicky and loose with type instructions.
//
// Note rust-numpy says "Currently, only integer/float/complex/object types are supported."
//   This means that this will not be usable to make string arrays?  This is unfortunate
//   since a lot of numpy output will want to show string output.
unsafe impl Element for TTimePrice {
  const IS_COPY: bool = true;  
  // const IS_COPY:-- Flag that indicates whether this type is trivially copyable.
  // Required methods
  fn get_dtype(py: Python<'_>) -> Bound<'_, PyArrayDescr> {
    //let outV = Python::with_gil(|py| {
    //let np = py.import_bound("numpy")?;
    let locals = [("np", get_array_module(py).unwrap())].into_py_dict(py).expect("Hey Unwrap np please");
    //let locals = PyDict::new(py);
    //let np = py.import_bound("numpy")?;
    //locals.set_item("np", np)?;
    
    // What String type does "eval_bound" require.
    //  Aparently py.eval used to take CString input
    //  However, in this case, it appears that i_str can be
    //  a &str and it works okay.
    //
    //let IStr = CString::new("dtypev=[('time','datetime64[ns]'),('price',np.float64)]; np.zeros(0, dtype=dtypev).dtype".to_string()).unwrap();
    // Worried about where
    let i_str = "import numpy as np; np.dtype([('time','datetime64[ns]'),('price',np.float64)])";
    println!("Element Implementation for TTimePrice, We wont need it but: i_str={}",
      i_str);

    //let IString = CString::new(format!());
    //as Option<&Bound<'py, PyDict>>, Some(locals), None as Option<&Bound<'_, PyDict>> )
    //.expect("Python Numpy FAIL Happened reason")
    //let pydtype = py
    //    .eval_bound(i_str, Some(&locals), None).unwrap()
    //    .downcast_into::<PyArrayDescr>().unwrap();
    let untyped_s: CString = CString::new("np.dtype([('time','datetime64[ns]'),('price','float64')])").expect("This CString should have allocated"); // Make this CStr?
    let untyped_stmt: &CStr = untyped_s.as_c_str();
    let result = py
       .eval(untyped_stmt, Some(locals).as_ref(), None);
    //  .map_err(|e| {
    //    e.print_and_set_sys_last_vars(py);
    //  });
    //let pydtype:<'_,PyArrayDescr> = result.downcast_into::<PyArrayDescr>().unwrap();
    //return Bound::new(py, )?;
    return result.unwrap().downcast_into::<PyArrayDescr>().unwrap();
    //});
    //return outV; 
    }
	    
    // Weird: because _py lifetime is unimportant, it asks that we call "_py" by underscore
    fn clone_ref(&self, _py: Python<'_>) -> Self {
	return TTimePrice{ ttime:(*self).ttime,tprice: (*self).tprice };	
    }
    /***********
    fn get_dtype<'py>(py: Python<'py>) -> &'py PyArrayDescr {
	  self.get_dtype_bound(py).unbind();
	
    }
	***********/
	// Weird,  -- by request, py is recorded by underscore "_py"
    fn vec_from_slice(_py: Python<'_>, slc: &[Self]) -> Vec<Self> { 
	  return slc.to_vec();
	}
	/*************
    fn array_from_view<D>(
        py: Python<'_>, view: ArrayView<'_, Self, D>,
      ) -> Array<Self, D>
    where D: Dimension {
	  let arr:Array<Self,D> =Array::new();
      for ii in [1..D] {
        arr[ii] = TTimePrice(view[ii].0, view[ii].1)
      }		  
      return(arr);   
	}
	***************/
}

// A Test can we create DTypes with Side included?
//   For fun we will move side to left.
unsafe impl Element for TTmPrSd {
  const IS_COPY: bool = true;  
  fn get_dtype(py: Python<'_>) -> Bound<'_, PyArrayDescr> {
    let locals = [("np", get_array_module(py).unwrap())].into_py_dict(py).expect("Hey TTmPrSd want error?");
      //let i_str = "import numpy as np; np.dtype([('time','datetime64[ns]'),('price',np.float64),('side','|S5')])";
      let untyped_s: CString = CString::new("np.dtype([('time','datetime64[ns]'),('price','float64'),('side','|S5')])").expect("Hey this is an easy string to allocate");
      let result = py.eval(untyped_s.as_c_str(), Some(locals).as_ref(), None);
      return result.unwrap().downcast_into::<PyArrayDescr>().unwrap();
    }
    // Weird: because _py lifetime is unimportant, it asks that we call "_py" by underscore
    fn clone_ref(&self, _py: Python<'_>) -> Self {
      return TTmPrSd{ttime:(*self).ttime,tprice: (*self).tprice,tside_s: (*self).tside_s.clone()};	
    }
    fn vec_from_slice(_py: Python<'_>, slc: &[Self]) -> Vec<Self> { 
      return slc.to_vec();
    }
}  

/*******
 // We thought of creating vecotr type of Time Price but not a good idea.
pub fn dynVtimePrice(n: u32) -> vTimePrice {
   let nV = vTimePrice::new();
   if n > 0 {
     for ii in [0..n] {
	   vTimePrice::irow.push(ii, ii * 100.0);	 
     }
   }
   return vTimePrice;
}
*******/

/*******************************************************
 ** Gemini Example: multiple array python usage
 ** Preumably Py_BuildValue can in C API can save user from makign more dynamic type
 
 // Create a structured dtype
    PyArray_Descr *dtype = NULL;
    PyObject *dtype_tuple = Py_BuildValue("[(s, s), (s, s)]", "name", "S10", "age", "i4");
    if (PyArray_DescrConverter(dtype_tuple, &dtype) < 0) {
        PyErr_Print();
        return 1;
    }
    Py_DECREF(dtype_tuple);
	npy_intp dims[1] = {2};
    PyObject *array = PyArray_SimpleNewFromDescr(1, dims, dtype);
    if (array == NULL) {
        PyErr_Print();
        return 1;
    }

    // Access and modify elements
    void *data = PyArray_DATA((PyArrayObject *)array);
    *((char *)(data + 0 * PyArray_STRIDE((PyArrayObject *)array, 0))) = 'A';
    *((int *)(data + 0 * PyArray_STRIDE((PyArrayObject *)array, 1))) = 20;
    *((char *)(data + 1 * PyArray_STRIDE((PyArrayObject *)array, 0))) = 'B';
    *((int *)(data + 1 * PyArray_STRIDE((PyArrayObject *)array, 1))) = 25;

    // Print the recarray
    PyObject_Print(array, stdout, 0);
    Py_DECREF(array);
	*******************************************************************/
