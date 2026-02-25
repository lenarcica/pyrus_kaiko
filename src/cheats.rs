//////////////////////////////////////////
/// cheats.rs
///
/// Alan Lenarcic 2025
///
/// GPL_v2 license, consider coding your own if open source is insufficient
/// (this is possibly not current because PyArrow Rust library makes changes to FFI
/// often)
///
/// When certain Pyarrow functionality seems to be broken, we implement some calls to
/// FFI_ArrowArrayStream.
///
/// My understanding is that at somepoint pyarrow might have a better method for converting
/// Recordbatches into an exportable PyResult.
use arrow::record_batch::{RecordBatch, RecordBatchIterator};
use pyo3::{PyObject, PyResult};
//use pyo3::{pymodule, types::PyModule, PyResult, Bound};
use pyo3::types::PyTuple;
use arrow::record_batch::RecordBatchReader;
use pyo3::Python;
use arrow::ffi_stream::{FFI_ArrowArrayStream};
//use arrow::ffi_stream::{ArrowArrayStreamReader, FFI_ArrowArrayStream};
use pyo3::ffi::Py_uintptr_t;
use pyo3::types::PyAnyMethods;

// A Cheat implementation, 
//   Merges the RecordBatch and RecordBatch reader implementations in
//   arrow::pyarrow crate, because there is some sort of "gated" restriction
//   Possibly the gate is triggered by trying to use pyo3->arrow in addition.
//
// In the end, we get to see a little of the constructions that are being called
// to allow the Record batch to be streamed into the Arrow interpretter and then bound
// as a tuple to the python level.  From the user level, it might be good
// to have some proficiency with this level, which seems to mimic the C/C++
// interactions with the arrow layer, since the rust libraries may want to 
// change their interfaces
pub fn recordbatch_to_pyarrow(rb:&RecordBatch, py:Python) -> PyResult<PyObject> {
  // Debugging comments, though apparently this does not crash so far.
  //println!("recordbatch_to_pyarrow -- Generating iterator.");
  let reader = RecordBatchIterator::new(vec![Ok(rb.clone())],rb.schema());
  let reader: Box<dyn RecordBatchReader + Send> = Box::new(reader);
  //println!("recordbatch_to_pyarrow -- Generating mut stream.");
  let mut stream = FFI_ArrowArrayStream::new(reader);
  let stream_ptr = (&mut stream) as *mut FFI_ArrowArrayStream;
  let module = py.import_bound("pyarrow")?;
  let class = module.getattr("RecordBatchReader")?;
  //println!("recordbatch_to_pyarrow -- args, PyTuple, new bound.");
  let args = PyTuple::new_bound(py, [stream_ptr as Py_uintptr_t]);
  let reader = class.call_method1("_import_from_c", args)?;

  let py_reader = PyObject::from(reader);
  py_reader.call_method0(py, "read_next_batch")
}
// From pyarrow though this doesn't import
/********
 // How to cheat get PyArrow to RecordBatch
 impl ToPyArrow for RecordBatch {
    fn to_pyarrow(&self, py: Python) -> PyResult<PyObject> {
        // Workaround apache/arrow#37669 by returning RecordBatchIterator
        let reader = RecordBatchIterator::new(vec![Ok(self.clone())], self.schema());
        let reader: Box<dyn RecordBatchReader + Send> = Box::new(reader);
        let py_reader = reader.into_pyarrow(py)?;
        py_reader.call_method0(py, "read_next_batch")
    }
}
// Note, this rabitt holes, since reader.into_pyarrow is also not implemented.
// There appears still to be fun in call_method0 which happens after the pyarrow conversion
 *************/
