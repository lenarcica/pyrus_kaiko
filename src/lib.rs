///~///////////////////////////////////////////////////////////////////////
///  lib.rs -- pyrus_kaiko package
///
///  Alan Lenarcic 2025-08-05
///
///  The following is an interaction between python and rust to process Kaiko full orderbook files
///  The results should be extractable for a full orderbook engine written in "pyrus_marketbook"
///
///  Kaiko "full orderbook" or "fob" files have a particular format showing an update of what
///  levels have new quantity posted. 
///
///  However, there can be some data feed issues and orders whose closing price goes unreported
///  We demonstrate a level verfication algorithm that can help clean that data and so it can
///  be used correctly in orderbook inspection algorithms.
// We don't need to configure hpsort for this library
//mod hpsort;

// mod verify_levs for verifying Kaiko/Simulated FOB bids against "top of book" (Best-Bid/Best-Ask)
// prices
mod verify_levs;
use verify_levs::{verify_buy, verify_sell, slow_verify_price, sort_v_pt, parallel_verify_buy, parallel_verify_sell};

//
mod b2v_struct;
mod sip_struct;
mod sip_algo;
// Verifying Levels is helpful when an exchange reporting to Kaiko has it's full order book feed
// depart from its trade or top-of-book feed.


use std::sync::Arc;
use std::vec::Vec;
use std::io::{BufReader, BufRead};
use std::fs::File;

use flate2::bufread::GzDecoder;

use pyo3::Python;
use pyo3::{pymodule, types::PyModule, PyObject, PyResult, Bound};

//rust-numpy package: Grab Rust
use numpy::ndarray::{ArrayD, ArrayViewD,ArrayViewMutD};
use numpy::{IntoPyArray, PyArrayDyn, PyArray1, PyReadonlyArrayDyn, PyArrayMethods};


//use pyo3_arrow::error::PyArrowResult;
//use pyo3_arrow::PyArray as arrowPyArray;
use pyo3_arrow::error::{PyArrowResult,PyArrowError};
use pyo3::exceptions::{PyValueError};
use pyo3_arrow::PyTable as arrowPyTable;
//use arrow_select;
use arrow::datatypes::Field;
use arrow::datatypes::Schema;
use arrow::error::ArrowError;
use arrow::record_batch::{RecordBatch};
use arrow_schema::SchemaRef;
use arrow::array::ArrayData;
use arrow::buffer::Buffer;
use pyo3_arrow::PyRecordBatch;
use pyo3_arrow::PyArray as arrowPyArray;
use arrow::array::{StringArray, Decimal64Array, Int8Array, UInt8Array, UInt64Array, Int64Array, Array};
use arrow::array::{Int32Array, Int16Array, Float64Array, BinaryArray,TimestampNanosecondArray};
use arrow::datatypes::{DataType, TimeUnit};

//use crate::hpsort::idx_sort_trio;
//use arrow_array::{Decimal64Array, Array};
//use crate::cheats::recordbatch_to_pyarrow;
// A Copied Cheat Cheat implementation: recordbatch conversion to pyarrow
//
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
//
/*
pub fn recordbatch_to_pyarrow(rb:&RecordBatch, py:Python) -> PyResult<PyObject> {
  // Debugging comments, though apparently this does not crash so far.
  //println!("recordbatch_to_pyarrow -- Generating iterator.");
  let reader = RecordBatchIterator::new(vec![Ok(rb.clone())],rb.schema());
  let reader: Box<dyn RecordBatchReader + Send> = Box::new(reader);
  //println!("recordbatch_to_pyarrow -- Generating mut stream.");
  let mut stream = FFI_ArrowArrayStream::new(reader);
  let stream_ptr = (&mut stream) as *mut FFI_ArrowArrayStream;
  //let module = py.import_bound("pyarrow")?;
  let module = py.import("pyarrow")?;
  let class = module.getattr("RecordBatchReader")?;
  //println!("recordbatch_to_pyarrow -- args, PyTuple, new bound.");
  //let args = PyTuple::new_bound(py, [stream_ptr as Py_uintptr_t]);
  let args = PyTuple::new(py, [stream_ptr as Py_uintptr_t]);
  let reader = class.call_method1("_import_from_c", args)?;

  let py_reader = PyObject::from(reader);
  py_reader.call_method0(py, "read_next_batch")
}
*/
///~////////////////////////////
/// get_n_prices_per_line:
///
/// Used in breaking down TOB files into prices, quantities
fn get_n_prices_per_line(fs: String) -> (u32, u8,u8) {
  let mut n_line: u32 = 0; 
  let mut nbar:u8 = 0;
  let mut in_bracket:u8 = 0;
  let mut n_comma = 0;
  let mut p_scale:u16 = 0;
  let mut q_scale:u16 = 0;
  let mut onp_scale:i16 = -1;
  let mut onq_scale:i16 = -1;
  let mut pq01 = 0;
  for c in fs.chars() {
    match c {
      '|'|';' => {
         if in_bracket > 0 {
           println!(" --- Error on reading fs = {}, we are on c={} but in_bracket={}",
           fs, c, in_bracket);
         }
         nbar = nbar + 1; },
      '[' => {
        if nbar < 2 {
          println!(" --- Error on reading fs = {}, we are on c={}, but we already have [ without 2 bars!",fs,c);
        } else if nbar > 4 {
          println!(" --- Error on reading fs = {} we have nbar = {}, but we shouldn't jump past 2 Buy/Sell states?",
           fs, c);
        }
        in_bracket = in_bracket + 1;  n_comma = 0; 
        pq01 = 0; onp_scale = -1; onq_scale = -1; },
      ']' => {
         if (n_comma != 1) && (in_bracket != 1) {
           println!(" -- Error on reading fs ={}, we are on c={} for outBracket but n_comma={}",
                fs, c, n_comma);
         } 
         in_bracket = in_bracket -1;
         if in_bracket == 1 { n_line = n_line + 1; 
           if onq_scale > (q_scale as i16) { q_scale = onq_scale as u16; }
         }
         onp_scale = -1; onq_scale = -1;
      },
      ',' => {
        if in_bracket == 2 {
          if onp_scale > (p_scale as i16) { p_scale = onp_scale as u16; }
          onp_scale = -1; pq01 = 1;
        }
        n_comma = n_comma + 1;
      },
      '.' => {
         if pq01 == 0 {  
           onp_scale = 0;
         } else {
           onq_scale = 0;
         }
      },
      '0'|'1'|'2'|'3'|'4'|'5'|'6'|'7'|'8'|'9' => {
         if (pq01 == 0) && (onp_scale >= 0) && (in_bracket == 2) { onp_scale = onp_scale + 1;
         } else if (pq01 == 1) && (onq_scale >=0) && (in_bracket==2) { onq_scale = onq_scale+1;
         }
      }
      _ =>  {}
    };
  }
  return (n_line, p_scale as u8, q_scale as u8);
}


#[pymodule]
fn pyrus_kaiko<'py>(_py: Python<'py>, m: &Bound<'py, PyModule>) -> PyResult<()> {


/***********************************************************************
from pyrus_kaiko import pyrus_kaiko; import polars as pl; import numpy as np;
import pandas as pd; import pyarrow as pa;
mypt = pyrus_kaiko.kaiko_fob("c:/Users/alanj/Dropbox/py/pyrus_kaiko/data/fakedata001.gz", 0);
mypd = mypt.to_pandas();
mypd = mypd.sort_values(by=['bs01','price','time']).reset_index(inplace=False, drop=True)
mclose = mypd['time'].shift(-1, fill_value=pd.NaT);
mclose[(mypd['bs01'] != mypd['bs01'].shift(-1, fill_value=2)) | (mypd['price'] != mypd['price'].shift(-1,fill_value=-1))] = pd.NaT;
mypd['close'] = mclose
mypd['time'] = mypd['time'].astype('datetime64[ns]')
mypd['close'] = mypd['close'].astype('datetime64[ns]')
mside = np.repeat('b', len(mypd))
mside[mypd['bs01'] == 1] = 's';
mypd['side'] = pd.Series(mside).astype(pd.StringDtype())

from pyrus_kaiko import pyrus_kaiko;
import pyrus_kaiko.kaiko as pko;
gzLoc = "c:/Users/alanj/Dropbox/py/pyrus_kaiko/data/fakedata001.gz";
pko.KaikoReadFOB(gzLoc=gzLoc, verbose = 0);
 ***********************************************************************/
///~///////////////////////////////////////////////////////////////////////////////
///  Format of Kaiko tob files
///
///           timestamp;    type;                          asks;              bids;
///     171353923043023;       u;[[100.0,0.0032],[101.0,0.004]];       [[99.9,10]];
///     171353923053932;       u;               [[392.0,0.004]];                [];
///     171353934039329;       s;[[10.09........
///
///  That is, timestamp will be large integer nanoseconds from Epoch (1970-01-01 UTC)
///  Type is either u for an update, s for a snapshot.  Asks come first, bids second
///  And nesting occurs between brackets, with a "bar" being a semicolon
///  Spaces are likely minimized
#[pyfn(m)]
#[pyo3(name = "kaiko_fob")]
fn kaiko_fob<'py>(py:Python<'py>, file_path: String,
    verbose: u8) -> PyResult<PyObject> {
  //let fpath: &str = file_path.into()?;  // Filename
  //let cow_fpath:Cow<String> = file_path.to_string_lossy();
  let fpath = file_path;
  //let cow_fpath:Cow<String> = file_path.to_cow(py).expect("Come On Cow this is a string!");
  //let fpath:String = cow_fpath.into_owned();  // Convert to a COW and then String 
  //let file_path = "your_file.gz"; // Replace with your .gz file path
  let mut total_out_lines = 0;
  let mut price_scale = 0;
  let mut qty_scale = 0;
  let precision:u8 = 15; // Need to know max precision of Decimal64, but this is pretty good
  if verbose >= 1 {
    println!("kaiko_fob_extract() -- Initiating for file {}", fpath);
  }
  // Ammusing Rust Scope trick.  Make things close by putting them in inner scope
  {
    let opened_file = File::open(fpath.clone())?;  // File Open
    let decoder = GzDecoder::new(BufReader::new(opened_file)); // GzDecoder
    // Create a BufReader around the GzDecoder for line-by-line reading
    let mut reader = BufReader::new(decoder);
    if verbose >= 1{
      println!("kaiko_fob_extract() -- we have opened first reader. ");
    }
    //let first_line = reader.lines().next(); // ignore first line
    reader.read_line(&mut String::new())?;
    for line_result in reader.lines() {
      let line = line_result?;
      let (nlines, p_scale, q_scale) = get_n_prices_per_line(line);
      total_out_lines = total_out_lines + nlines;
      if p_scale > price_scale { 
        price_scale = p_scale;
      }
      if q_scale > qty_scale {
        qty_scale = q_scale; 
      }
    }
  }

  if verbose >= 1 {
    println!("kaiko_fob_extract(): Note we read {} lines, with a price_scale={}, qty_scale={} ",
       total_out_lines, price_scale, qty_scale);
  }
  // If we can determine the max scale for every price and quantity we can create a Decmial64
  //  array.  Unfortunately, this does not guarantee all of the same prices will have same number
  //  of decimals so we still need to read number of decimals every time.
  if verbose >= 1 {
    println!("kaiko_fob_extract() -- Creating time/bs01/price/qty vector length {}.", total_out_lines);
  }
  let mut col00_time_vec = vec![0 as u64;total_out_lines as usize];
  let mut col01_us01_vec = vec![0 as u8;total_out_lines as usize];
  let mut col02_bs01_vec = vec![0 as u8;total_out_lines as usize];
  let mut col03_price_vec = vec![0 as u64;total_out_lines as usize];
  let mut col04_qty_vec = vec![0 as u64;total_out_lines as usize];

  let mut onii = 0;
  let radix:u32 = 10;
  {
    let opened_file_2 = File::open(fpath)?;  // File Open
    let decoder = GzDecoder::new(BufReader::new(opened_file_2)); // GzDecoder
    // Create a BufReader around the GzDecoder for line-by-line reading
    let mut reader = BufReader::new(decoder);
    reader.read_line(&mut String::new())?;
    for line_result in reader.lines() {
      let line = line_result?;
      if verbose >= 2 {
        println!("  -- on line: {}", &line);
      }
      // Always re-initiate these things on loop
      let mut on_time: u64 = 0;
      let mut on_price:u64 = 0;
      let mut on_qty:u64 = 0;
      let mut n_bars:u32 = 0;
      let mut in_bracket:u32 = 0;
      let mut n_comma = 0; let mut p_loc:i16 = -1; let mut q_loc:i16 = -1;
      let mut on_us01: u8 = 2;
      for c in line.chars() {
        // Fun little match, calculate qty
        //println!("on c = {}, on_price={}, on_qty={}, n_bars={}, n_comma={}, in_bracket={}",&c, &on_price, &on_qty, &n_bars, &n_comma, &in_bracket);
        match c {
          '|'|';' => { n_bars = n_bars + 1; },
          '[' => { in_bracket = in_bracket + 1; n_comma = 0; },
          ']' => {
            if ((n_bars==2) || (n_bars==3)) && (in_bracket==2) {
              col00_time_vec[onii] = on_time;
              col01_us01_vec[onii] = on_us01;
              col02_bs01_vec[onii] = if n_bars == 2 { 1 as u8 } else { 0 };  // asks, then bids
              while (q_loc >= 0) && ( (q_loc as u8) < qty_scale) {
                on_qty = on_qty * 10;  q_loc = q_loc + 1; // rescale quantity if we saw less decimals than qty_scale
              }
              while (p_loc >= 0) && ( (p_loc as u8) < price_scale) {
               on_price = on_price*10; p_loc = p_loc + 1; // Rescale price if it has less decimals than normal scale
              }  
              col03_price_vec[onii] = on_price;
              col04_qty_vec[onii] = on_qty;
              onii = onii + 1; p_loc = -1; on_price=0; on_qty = 0; p_loc = -1; q_loc = -1; n_comma = 0;
            }
            in_bracket = in_bracket -1;
          },
          'S'|'s' => { if (n_bars == 1) && (on_us01==2) {on_us01=1;}},
          'U'|'u' => { if (n_bars == 1) && (on_us01==2) {on_us01=0;}},
          '0'|'1'|'2'|'3'|'4'|'5'|'6'|'7'|'8'|'9' => {
            if n_bars == 0 {
              // c.to_digit(radix).unwrap()
              // Time parsing, I guess we could multiple 00:00:00.0000 into numbers
              // That said, Kaiko default format is Nanoseconds from UTC 1970-01-01
              on_time = on_time*(10 as u64) + (c.to_digit(radix).expect("Hey Trying to unwrap on_time digit to a digit") as u64);
            } else if n_bars >= 2 {
              if n_comma == 0 {
                on_price = on_price*10 + (c.to_digit(radix).expect("Hey trying to unwrap on_price digit to a digit") as u64);
                if p_loc >= 0 { p_loc = p_loc + 1; }
              } else if n_comma == 1 {
                on_qty = on_qty*10 + (c.to_digit(radix).expect("Hey Trying to unwrap on_qty to a digit") as u64);
                if q_loc >= 0 { q_loc = q_loc + 1; }
              }
            }
          },
          '.'=> {
            // Found our period, hope there is none;  Are we updating price or qty now?
            if (n_bars >= 2) && (n_comma == 0) {
              p_loc = 0;
            } else {
              q_loc = 0;
            }
          },
          ','=> {
             if (n_comma == 0) && (n_bars >= 2) {
               //while p_loc >= 0 {
               //  // Of course we could create decimal type, but we don't feel like it.  Dividing by
               //  // 10 probably does it, though accuracy might go away;
               //  on_price = on_price * 0.1;  p_loc = p_loc -1;
               //}
               n_comma = n_comma + 1;
             }
          },
          _ => {
            // Default do nothing
          }
        }
      }  // End loop for c in line
    } // End loop thru lines
  } // End lifetime of file object #2

  // Lets write the columns to an arrow table
  //let col00_time = Arc::new(UInt64Array::from(col00_time_vec)) as _;
  //let col01_us01 = Arc::new(UInt8Array::from(col01_us01_vec)) as _;
  //let col02_bs01 = Arc::new(UInt8Array::from(col02_bs01_vec)) as _;
  let col00_time = Arc::new(UInt64Array::from(col00_time_vec));
  let col01_us01 = Arc::new(UInt8Array::from(col01_us01_vec));
  let col02_bs01 = Arc::new(UInt8Array::from(col02_bs01_vec));
  let buffer_price = Buffer::from_slice_ref(&col03_price_vec);

  // Attempting to write Decimals
  let array_price_data = ArrayData::builder(arrow::datatypes::DataType::Decimal64(precision, price_scale as i8))
    .len(total_out_lines as usize).add_buffer(buffer_price).build().unwrap();
  let buffer_qty = Buffer::from_slice_ref(&col04_qty_vec);
  let array_qty_data = ArrayData::builder(arrow::datatypes::DataType::Decimal64(precision, qty_scale as i8))
    .len(total_out_lines as usize).add_buffer(buffer_qty).build().unwrap();
  let col04_qty_array = Arc::new(Decimal64Array::from(array_qty_data));
  let col03_price_array = Arc::new(Decimal64Array::from(array_price_data));

  if verbose >= 1 {
    println!("lib.rs->kaiko_tob -- completing generation of TOB arrow array");
  }
  let schema = SchemaRef::new(Schema::new(vec![
     Field::new("time", arrow::datatypes::DataType::UInt64, false),
     Field::new("us01", arrow::datatypes::DataType::UInt8, false),
     Field::new("bs01", arrow::datatypes::DataType::UInt8, false),
     Field::new("price", arrow::datatypes::DataType::Decimal64(precision, price_scale as i8), false),
     Field::new("qty", arrow::datatypes::DataType::Decimal64(precision, qty_scale as i8), false)]));

    let batch = RecordBatch::try_new(schema, 
      vec![col00_time,col01_us01,col02_bs01,col03_price_array,col04_qty_array]).expect(" -- ERROR Why didnt RecordBatch create?");
    //let batch = RecordBatch::try_from_iter([("time", col00_time), 
    //  ("us01", col01_us01),("bs01",col02_bs01),
    //  ("price",col03_price_array),
    //  ("qty", col04_qty_array)
    //]).unwrap();
    // Yes we tried a buncy of ways to convert the batch to pyarrow
    //   But there appears to be a "gated" off warning when trying
    //   to directly access to_pyarrow() trait.
    //return arrowPyTable::from(batch);
    //return batch.to_pyarrow();
    if verbose >= 1 {
      println!("lib.rs->kaiko_tob -- exporting to recordbatch.");
    }
    //return recordbatch_to_pyarrow(&batch, py);
    let pybatch:PyRecordBatch = PyRecordBatch::new(batch);
    return pybatch.to_pyarrow(py);
  }

fn dec_print(vc: &mut Vec<char>, m_dec: u8, idn: u64) {
  if idn == 0 {
    vc.push('0'); return;
  }
  let binding = idn.to_string();
  let an = binding.chars();
  let l_an = an.clone().count();
  let tgt:usize = ((l_an as i64) - (m_dec as i64)).try_into().expect("dec_print, this cannot have m_dec less than l_an");
  for (index, c) in an.enumerate() {
    if tgt == index {
      vc.push('.');
    }
    vc.push(c);
  }
}
/************************************************************************
 ## Test kaiko_make_u_fob
from pyrus_marketbook import pyrus_marketbook; from pyrus_marketbook import ob;
from pyrus_kaiko import pyrus_kaiko;
import pyrus_kaiko.kaiko as pkk;
from pyrus_kaiko.kaiko import genASnapShot;
import numpy as np; import pandas as pd;
import pyarrow as pa;
mktbook = pkk.Gen_Fake_Mktbook(nMidpoints=10000, nGen=20000, nVen=10, verbose=1,
  stTime = '2025-07-25 10:00', deltaSec=1, startMid = 1000)
fobbook = pkk.Kaiko_FOB_format_mktbook(mktbook, nSnaps=20, verbose =1, 
  printEvery=1000, dec_price=2, dec_qty=0) 
pkk.Kaiko_FOB_save_fobbook(fobbook, saveLoc = "c:/users/alanj/Dropbox/py/pyrus_kaiko/data/data_save_001.gz") 
 ************************************************************************/
#[pyfn(m)]
#[pyo3(name = "kaiko_make_u_fob")]
fn kaiko_make_u_fob<'py>(py:Python<'py>, price_dec: u8, qty_dec:u8, verbose: u8,
    print_every:u64,
    topen: PyReadonlyArrayDyn<'py, i64>, 
    bs01: PyReadonlyArrayDyn<'py, u8>,
    idprice: PyReadonlyArrayDyn<'py,u64>,idqty: PyReadonlyArrayDyn<'py,u64>
  ) -> PyResult<PyObject> {
  // Here we only create "update" scripts.
  // It is easy enough to make "snapshots", as those are a minority of data, and can append/sort
  // anywhere
  let mut total_unique_times:u64 = 0;
  let topen = topen.as_array();
  let bs01 = bs01.as_array();
  let idqty = idqty.as_array();
  let idprice = idprice.as_array();
  if topen.len() <= 1 {
    println!("kaiko_make_fob, nothing to do for length 1 or less.");
  }
  for ii in 1..topen.len() {
    if topen[ii] != topen[ii-1] {
      total_unique_times = total_unique_times + 1;
    }
  }
  let total_unique_times:usize = total_unique_times as usize;
  if verbose >= 1 {
    println!("kaiko_make_fob, we have to create material");
  }
  let mut ov_time:Vec<i64> = Vec::new();
  // Let's be lazy vector of 'u'!
  let mut ov_type:Vec<char> = Vec::new(); 
  let mut ov_asks:Vec<String> = Vec::new();
  let mut ov_bids:Vec<String> = Vec::new();
  let precision:u8 = 15; // Need to know max precision of Decimal64, but this is pretty good
  if verbose >= 1 {
    println!("kaiko_make_fob_t() -- Initiating for table of length {}", total_unique_times);
  }
  let mut on_time = topen[0];
  //let mut on_w_ii: u64 = 0;
  let mut askv: Vec<char> = Vec::new();
  let mut bidv: Vec<char> = Vec::new();
  let tlen = topen.len();
  askv.push('['); bidv.push('[');
  for  ii in 0..tlen {
    if (verbose >= 2) && ( (0 as u64) == (ii as u64) % (print_every as u64)) {
      println!("{}/{}  -- ont={}, askv_len={}, bidv_len={}",
        ii, tlen, on_time, askv.len(), bidv.len());
    }
    if (on_time != topen[ii]) || (ii == tlen-1)  {
      ov_time.push(on_time); ov_type.push('u');
      askv.push(']'); bidv.push(']');
      ov_asks.push(askv.clone().into_iter().collect()); ov_bids.push(bidv.clone().into_iter().collect());
      askv.clear(); askv.push('['); bidv.clear(); bidv.push('[');
      on_time = topen[ii];
    }
    if bs01[ii] == 0 {
      if bidv.len() > 1 { bidv.push(','); }  
      bidv.push('[');
      dec_print(&mut bidv,price_dec, idprice[ii]);
      bidv.push(',');
      dec_print(&mut bidv, qty_dec, idqty[ii]);
      bidv.push(']');
    } else {
      if askv.len() > 1 { askv.push(','); }  
      askv.push('[');
      dec_print(&mut askv,price_dec, idprice[ii]);
      askv.push(',');
      dec_print(&mut askv, qty_dec, idqty[ii]);
      askv.push(']');
    }
  }
  if verbose >= 1 {
    println!(" -- At End of algorithm, Length ov_time={}, ov_type={}, ov_asks={}, ov_bids={}",
      ov_time.len(), ov_type.len(), ov_asks.len(), ov_bids.len());
  }
  let col00_time = Arc::new(Int64Array::from(ov_time));
  // vec![Some(&my_string)]
  //  This seems to be cheat way, turn Vec<char> in silly way to Vec<String>, and see Arrow to
  //  work.
  let ov_type_string_vec: Vec<String> = ov_type 
        .into_iter()
        .map(|c| c.to_string())
        .collect();
  //let ov_type_string_array:StringArray = StringArray::from(vec![ov_type_string.as_str()]);
  let col01_type = Arc::new(StringArray::from(ov_type_string_vec));
  let col02_asks = Arc::new(StringArray::from(ov_asks));
  let col03_bids = Arc::new(StringArray::from(ov_bids));

  if verbose >= 1 {
    println!("lib.rs->kaiko_make+fpb -- completing generation of FOB arrow array");
    println!("  -- note col00_time is length {}, col01_type has length {}, col02/03 asks/bids={},{}", 
      col00_time.len(), col01_type.len(), col02_asks.len(), col03_bids.len());
  }
  let schema = SchemaRef::new(Schema::new(vec![
     Field::new("timestamp", arrow::datatypes::DataType::Int64, false),
     Field::new("type", arrow::datatypes::DataType::Utf8, false),
     Field::new("bids", arrow::datatypes::DataType::Utf8, false),
     Field::new("asks", arrow::datatypes::DataType::Utf8, false)]));

    let batch = RecordBatch::try_new(schema, 
      vec![col00_time,col01_type,col02_asks,col03_bids]).expect("Make Fob -- ERROR Why didnt RecordBatch create?");
    //let batch = RecordBatch::try_from_iter([("time", col00_time), 
    //  ("us01", col01_us01),("bs01",col02_bs01),
    //  ("price",col03_price_array),
    //  ("qty", col04_qty_array)
    //]).unwrap();
    // Yes we tried a buncy of ways to convert the batch to pyarrow
    //   But there appears to be a "gated" off warning when trying
    //   to directly access to_pyarrow() trait.
    //return arrowPyTable::from(batch);
    //return batch.to_pyarrow();
    if verbose >= 1 {
      println!("lib.rs->kaiko_make_fob -- exporting to recordbatch.");
    }
    //return recordbatch_to_pyarrow(&batch, py);
    let pybatch:PyRecordBatch = PyRecordBatch::new(batch);
    return pybatch.to_pyarrow(py);
}

// Typically this can be reduced to an integer problem with effort.
//
//  Note, the major challenges came from trying to unwrap and get n.
//  Length unwrapping is hard!  Presumably if t0 comes out as a null this has a panic.
// AS JOIN a critical sort join for time series data
//   Its apparently pretty easy to implement with just above.
/**********************************
from pyhr import pyhr; import numpy as np;
n0 = 100; n1 = 200;  TMax = 10000000
t0 = np.sort(np.random.randint(0,TMax,n0)).astype(np.int64)
t1 = np.sort(np.random.randint(0,TMax,n1)).astype(np.int64)
## t0/t1 will be not unique but both sorted.
out0 = pyhr.asjoin_py(t0,t1)

import pandas as pd;
pdJ = pd.DataFrame({"t0":t0, "t1N": t1[out0], "out0":out0});
np.max(t1[t1 <= t0[0]]) ## Quickly check and verify that pdJ mach is correct.

******************************************/
  #[pyfn(m)]
  #[pyo3(name = "asjoin_py")]
  fn asjoin_py<'py>(py:Python<'py>, t0: PyReadonlyArrayDyn<'py,i64>,
    t1: PyReadonlyArrayDyn<'py,i64>)
    -> Bound<'py, PyArray1<i64>> {
    //-> Bound<'py, numpy::PyArray<i64, Dim<IxDynImpl>>> {
    let t0 = t0.as_array();
    let t1 = t1.as_array();
    let n:u32 = t0.len().try_into().unwrap();

    //let iDyn = IxDynImpl(&dVec);
    //let DimVec = dVec.into_dimension();
    //Bound<'py, PyArrayDyn<f64>> 
    let out_array0 = PyArray1::<i64>::zeros(py, (n) as usize, false);
    //let mut o0 = t0.clone();  // Output will be copy same size as t0;
    // Of course, we'd love index output to be usize, but we negatives/NANS for
    //   non matches
    let mut on1:i64 = 0;
    let mut on0:i64 = 0;
    while on0 < (t0.len() as i64) {
      while (on1 < (t1.len()-1) as i64) && (t1[(on1+1) as usize] < t0[on0 as usize]) {
        on1 = on1 + 1;
      }
      if t1[on1 as usize] <= t0[on0 as usize] {
        // Why do these not work but the zeros bound below works?
        unsafe { *(out_array0.data().wrapping_add(on0 as usize)) = on1; }
        //unsafe{ *(out_qrray0.data() + on0) = on1 }; // Fails now?
      } else {
        unsafe{ *(out_array0.data().wrapping_add(on0 as usize))=-1; }
        //unsafe{ *(out_array0.data() + on0) = on1};
      }
      on0 = on0+1;
    }
    out_array0
  }

///~////////////////////////////////////////////////////////////////////////////
/// unsorted_asjoin_py
///
/// Yes, we cannot always assume data for as of join is sorted.  To simply this use case
///  and help in cases where it is difficult or annoying to save sorted data, this generates
///  an out_array that gives an as of join based upon sorted Rust data
  #[pyfn(m)]
  #[pyo3(name = "unsorted_asjoin_py")]
  fn unsorted_asjoin_py<'py>(py:Python<'py>, t0: PyReadonlyArrayDyn<'py,i64>,
    t1: PyReadonlyArrayDyn<'py,i64>)
    -> Bound<'py, PyArray1<i64>> {
    //-> Bound<'py, numpy::PyArray<i64, Dim<IxDynImpl>>> {
    let t0 = t0.as_array();
    let mut t_ord = vec![0 as usize; t0.len()];
    for ii in 0..t0.len() { t_ord[ii] = ii; }
    t_ord.sort_by_key(|&jj| &t0[jj]); 
    let t1 = t1.as_array();
    let n:u32 = t0.len().try_into().unwrap();

    let out_array0 = PyArray1::<i64>::zeros(py, (n) as usize, false);
    //let mut o0 = t0.clone();  // Output will be copy same size as t0;
    // Of course, we'd love index output to be usize, but we negatives/NANS for
    //   non matches
    let mut on1:i64 = 0;
    let mut on0:i64 = 0;
    while on0 < (t0.len() as i64) {
      while (on1 < (t1.len()-1) as i64) && (t1[(on1+1) as usize] < t0[t_ord[on0 as usize]]) {
        on1 = on1 + 1;
      }
      if t1[on1 as usize] <= t0[t_ord[on0 as usize]] {
        // Why do these not work but the zeros bound below works?
        unsafe { *(out_array0.data().wrapping_add(t_ord[on0 as usize])) = on1; }
        //unsafe{ *(out_qrray0.data() + on0) = on1 }; // Fails now?
      } else {
        unsafe{ *(out_array0.data().wrapping_add(t_ord[on0 as usize]))=-1; }
        //unsafe{ *(out_array0.data() + on0) = on1};
      }
      on0 = on0+1;
    }
    out_array0
  }


/*****************
unsafe impl Element for char {
  const IS_COPY: bool = false;  
  fn get_dtype_bound(py: Python<'_>) -> Bound<'_, PyArrayDescr> {
    let locals = [("np", get_array_module(py).unwrap())].into_py_dict_bound(py);
      //let i_str = "import numpy as np; np.dtype([('time','datetime64[ns]'),('price',np.float64),('side','|S5')])";
      let result = py
        .eval_bound("np.dtype('x','|S1')", Some(&locals), None);
      return result.unwrap().downcast_into::<PyArrayDescr>().unwrap();
    }
    // Weird: because _py lifetime is unimportant, it asks that we call "_py" by underscore
    fn clone_ref(&self, _py: Python<'_>) -> Self {
      return self;
    }
    fn vec_from_slice(_py: Python<'_>, slc: &[Self]) -> Vec<Self> { 
      return slc.to_vec();
    }
}  
***********************/
/********************************
Scratch function for testing PyArrow as it sends data into Rust.
Consider pyarrow tables, test_pa will try and read some data from the table as needed.

import numpy as np; from pyrus_test import pyrus_test; import pyarrow as pa;
import pandas as pd;
v_s = np.array(['b','b','b','b','s','b','b','b','b','b','s','s','b','b','b']).astype('|S1')
v_t = np.asarray([0,1,2,3,4,5,6,7,8,9,10,11,12,13,-1]).astype(np.int64);
v_a = np.zeros(len(v_t), np.int64) + 1;
v_p = np.asarray([9,9,1,1,1,1,9,9,9,9,9,8,8,1,1]).astype(np.float64);
regorder = [6,7,0,1,2,3,8,9,10,11,13,12,5,4]
othorder = [2,3,4,5,13,12,0,1,6,7,8,9,12,11]
df = pd.DataFrame({'t':v_t, 's':v_s, 'p':v_p,'a':v_a,'t':v_t});
patb = pa.table(df)
pyrus_test.test_pa(patb);

df.iloc[othorder]
 ********************************/
#[pyfn(m)]
#[pyo3(name = "test_pa")]
fn test_pa_py<'py>(_py:Python<'py>, intable: arrowPyTable) { 
  let (inner_vec_record_batch, _inner_schema): (Vec<RecordBatch>, SchemaRef) = intable.into_inner();
  let inner_record_batch_0:RecordBatch= inner_vec_record_batch[0].clone();

  let vin_01_side = inner_record_batch_0.column(1).clone();
  let vin_01_side = vin_01_side.as_any().downcast_ref::<BinaryArray>().unwrap();
  //let vin_01_side = vin_01_side.into_data() as Vec<u32>;
  println!("len vin_01_side = {}, vin_01_side[0] hard to get.", vin_01_side.len());
  let v00_ch: char = vin_01_side.value(0)[0].try_into().unwrap();
  println!(" v00_ch is {}", v00_ch);
  let v00_ch4: char = vin_01_side.value(4)[0].try_into().unwrap(); 
  println!(" v00_ch4 is {}", v00_ch4);
}

 ///~///////////////////////////////////////////
 /// verify prices algorithm, a Rust accelerated algorithm for multi-thread validation
 ///
 ///  Operation of algorithm.  (assume we are validating bids using a midpoint or best ask "maximum
 ///  allowed price")
 ///
 ///  This algorithm is single threaded (one can trivially parallelize this by running multiple
 ///  independent price validations at a time.
 ///
 ///  We order the observed bids in order of lowest price, earliest time to highest price, highest
 ///  time.
 ///
 ///  We start at lowest observed best ask price (minus a desired delta called "pdiff", which allows
 ///   us to investigate for "all prices that come within 10 cents of best ask" or "all prices
 ///   that eventually exceed best ask by at least 2 cents")
 ///
 ///   We look at all "crossing states" where the best ask goes into, out of, or rests for a time
 ///   at the barier price of interest.  In general this is a number c[p] much less than all BBO
 ///   price change events.
 ///
 ///   We test all bids at this price level to see if they are overlapped with a crossing state.
 ///     This is only the prices at a given level NB[p] is a subset of NB (Num Bids)
 ///     This test is only roughly NB[p]*c[p] versus some sort of other iterative test.
 ///
 ///   We then move the crossing states price 1 tick upwards.  The new set of crossing states can
 ///   be determined from the previous list by eliminating any dips that don't exist and adding
 ///   dips that just come up to this level.  This update is length c[p] + NBO[p]
 #[pyfn(m)]
 #[pyo3(name = "pyrus_verify_prices")]
 pub fn pyrus_verify_prices( py: Python, pdiff: f64, order_records: arrowPyTable, 
     nbbo_records: arrowPyTable, eqstate:bool, sameside:bool, 
     verbose:i8, doslow:i8,timetouch:bool)
   -> PyArrowResult<PyObject> {

   let (order_vec_record_batch, order_schema): (Vec<RecordBatch>, SchemaRef) = order_records.into_inner();
   let order_schema = order_schema.clone();
   if order_schema.fields().len() < 5 {
     println!("pyrus_verify_prices: Order Schema of only length {} supplied.", order_schema.fields().len());
     //return Err(String::from("Data is insufficient size"));
     return Err(PyArrowError::from(ArrowError::InvalidArgumentError(
                "Order Schema does not have 4 elements.".to_string()
            )));
   }
   let order_vec_record_batch_0 = order_vec_record_batch[0].clone();
   let vstr = format!("pyrus_verify_prices({}): ", verbose);
   if verbose >= 3 { println!("\n\n--------------------------------------------------------------------------------------------"); }
   if verbose >= 1 { println!("{} -- Initiate. pdiff={}", vstr, pdiff); }

   macro_rules!col_to_aref{
     ( $rbname: expr, $ii: expr, $a_type:ty, $v_type: ty, $str_t:expr) => {
        $rbname.column($ii).as_any().downcast_ref::<$a_type>()
           .ok_or_else(|| PyValueError::new_err(format!("Failed to extract Batch column {}:{}",$ii, $str_t)))?.values() as &[$v_type]
     }
   }
   // To fix this we have to convince people to send i8 0/1 as preferred array for this.
   // I suppose I can invent other way to copy, but still probably involves a heap allocation for
   // when the type is BinaryArray/Character.
   let vin_side = order_vec_record_batch_0.column(0).clone();
   let vin_side = vin_side.as_any().downcast_ref::<BinaryArray>().unwrap();
   let vstr = format!("pyrus_verify_prices(vrb={},n_v={})", verbose, vin_side.len());
   if verbose >= 1 { println!("{} -- we have vinside ", vstr);} 
   let mut vside:Vec<i8> = vec![-1 as i8; vin_side.len() as usize];
   for ii in 0..vin_side.len() {
     vside[ii] = if 'b' == vin_side.value(ii)[0].try_into().unwrap() { 0 } else {1};
   }
   let vin_p = col_to_aref!(order_vec_record_batch_0,1,Float64Array, f64, "price");
   if verbose >= 1 { println!("{} We have vin_p now, [{},{}...]", vstr,
        if vin_p.len() > 0 { vin_p[0] } else {-100.0},
        if vin_p.len() > 1 { vin_p[1] } else {-100.0}); }
   let vin_ot:&[i64] = match  order_schema.fields().get(3).expect("We confirmed this is length 4+").data_type()
      {
       DataType::Timestamp(TimeUnit::Nanosecond,_) =>  order_vec_record_batch_0.column(3).as_any().downcast_ref::<TimestampNanosecondArray>()
          .ok_or(ArrowError::ParseError(
             "Failed to downcast to TimestampNanosecondArray".to_string(),
          ))?.values(),
       _  => col_to_aref!(order_vec_record_batch_0,3, Int64Array, i64, "open") 
      };
 
   let vin_ct:&[i64] = match order_schema.fields().get(4).expect("i64 from timestamp col 4").data_type() {
       DataType::Timestamp(TimeUnit::Nanosecond,_) =>  order_vec_record_batch_0.column(4).as_any().downcast_ref::<TimestampNanosecondArray>()
          .ok_or(ArrowError::ParseError(
             "Failed to downcast to TimestampNanosecondArray".to_string(),
          ))?.values(),
       _  => col_to_aref!(order_vec_record_batch_0, 4, Int64Array, i64, "close") 
   };
   if verbose >= 1 { println!("{}:  We have collected, side, price, open, close.", vstr); }
   
   let (nbbo_record_batch, nbbo_schema): (Vec<RecordBatch>, SchemaRef) = nbbo_records.into_inner();
   let nbbo_schema = nbbo_schema.clone();
   let nbbo_record_batch_0 = nbbo_record_batch[0].clone();

   if nbbo_schema.fields().len() < 3 {
     println!("pyrus_verify_prices: NBBO Schema of only length {} supplied.", nbbo_schema.fields().len());
     //return Err(String::from("Data is insufficient size"));
     return Err(PyArrowError::from(ArrowError::InvalidArgumentError(
                "NBBO Schema does not have 3 elements.".to_string()
            )));
   }
   /***
   We get by with improved Reference call since we don't strictly need vectors to do this.
   macro_rules!nbbo_col_to_avec{
     ( $ii: expr, $a_type:ty, $v_type: ty, $str_t:expr) => {
        Vec::<$v_type>::from(nbbo_record_batch_0.column($ii).as_any().downcast_ref::<$a_type>()
           .ok_or_else(|| PyValueError::new_err(format!("Failed to extract NBBO column {}:{}",$ii, $str_t)))?.values() as &[$v_type])
     }
   }
   ***/
   let v_nbbot:&[i64] = match  nbbo_schema.fields().get(0).expect("We confirmed this nbbo table is length 3+").data_type()
      {
       DataType::Timestamp(TimeUnit::Nanosecond,_) =>  nbbo_record_batch_0.column(0).as_any().downcast_ref::<TimestampNanosecondArray>()
          .ok_or(ArrowError::ParseError(
             "Failed to downcast to TimestampNanosecondArray".to_string(),
          ))?.values(),
       _  => col_to_aref!(nbbo_record_batch_0, 0, Int64Array, i64, "time") 
      };
   let vstr = format!("pyrus_verify_prices(vrb={},n_v={},n_b={}", verbose, vin_p.len(), v_nbbot.len());
   if verbose >= 1 { println!("{} we have v_nbbot and are completing nbbo extraction.", vstr); }
   // We can likely assume
   let v_nbb:&[f64] = col_to_aref!(nbbo_record_batch_0, 1,Float64Array, f64, "nbb");
   let v_nbo:&[f64] = col_to_aref!(nbbo_record_batch_0, 2,Float64Array, f64, "nbo");
   if verbose >= 1 { println!("{}  -- initiate verify_buy.", vstr); }
   let mut vec_out:Vec<u64> = vec![v_nbbot.len() as u64; vin_p.len()];  
   if vin_p.len() <= 0 {
   } else if doslow == 0 {
     slow_verify_price(&mut vec_out, pdiff, &vside[..], vin_p, 
                          vin_ot, vin_ct,
                          if sameside {&v_nbb  } else { v_nbo }, v_nbbot,
                          eqstate, verbose, timetouch);
   } else if doslow == 1 {
     let st_v = sort_v_pt(&vside, vin_p, vin_ot);
     let confirm_v_buy = verify_buy(&mut vec_out, &st_v, pdiff, &vside[..], vin_p, 
                          vin_ot, vin_ct,
                          if sameside { v_nbb  } else { v_nbo }, v_nbbot,
                          eqstate, verbose, timetouch);
     if verbose >= 1 { println!("{} -- We received confirm_v_buy={}", vstr, confirm_v_buy); }
     if verbose >= 3 { println!("\n\n---------------------------------------------------------------------------------------"); } 
     if verbose >= 1 { println!("{} -- Now inititate verify_sell.", vstr); }
     let confirm_v_sell =  verify_sell(&mut vec_out, &st_v, pdiff, &vside[..], vin_p,
                          vin_ot,  vin_ct,  
                          if sameside { v_nbo  } else { v_nbb }, v_nbbot,
                          eqstate, verbose, timetouch);
     if verbose >= 1 { println!("{} -- We receive confirm_v_sell = {}", vstr, confirm_v_sell); }
   } else {
     let st_v = sort_v_pt(&vside, vin_p, vin_ot);
     let confirm_v_buy = parallel_verify_buy(&mut vec_out, &st_v, pdiff, &vside[..], vin_p, 
                          vin_ot, vin_ct,
                          if sameside { v_nbb  } else { v_nbo }, v_nbbot,
                          eqstate, verbose, timetouch);
     if verbose >= 1 { println!("{} -- We received confirm_v_buy={}", vstr, confirm_v_buy); }
     if verbose >= 3 { println!("\n\n---------------------------------------------------------------------------------------"); } 
     if verbose >= 1 { println!("{} -- Now inititate verify_sell.", vstr); }
     let confirm_v_sell =  parallel_verify_sell(&mut vec_out, &st_v, pdiff, &vside[..], vin_p,
                          vin_ot,  vin_ct,  
                          if sameside { v_nbo  } else { v_nbb }, v_nbbot,
                          eqstate, verbose, timetouch);
     if verbose >= 1 { println!("{} -- We receive confirm_v_sell = {}", vstr, confirm_v_sell); }
   }
   //v_buy.append(&mut v_sell);
   //return PyArray::from_vec_bound(py,v_buy);
   let uint64_array: UInt64Array = vec_out.into(); // Direct conversion
   let uint64_field = Field::new("keep0kill1", DataType::UInt64, false);
   //let int8_schema = Schema::new(vec![int8_field]);
   //let int8_schemaref:SchemaRef = Arc::new(int8_schema);
   Ok(arrowPyArray::new(Arc::new(uint64_array), Arc::new(uint64_field)).to_arro3(py)?.into())
 }

  ///~///////////////////////////////////////////////////////////
  /// sip_arrow[] Rust function
  ///
  /// Demonstrates a basic SIP NBBO capability, joining multiple simultanous Best-bid/ask streams
  ///  together
  ///
  /// Taking Table input possible but difficult in this case.
  ///  It seems that Arrow is preferable for column based output, given
  ///  That it's interface will join disparate columns together
  ///
  ///
  ///  But for safety/instruction, it is useful to maintain pyo3 code
  ///  in the lib.rs file alone, and all linked libaries/functions be vanila rust.
  #[pyfn(m)]
  #[pyo3(name = "sip_arrow")]
  fn sip_arrow<'py>(py:Python<'py>, n_v: u8, vt_in: PyReadonlyArrayDyn<'py,i64>,
    vbp_in: PyReadonlyArrayDyn<'py,f64>, vsp_in: PyReadonlyArrayDyn<'py,f64>,
    vbq_in: PyReadonlyArrayDyn<'py,f64>, vsq_in: PyReadonlyArrayDyn<'py,f64>,
    vbv_in: PyReadonlyArrayDyn<'py,u8>, vsv_in: PyReadonlyArrayDyn<'py,u8>,
    verbose: u8) -> PyResult<PyObject> {
    let vt_in = vt_in.as_array();
    let vbp_in = vbp_in.as_array();
    let vbq_in = vbq_in.as_array();
    let vbv_in = vbv_in.as_array();
    let vsp_in = vsp_in.as_array();
    let vsq_in = vsq_in.as_array();
    let vsv_in = vsv_in.as_array();
    let (vt_in, vt_ix) = vt_in.into_owned().into_raw_vec_and_offset();
    let (vbp_in, vbp_ix) = vbp_in.into_owned().into_raw_vec_and_offset();
    let (vbq_in, vbq_ix) = vbq_in.into_owned().into_raw_vec_and_offset();
    let (vbv_in, vbv_ix) = vbv_in.into_owned().into_raw_vec_and_offset();
    let (vsp_in, vsp_ix) = vsp_in.into_owned().into_raw_vec_and_offset();
    let (vsq_in, vsq_ix) = vsq_in.into_owned().into_raw_vec_and_offset();
    let (vsv_in, vsv_ix) = vsv_in.into_owned().into_raw_vec_and_offset();
    if verbose >= 2  {
      println!("lib.rs->sip_arrow(Vb={}, t_ix{}, b[p,q,v]=[{},{},{}], s[p,q,v]=[{},{},{}])",
        verbose, ps_ix(vt_ix), ps_ix(vbp_ix), ps_ix(vbq_ix), ps_ix(vbv_ix), 
        ps_ix(vsp_ix), ps_ix(vsq_ix), ps_ix(vsv_ix));
    }
    let n_t = vt_in.len();
    if verbose >= 1 {
      println!("lib.rs->sip_arrow(Vb={},vt_in.len={}) --- Joined.", verbose, n_t);
    }
    let n_bp = vbp_in.len();
    let n_bq = vbq_in.len();
    let n_bv = vbv_in.len();
    let n_sp = vsp_in.len();
    let n_sq = vsq_in.len();
    let n_sv = vsv_in.len();

    if (n_bp != n_bq) || (n_bp != n_bv) || (n_t != n_bp) || 
       (n_t != n_sp) || (n_t != n_sq) || (n_t != n_sv) {
      println!("lib.rs->sip_arrow(Vb={}, (T,bP,bQ,bV,sP,sQ,sV)=[{}|{},{},{}|{},{},{}] ",
        verbose, n_t, n_bp, n_bq, n_bv, n_sp, n_sq, n_sv);
      //let out_error = PyArray1::<PublishSip>::zeros_bound(py, (0) as usize, false); 
      //return out_error;
    }
    let out_rec_sip:RecSip = compute_sip_nbbo_rec_sip(n_v as u8, vt_in,
      vbq_in, vsq_in, vbp_in, vsp_in, vbv_in, vsv_in, verbose as u8).unwrap();
    if verbose >= 1 {
      println!("lib.rs->sip_arrow -- complete with out_rec_sip of length {}.",
         out_rec_sip.i_r);
    }
    let col00_time = Arc::new(Int64Array::from(out_rec_sip.vt[0..out_rec_sip.i_r].to_vec())) as _;
    let col01_nbb = Arc::new(Float64Array::from(out_rec_sip.vb_p[0..out_rec_sip.i_r].to_vec())) as _;
    let col02_nbo = Arc::new(Float64Array::from(out_rec_sip.vs_p[0..out_rec_sip.i_r].to_vec())) as _;
    let col03_nbbq = Arc::new(Float64Array::from(out_rec_sip.vb_q[0..out_rec_sip.i_r].to_vec())) as _;
    let col04_nboq = Arc::new(Float64Array::from(out_rec_sip.vs_q[0..out_rec_sip.i_r].to_vec())) as _;
    let col05_b_ii = Arc::new(UInt64Array::from(out_rec_sip.vb_u[0..out_rec_sip.i_r].to_vec())) as _;
    let col06_s_ii = Arc::new(UInt64Array::from(out_rec_sip.vs_u[0..out_rec_sip.i_r].to_vec())) as _;

    let col07_iline = Arc::new(UInt64Array::from(
       to_v64(out_rec_sip.vi_line[0..out_rec_sip.i_r].to_vec() )  
    )  ) as _;
    let col08_b_nwv = Arc::new(UInt8Array::from(out_rec_sip.vb_nwv[0..out_rec_sip.i_r].to_vec())) as _;
    let col09_s_nwv = Arc::new(UInt8Array::from(out_rec_sip.vs_nwv[0..out_rec_sip.i_r].to_vec())) as _;

    if verbose >= 1 {
     println!("lib.rs->sip_arrow -- about to try to create RecordBatch");
    }
    let batch = RecordBatch::try_from_iter([("time", col00_time), 
      ("nbb", col01_nbb),("nbo",col02_nbo),
      ("nbbq", col03_nbbq),("nboq",col04_nboq),
      ("b_II", col05_b_ii),("s_II",col06_s_ii),
      ("iline", col07_iline),
      ("num_b", col08_b_nwv),("num_s",col09_s_nwv)
    ]).unwrap();
    // Yes we tried a buncy of ways to convert the batch to pyarrow
    //   But there appears to be a "gated" off warning when trying
    //   to directly access to_pyarrow() trait.
    //return arrowPyTable::from(batch);
    //return batch.to_pyarrow();
    if verbose >= 1 {
      println!("lib.rs->sip_arrow -- ready to try pyarrow.");
    }
    return recordbatch_to_pyarrow(&batch, py);
  }

  Ok(())
}
