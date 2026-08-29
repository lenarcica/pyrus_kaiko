/// cumulate.rs
///   Cumulation Algorithm
///   Preperatory to Orderbook totalling
///
/// Typically if order messages arrive in tabular form
///
///
/// order-id | seq |side | venue |time_origin  | rest_price | rest_qty | time_conclusion  | 
/// --------------------------------------------------------------------------------
///    90001 | 001 |   B | NY    |10:00:00.000 |    $100.00  |      100 |   10:00:03.433  |
///    90001 | 002 |   B | NY    |10:00:03.433 |    $100.00  |       50 |   10:00:05.343  |
///    90002 | 003 |   B | BOS   |10:00:04.323 |     $99.99  |      200 |   10:00:06.323  |
///    90001 | 004 |   B | NY    |10:00:05.343 |    $100.00  |        0 |   xx:xx:xx:xxx  |
///    90003 | 005 |   S | NJ    |10:00:05.343 |    $100.01  |       50 |   10:00:05.353  |
///    ....
///
/// However, for orderbook totals we want to transform this data into a tabulated
///  form where the size of resting quantity at every side/price/venue combination in order 
///  of time.  That is we want to know an aggregated time series
///
///  side | rest_price | venue | time         |     rqty |   mqty|
///  -------------------------------------------------------------
///     B |       .001 |    NY | 10:00:00.323 |        20|     20|
///     B |       .001 |    NY | 11:40:32.432 |         0|    200|
///     B |       .001 |   BOS | 10:32:03.43  |       200|    220|
///     B |       .001 |   BOS | 12:30:03.43  |         0|      0|
///     B |       1.01 |    NY | 09:59:32.343 |      1000|   1000|
///     B |       1.01 |    NY | 10:00:04.632 |      1400|   1400|
///     B |       1.01 |    NY | 10:00:06.343 |      2343|   2343|
///     B |       1.01 |    NY | 10:00:32:432 |         0|      0|
///     .....
///
/// This is the format our orderbook algorithm will be able to tackle
///   The challenge of generating a cumulation algorithm however is that involves
///   careful process to tally the data.
/// Note that if the data is originally supplied by order-id/time, the "time_conclusion"
///   is just a shift of the time problem.
///

//use std::mem;
use crate::Arc;
//use crate::ord_struct::{ TBS };
use crate::fill_unfilled_seq::{ TTime };
//use crate::fill_unfilled_seq::{TPrice, TTimePrice,  TSideS, TTmPrSd};
use arrow::record_batch::RecordBatch;
//use cheats::RecordBatchToPyArrow; // If Traits don't help me I cheat
// This will go to lib.rs probably.
//use crate::cheats::recordbatch_to_pyarrow;
//use crate::ord_struct::{TQ, TP, TBS, TPi, TNRi};
pub type TQ = f64;
pub type TP = f64;
pub type TNRi = u8; pub type TPi = usize;
use std::fmt;
#[derive(PartialEq, Eq, Copy, Clone)]
pub enum TBS {
  B,S
}
impl fmt::Display for TBS {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    match self {
      TBS::B => write!(f, "B"),
      TBS::S => write!(f, "S"),
    }
  }
}
use arrow::array::{Int64Array, Float64Array, BinaryArray,TimestampNanosecondArray};
//use arrow::datatypes::{DataType, TimeUnit};

// Could not get "zerocopy" to import, so we hacked instead on "B"/"S" items
// necessary to convert vec<char> to &[u8]?
//use zerocopy::AsBytes;

//
//use arrow::array::{UInt64Array, UInt8Array, Array, BinaryArray};
//let mut x = 5;
//let mut y = 42;

//mem::swap(&mut x, &mut y);

// Note triomax01 doesn't work because v_p/v_s/v_t aren't in "macro scope"
//  This appears quite different from C macro in terms of writing code where desired.
/***
macro_rules!triomax01 {
  ($i0:expr, $i1:expr) => {
    (if v_s[$i0] < v_s[$i1] { 
    } else if v_s[$i0] > v_s[$i1] { 0
    } else if v_p[$i0] < v_p[$i1] { 1
    } else if v_p[$i0] > v_p[$i1] { 0
    } else if v_p[$i0] < v_p[$i1] { 1
    } else if v_v[$i0] > v_v[$i1] { 0
    } else if v_v[$i0] < v_v[$i1] { 1
    } else if v_t[$i0] < v_t[$i1] { 1 
    } else {0})
  }
}
***/

pub fn idx_sort_trio(v_s: &[i8], v_p: &[f64], v_v: &[i64], v_t:&[i64]) -> Vec<usize> {
  // It appears that to make use of variables namespace for this loop, we need macro inside function
  // I suppose this is easier than passing 5 expressions to triomax02 and remembering
  // which entries work
  macro_rules!triomax02 {
    ($i0:expr, $i1:expr) => {

    (      if v_s[ $i0 ] < v_s[ $i1 ] { 1 
    } else if v_s[ $i0 ] > v_s[ $i1 ] { 0
    } else if v_p[ $i0 ] < v_p[ $i1 ] { 1
    } else if v_p[ $i0 ] > v_p[ $i1 ] { 0
    } else if v_p[ $i0 ] < v_p[ $i1 ] { 1
    } else if v_v[ $i0 ] > v_v[ $i1 ] { 0
    } else if v_v[ $i0 ] < v_v[ $i1 ] { 1
    } else if v_t[ $i0 ] < v_t[ $i1 ] { 1 
    } else {0})
  
    }
  }
  let n: usize = v_s.len() as usize;
  let mut i: usize; 
  let mut j: usize; let mut l: usize = (n >> 1) + 1;
  let mut rra_ix:usize;
  let mut ir: usize = n-1;
  let mut out_ix:Vec<usize> = vec![0 as usize;n];
  for ii in 0..n { out_ix[ii as usize] = ii as usize; } 

  if n < 2 { return out_ix; }
  let mut nloop = 0;
  loop {
    if l > 0 {
      l = l-1;
      rra_ix = out_ix[l];
    } else {
      rra_ix = out_ix[ir];
      out_ix[ir] = out_ix[0];
      ir = ir -1;
      if ir == 0 {
        out_ix[0] = rra_ix;
        break;
      }
    }
    i = l;
    j = (l+1) + (l+1)-1;
    while j  <= ir {
      if (j < ir) && (1==(triomax02!( out_ix[j], out_ix[j+1] ))) {
        j=j+1;
      }
      if 1 == (triomax02!( rra_ix, out_ix[j]   )) {
        out_ix[i] = out_ix[j];
        i=j;
        j = ((j+1) << 1)-1;
      } else {
        break;
      }
      nloop = nloop + 1;
      if nloop > n*n {
        println!("-- idx_sort_trio: some error nloop={}, n={} internally, j={}, ir={}, l={}",
          nloop, n, j, ir, l);
        break;
      }
    }
    out_ix[i] = rra_ix;
    nloop = nloop + 1;
    if nloop > n*n {
      println!("-- Forced to break from Loop on nLoop = {}, n={}, possibly algo fail",
        nloop, n);
      break;
    }
  }
  return out_ix;
}


pub fn idx_sort_ven_after_t(v_s: &[i8], v_p: &[f64], v_t:&[i64], v_ven:&[i64], n: usize) -> Vec<usize> {
  // It appears that to make use of variables namespace for this loop, we need macro inside function
  // I suppose this is easier than passing 5 expressions to triomax02 and remembering
  // which entries work
  macro_rules!duomax02 {
    ($i0:expr, $i1:expr) => {

    (      if v_s[ $i0 ] < v_s[ $i1 ] { 1 
    } else if v_s[ $i0 ] > v_s[ $i1 ] { 0
    } else if v_p[ $i0 ] < v_p[ $i1 ] { 1
    } else if v_p[ $i0 ] > v_p[ $i1 ] { 0
    } else if v_t[ $i0 ] < v_t[ $i1 ] { 1 
    } else if v_t[ $i0 ] > v_t[ $i1 ] { 0
    } else if v_ven[ $i0 ] < v_ven[ $i1 ] { 1
    } else {0})
  
    }
  }
  if n > v_s.len() {
    println!("Error idx_sort_duo, can't sort because n given = {} but v_s.len = {}", n, v_s.len());
  }
  let mut i: usize; let mut ir: usize;
  let mut j: usize; let mut l: usize = (n >> 1) + 1;
  let mut rra_ix;
  ir = n-1;
  let mut out_ix:Vec<usize> = vec![0 as usize;n];
  for ii in 0..n { out_ix[ii as usize] = ii as usize; } 

  if n < 2 { return out_ix; }
  let mut nloop = 0;
  loop {
    if l > 0 {
      l = l-1;
      rra_ix = out_ix[l];
    } else {
      rra_ix = out_ix[ir];
      out_ix[ir] = out_ix[0];
      ir = ir -1;
      if ir == 0 {
        out_ix[0] = rra_ix;
        break;
      }
    }
    i = l;
    j = (l+1) + (l+1)-1;
    while j  <= ir {
      if (j < ir) && (1==(duomax02!( out_ix[j] ,  out_ix[j+1]  ))) {
        j=j+1;
      }
      if 1 == (duomax02!( rra_ix, out_ix[j]  )) {
        out_ix[i] = out_ix[j];
        i=j;
        j = ((j+1) << 1)-1;
      } else {
        break;
      }
      nloop = nloop + 1;
      if nloop > n*n {
        println!("-- idx_sort_trio: some error nloop={}, n={} internally, j={}, ir={}, l={}",
          nloop, n, j, ir, l);
        break;
      }
    }
    out_ix[i] = rra_ix;
    nloop = nloop + 1;
    if nloop > n*n {
      println!("-- Forced to break from Loop on nloop = {}, n={}, possibly algo fail",
        nloop, n);
      break;
    }
  }
  return out_ix;
}

/**************************
 *  Algo above adapted From C recipies (had to solve for 0-index)
 * https://www.foo.be/docs-free/Numerical_Recipe_In_C/c8-3.pdf
 * Not actually that different from other implementations of heapsort
{
  unsigned long i,ir,j,l;
  float rra;
  if (n < 2) return;
  l=(n >> 1)+1; // l initialized at midpoint of n
  ir=n;
  for (;;) {
    if (l > 1) { // Still in hiring phase.
      rra=ra[--l];
    } else { // In retirement-and-promotion phase.
      rra=ra[ir]; // Clear a space at end of array.
      ra[ir]=ra[1]; // Retire the top of the heap into it.
      if (--ir == 1) { // Done with the last promotion.
        ra[1]=rra; //The least competent worker of all!
        break;
      }
    }
    i=l;  // Whether in the hiring phase or 
          // promotion phase, we here set up to sift down 
          // element rra to its proper level.
    j=l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) {
        j++; //Compare to the better underling
      }
      if (rra < ra[j]) {                  //Demote rra.
        ra[i]=ra[j];
        i=j;
        j <<= 1;
      } else {
        break;   //Found rra’s level. Terminate the sift-down.
      }
    }
    ra[i]=rra; Put rra into its slot.
  }
}
********************************************************/

/***********************************************************************************/

///  cumulation_algo()  
///   Cumulation Algorithm
///   Preperatory to Orderbook totalling
///
/// Typically if order messages arrive in tabular form
///
///
/// order-id | seq |side | venue |time_origin  | rest_price | rest_qty | time_conclusion  | 
/// --------------------------------------------------------------------------------
///    90001 | 001 |   B | NY    |10:00:00.000 |    $100.00  |      100 |   10:00:03.433  |
///    90001 | 002 |   B | NY    |10:00:03.433 |    $100.00  |       50 |   10:00:05.343  |
///    90002 | 003 |   B | BOS   |10:00:04.323 |     $99.99  |      200 |   10:00:06.323  |
///    90001 | 004 |   B | NY    |10:00:05.343 |    $100.00  |        0 |   xx:xx:xx:xxx  |
///    90003 | 005 |   S | NJ    |10:00:05.343 |    $100.01  |       50 |   10:00:05.353  |
///    ....
///
/// However, for orderbook totals we want to transform this data into a tabulated
///  form where the size of resting quantity at every side/price/venue combination in order 
///  of time.  That is we want to know an aggregated time series
///
///  side | rest_price | venue | time         |   cum_qty 
///  -----------------------------------------------------
///     B |       .001 |    NY | 10:00:00.323 |        20
///     B |       .001 |    NY | 11:40:32.432 |         0
///     B |       .001 |   BOS | 10:32:03.43  |       200
///     B |       .001 |   BOS | 12:30:03.43  |         0
///     B |       1.01 |    NY | 09:59:32.343 |      1000
///     B |       1.01 |    NY | 10:00:04.632 |      1400
///     B |       1.01 |    NY | 10:00:06.343 |      2343
///     B |       1.01 |    NY | 10:00:32:432 |         0
///     .....
///
/// This is the format our orderbook algorithm will be able to tackle
///   The challenge of generating a cumulation algorithm however is that involves
///   careful process to tally the data.
/// Note that if the data is originally supplied by order-id/time, the "time_conclusion"
///   is just a shift of the time problem.
///
///    This algorithm is designed to accumulate quantity amongst multiple resting
///    orders at same side/price/aggregation value over the timestamps in the day.
///    v_side/v_price/v_ven are the inputs that represent disting groups
///    v_t0 and v_t1 are open and close times for each order message state.
///
///    Depending on how data is sorted, "v_t1" could be implied by later records or
///    are easier to generate for with ".shift()" calculations in polars/pandas. 
pub fn cumulation_algo(v_b0s1: &[i8], v_price: &[f64], v_ven: &[i64],
                        v_qty: &[TQ], v_t0:&[TTime], v_t1:&[TTime], verbose: u8) -> RecordBatch {
  let nn = v_b0s1.len();

  // For some reason Rust hates it when we initiate mutable elements with values but
  //  don't use initial case.
  let mut on_t: TTime = 0 as i64; let mut on_s: i8 = -1;
  let mut on_p: f64 = 0.0 as f64;  let mut on_ven: i64 = -1 as i64;
  let mut prev_q: TQ;
  //let mut on_p: f64 = 0.0 as f64; let mut on_ven: i64 = 0 as i64;
  let mut c_q: TQ; //let mut prev_q: TQ = 0.0 as TQ;
  
  //let mut c_mq: TQ = 0.0 as TQ;
  // Not sure why we need to test the trio sort.
  //let mut indx0 = (0..v_b0s1.len()).collect::<Vec<_>>();
  //indx0.sort_by(|ai, bi| v_b0s1[ai].cmp(v_b0s1[bi]) );
  if verbose >= 1 {
    println!("cumulation_algo(V={},len={}) -- Start with sorts. ", verbose, nn);
  }
  let s0:Vec<usize> = idx_sort_trio(&v_b0s1, &v_price, &v_ven, &v_t0);
  let s1:Vec<usize> = idx_sort_trio(&v_b0s1, &v_price, &v_ven, &v_t1);


  let mut i0:usize = 0; let mut i1:usize = 0;
  //let mut ii0:usize = 0; let mut ii1: usize = 0;
  let mut tt_0: TTime = v_t0[s0[i0]]; let mut ss_0: i8 = v_b0s1[s0[i0]];
  let mut pp_0: f64 = v_price[s0[i0]]; let mut vv_0: i64 = v_ven[s0[i0]];

  let mut tt_1: TTime = v_t1[s1[i1]]; let mut ss_1: i8 = v_b0s1[s1[i1]];
  let mut pp_1: f64 = v_price[s1[i1]]; let mut vv_1: i64 = v_ven[s1[i1]];

  let mut new_level = 0;
  //let mut out_s:Vec<char> = vec!['u' as char; nn*2];
  let mut out_s:Vec<i8> = vec![-1 as i8;nn*2];
  let mut out_p:Vec<f64> = vec![0.0 as f64; nn*2];
  let mut out_v:Vec<i64> = vec![0 as i64; nn*2];
  let mut out_t:Vec<TTime> = vec![0 as TTime; nn*2];
  let mut out_rq:Vec<TQ> = vec![0.0 as TQ; nn*2];
  //let mut out_mq:Vec<TQ> = vec![0.0 as TQ; nn*2];
  let mut ii_print:usize = 0;
  let mut t_loops = 0;
  /***************************************************************************/
  // Helper Macros need to be in function scope
  //  These macros simplify typing and help the algorithm to look smoother
  //  For simplicity's sake, we are putting the functions in scope so they see
  //  other data.
  macro_rules!t0_or_t1 {
    () => {
      if i0 >= nn { 1
      } else if i1 >= nn { 1
      } else if (v_b0s1[s0[i0]] < v_b0s1[s1[i1]]) { 0 
      } else if v_b0s1[s0[i0]] > v_b0s1[s1[i1]] { 1
      } else if v_price[s0[i0]] < v_price[s1[i1]] { 0 
      } else if v_price[s0[i0]] > v_price[s1[i1]] { 1
      } else if v_ven[s0[i0]] < v_ven[s1[i1]] { 0
      } else if v_ven[s0[i0]] > v_ven[s1[i1]] { 1 
      } else if v_t0[s0[i0]] < v_t1[s1[i1]] { 0 
      } else if v_t0[s0[i0]] > v_t1[s1[i1]] { 1 } else {2}
    }
  }
  macro_rules!next_t0_is_same {
     () => {
     if i0 >= nn-2 { false 
     } else if v_t0[s0[i0+1]] != on_t { false
     } else if v_ven[s0[i0+1]] != on_ven { false
     } else if v_price[s0[i0+1]] != on_p { false
     } else if v_b0s1[s0[i0+1]] != on_s { false
     } else { true }
  }}
  macro_rules!next_t1_is_same {
     () => {
     if i1 >= nn-2 { false 
     } else if v_t1[s1[i1+1]] != on_t {false 
     } else if v_ven[s1[i1+1]] != on_ven {false 
     } else if v_price[s1[i1+1]] != on_p {false 
     } else if v_b0s1[s1[i1+1]] != on_s { false 
     } else { true }
  }}
  let mut t0tot1;
  macro_rules!bool_print {
    () => {
      if (new_level > 0) && (c_q == 0.0 as TQ) { false 
      } else if  t0tot1 > 1  { false 
      } else if ( t0tot1 == 0 ) && ( next_t0_is_same!() ) 
             && ( i0 < nn-1 ) && (v_t0[s0[i0+1]] == on_t) {  false
      } else if ( t0tot1 == 1 ) && ( next_t1_is_same!() ) 
              && ( i1 < nn-1 ) && (v_t1[s1[i1+1]] == on_t) { false 
      } else if  c_q != prev_q  { true } else { false }
  }}
  
  macro_rules!v0_same_group {
     () => { if (ss_0 == on_s) && (pp_0 == on_p) && (vv_0 == on_ven) { true } else { false } }
  }
  macro_rules!v1_same_group { 
     () => { if (ss_1 == on_s) && (pp_1 == on_p) && (vv_1 == on_ven) { true } else { false } }
  }
  macro_rules!print_to_arr {
     () => {
       out_s[ii_print] = on_s; out_v[ii_print] = on_ven;
       out_rq[ii_print] = c_q; out_p[ii_print] = on_p;
       out_t[ii_print] = on_t;
       ii_print += 1
  }}

  // Strategy, first accumulate the Side/Price/Vendor/Time methodology.
  // Then we accumulate all Vendors together with a second, easier sort and accumulate.
   // Actual Algorithm for cumulation:
  prev_q = -1.0; c_q = 0.0;
  while (i0 < nn) || (i1 < nn) {
     t0tot1 = t0_or_t1!();
     if verbose >= 3 {
       println!("--- Cumulate (i0={}/{}, i1={}/{}, on_s={}, on_p={}, on_ven={}, on_t={},: t0_to_t1={}",
         i0, nn, i1, nn, on_s, on_p, on_ven, on_t, t0tot1);
     }
     //prev_q = c_q;  // this value may only be valuable in debugging previous change.
     if verbose >= 1 {
       if (on_s != ss_0) && (on_s != ss_1) {
         println!("--- Cumulate (i0={}/{}, i1={}/{}, on_s={}, on_p={}, on_ven={},on_t={} --  ss_0={}, ss_1={}: expect a change of sign.",
           i0, nn, i1, nn, on_s, on_p, on_ven, on_t, ss_0, ss_1);
       }
     }
     if t0tot1 != 1 {
       c_q += v_qty[s0[i0 as usize] as usize];  on_t = tt_0; on_s = ss_0; on_ven = vv_0; on_p = pp_0;
     } else {
       c_q -= v_qty[s1[i1 as usize] as usize];  on_t = tt_1; on_s = ss_1; on_ven = vv_1; on_p = pp_1;
     }
     if bool_print!() {
       if verbose >= 2 {
         println!("    ---- Print ii_print={}, c_q = {}, on_s={},on_p={}, on_ven={}, on_t={}",
           ii_print, c_q, on_s, on_p, on_ven, on_t);
       }
       print_to_arr!();
       prev_q = c_q;
     } else {
        if verbose >= 4 {
           if (new_level > 0)  && (c_q == 0.0 as TQ) { println!(" --- No Print i0={}, i1={}, new_level={} c_q={}",
             i0, i1, new_level, c_q); 
           } else if t0tot1 > 1 { println!(" --- NoPrint, t0t1={}, i0={}, i1={}", t0tot1, i0, i1);
           } else if (t0tot1 == 0) && (next_t0_is_same!()) { println!(" --- No Print, t0t1=0, next_t0_is_same={}", next_t0_is_same!());
           } else if (t0tot1 == 1) && (next_t1_is_same!()) { println!(" --- No Print, t0t1=1, next_t1_is_same={}", next_t1_is_same!());
           } else { println!(" --- No Print because c_q={}, prev_q={}, i0={}/{}, i1={}/{}", c_q, prev_q, i0, nn, i1, nn);
           }
        }
     }
     if t0tot1 != 1 {
       i0 += 1;
       if i0 < nn {
         tt_0 = v_t0[s0[i0 as usize] as usize]; vv_0 = v_ven[s0[i0 as usize] as usize]; 
         ss_0 = v_b0s1[s0[i0 as usize] as usize]; pp_0 = v_price[s0[i0 as usize] as usize];
       }
     } else {
       i1 += 1;
       if i1 < nn {
         tt_1 = v_t1[s1[i1 as usize] as usize]; vv_1 = v_ven[s1[i1 as usize] as usize]; 
         ss_1 = v_b0s1[s1[i1 as usize] as usize]; pp_1 = v_price[s1[i1 as usize] as usize];
       }
     }
  
    if ((i0 < nn) && (v0_same_group!())) || ((i1 < nn) && (v1_same_group!())) {
      new_level = 0;
    } else {
      new_level = 1; prev_q = -1000.0; c_q = 0.0;
      if c_q > 0.0001 {
        println!("--- Issue at t_loops={}/{}, i0={},i1={}-on_t={}, new_level but c_q={} was not 0, prev_q was {}, on_s={}.",
          t_loops, nn, i0, i1, on_t, c_q, prev_q, if on_s == 0 {'b'} else {'s'});
      }
    }
    t_loops += 1;
    if t_loops >= nn*2 {
      break;
    }
  }

  if verbose >= 1{
    println!(" cumulate.rs -- cumulation finishes first step with iiprint={}/{}",  ii_print, nn);
  }
  // Market Second level accumulate should be "easier"
  // First we sort by just side/price/time (ignore Venue key now).
  // It is straightforward to accumulate right now and fill in the last column.
  //
  // The hard part is we need to know if we are still using the same price/side we used before
  // as well we need to know if the current index is a "new quantity print" starting from 0 or
  // whether it is accumulated from the value above it.
  let ss_t0:Vec<usize> = idx_sort_ven_after_t(&out_s, &out_p, &out_t, &out_v, ii_print);
  let mut out_mq = vec![0.0 as TQ; ii_print];
  let mut on_mq = out_rq[ss_t0[0]]; 
  // on_p = out_p[ss_t0[0]]; on_s = out_s[ss_t0[0]];
  out_mq[ss_t0[0]] = on_mq;
  for ii in 1..ii_print {
    if (out_s[ss_t0[ii-1]] == out_s[ss_t0[ii]]) && (out_p[ss_t0[ii-1]] == out_p[ss_t0[ii]]) {
       if ss_t0[ii] == 0 {
          on_mq += out_rq[ss_t0[ii]];
       } else if (out_s[ss_t0[ii]] == out_s[ss_t0[ii]-1]) && (out_p[ss_t0[ii]] == out_p[ss_t0[ii]-1]) &&
                 (out_v[ss_t0[ii]] == out_v[ss_t0[ii]-1]) {
         on_mq = on_mq + out_rq[ss_t0[ii]] - out_rq[ss_t0[ii]-1];
       } else {
         on_mq = on_mq + out_rq[ss_t0[ii]];
       }
    } else {
      on_mq = out_rq[ss_t0[ii]];
    }
    out_mq[ss_t0[ii]] = on_mq;
  }
   
  // Permute Side Price Time, no venue:
  //   Unfortunately permutators would require external library.  We already have ss_t0 so we don't
  //   need to sort. 
  // Note: Laziest way to convert Vec<char> to Vec<&[u8]> which BinaryArray actually supports.
  let mut side_values: Vec<&[u8]> =  vec![b"u"; ii_print];
  for ii in 0..ii_print { 
    //side_values[ii] = if (out_s[ss_t0[ii]] == 'B') || (out_s[ss_t0[ii]] == 'b') { b"b" } else if (out_s[ss_t0[ii]] == 's') || (out_s[ss_t0[ii]] == 'S') { b"s" } else { b"u"}; 
    side_values[ii] = if out_s[ss_t0[ii]] == 0 { b"b" } else if out_s[ss_t0[ii]] == 1 { b"s" } else { b"u"}; 
  }
  let mut v_p = vec![0.0 as f64; ii_print];
  for ii in 0..ii_print { v_p[ii] = out_p[ss_t0[ii]]; }
  let mut v_t = vec![0.0 as i64; ii_print];
  for ii in 0..ii_print { v_t[ii] = out_t[ss_t0[ii]]; }
  let mut v_ven = vec![0.0 as i64; ii_print];
  for ii in 0..ii_print { v_ven[ii] = out_v[ss_t0[ii]]; }
  let mut v_rq = vec![0.0 as f64; ii_print];
  for ii in 0..ii_print { v_rq[ii] = out_rq[ss_t0[ii]]; }
  let mut v_mq = vec![0.0 as f64; ii_print];
  for ii in 0..ii_print { v_mq[ii] = out_mq[ss_t0[ii]]; }

  let col00_side = Arc::new(BinaryArray::from_vec(side_values)) as _;
  let col01_price = Arc::new(Float64Array::from(v_p)) as _;
  let col02_time = Arc::new(TimestampNanosecondArray::from(v_t)) as _;  // Output now Nanosecond type
  //let col02_time = Arc::new(Int64Array::from(v_t)) as _;
  let col03_agg_id = Arc::new(Int64Array::from(v_ven)) as _;
  let col04_crqty = Arc::new(Float64Array::from(v_rq)) as _;
  let col05_cmqty = Arc::new(Float64Array::from(v_mq)) as _;
  // Previously we didn't want to permute but now we do.
  //let col00_side = Arc::new(BinaryArray::from_vec(side_values)) as _;
  //let col01_price = Arc::new(Float64Array::from(out_p[0..ii_print].to_vec().sort_by)) as _;
  //let col02_time = Arc::new(Int64Array::from(out_t[0..ii_print].to_vec())) as _;
  //let col03_agg_id = Arc::new(Int64Array::from(out_v[0..ii_print].to_vec())) as _;
  //let col04_crqty = Arc::new(Float64Array::from(out_rq[0..ii_print].to_vec())) as _;
  //let col05_cmqty = Arc::new(Float64Array::from(out_mq)) as _;


  if verbose >= 1 {
    println!("cumulator.rs -- about to try to create RecordBatch");
  }
  let batch = RecordBatch::try_from_iter([("side", col00_side), 
      ("price", col01_price),("time", col02_time),
      ("related_id",col03_agg_id),
      ("cr_qty",col04_crqty), ("cm_qty", col05_cmqty)
   ]).unwrap();
  return batch;
} 
pub struct CData {
  pub vt: Vec<TTime>,
  pub vpi: Vec<TPi>,
  pub vq: Vec<TQ>,
  pub nn: usize
}
impl CData {
  pub fn sort_tp(self: &mut Self) {
    let nn = self.nn;
    let mut idx = vec![0 as usize;nn];
    for ii in 0..nn { idx[ii] = ii; }
    idx.sort_by(|ai,bi| if self.vt[*ai]==self.vt[*bi] {self.vpi[*ai].cmp(&self.vpi[*bi])} else { self.vt[*ai].cmp(&self.vt[*bi]) } );
    self.vpi = idx.iter().map(|&i| self.vpi[i]).collect();
    self.vq = idx.iter().map(|&i| self.vq[i]).collect();
    self.vt = idx.iter().map(|&i| self.vt[i]).collect();
  }
}
pub struct InputSideStruct {
  pub bs: TBS,
  pub u_p: Vec<TP>,
  pub vpi: Vec<TPi>,
  pub vt:  Vec<TTime>,
  pub vq:  Vec<TQ>,
  pub vr:  Vec<TNRi>
}
impl InputSideStruct {
  pub fn new(in_bs: TBS, in_u_p: Vec<TP>, in_vpi: Vec<TPi>, in_vt: Vec<TTime>, in_vq: Vec<TQ>) -> Self {
    return InputSideStruct{bs:in_bs, u_p:in_u_p, vpi:in_vpi, vt: in_vt, vq: in_vq, vr:vec![0 as TNRi;0]};
  }
}
pub struct InputStruct {
  pub b: InputSideStruct, pub s: InputSideStruct
}
impl InputStruct {
  pub fn new(oss_b:InputSideStruct, oss_s:InputSideStruct) -> Self{
    return InputStruct{b:oss_b,s:oss_s};
  }
}
pub fn idx_spt(vpi: &Vec<TPi>,vt: &Vec<TTime>) -> Option<Vec<usize>> {
  if vpi.len() != vt.len() { return None; }
  let nn: usize = vpi.len(); 
  let mut idx = vec![0 as usize;nn];
  for ii in 0..nn {  idx[ii] = ii; }
  //|ai, bi| v_b0s1[ai].cmp(v_b0s1[bi])  
  idx.sort_by(|ai,bi| if vpi[*ai]==vpi[*bi] { vt[*ai].cmp(&vt[*bi]) } else { vpi[*ai].cmp(&vpi[*bi]) } );
  return Some(idx);
}
pub fn make_c_data(idx: Vec<usize>, vpi: &Vec<TPi>, vt: &Vec<TTime>, vq: &Vec<TQ>) -> Option<CData> {
  let nn:usize = vpi.len();
  if  (nn != vt.len()) || (nn != vq.len()) {
    //Err(String::from(format!["make_c_data: Hey Error (idx={}, vpi={}, vt={}, vq={})",idx.len(),vpi.len(),vt.len(), vq.len()]));
    return None;
  }
  let mut o_vpi:Vec<TPi> = vec![0 as TPi;nn];
  let mut o_vt:Vec<TTime> = vec![0 as TTime;nn];
  let mut o_vq:Vec<TQ> = vec![0 as TQ;nn];
  for ii in 0..nn {
    o_vpi[ii] = vpi[idx[ii]];
    o_vt[ii] = vt[idx[ii]];
    o_vq[ii] = vq[idx[ii]];
  }
  return Some(CData{vt:o_vt,vpi:o_vpi,vq:o_vq,nn:nn});
} 
pub fn cumulate_dside(bs: TBS, u_p: Vec<TP>, vo: Vec<TTime>,vc: Vec<TTime>,vq: Vec<TQ>,vpi: Vec<TPi>) -> Result<InputSideStruct,String> {
  let nn = vo.len();
  if (nn != vo.len()) || (nn != vc.len()) || (nn != vq.len()) || (nn != vpi.len()) {
    return Err(String::from(format!["cumulate_dside, length (o={},c={},q={},pi={})=", vo.len(), vc.len(), vq.len(), vpi.len()]));
  }
  if nn == 0 {
    //Size Zero table, I guess thats best we can do.  More useful than a None
    return Ok(InputSideStruct{bs:bs,u_p:u_p.clone(), vpi:vec![0 as TPi;0],vt:vec![0 as TTime;0],vq:vec![0 as TQ;0],vr:vec![0 as TNRi;0]});
  }
  let idx_o:Vec<usize> = idx_spt(&vpi,&vo).expect("cumulate_dside, idx_spt(vo,t) should be even size, honestly already tested.");
  let idx_c:Vec<usize> = idx_spt(&vpi,&vc).expect("cumulate_dside, idx_spt(vc,t) should be even size, honestly already tested.");
  let c_data_c:CData = make_c_data(idx_o, &vpi, &vo, &vq).expect("cumulate_dside: make_c_data(idx_o,v_o,vt,vq) should have worked");
  let c_data_o:CData = make_c_data(idx_c, &vpi, &vc, &vq).expect("cumulate_dside: make_c_data(idx_c,v_c,vt,vq) should have worked");

  let mut ip:usize=0; let n2:usize = nn.checked_add(nn).expect("cumulate_dsize, nn so large that 2n larger than usize");
  let mut i_o:usize=0; let mut i_c:usize=0; let mut on_q: TQ;
  let mut o_vt:Vec<TTime> = vec![0 as TTime;n2]; let mut o_vq:Vec<TQ> = vec![0 as TQ;n2]; let mut o_vpi:Vec<TPi> = vec![0 as TPi;n2];
  o_vt[ip] = c_data_o.vt[i_o]; o_vq[ip] = c_data_o.vq[i_o]; o_vpi[ip] = c_data_o.vpi[i_o]; on_q = o_vq[0]; ip = 1;
  let mut on_pi:usize = o_vpi[0];
  for ii in 1..n2 {
    if (i_o >=nn) || (c_data_c.vpi[i_c] < c_data_o.vpi[i_o]) || (c_data_c.vt[i_c] < c_data_o.vt[i_o]) {
       if on_q - c_data_c.vq[i_c] < (0 as TQ) {
          return Err(String::from(format!["cumulate_dside, we have an on_q error on ii={}/{}, on_q={}, c_data_c.vq[{}]={}", ii,n2,on_q,i_c,c_data_c.vq[i_c]]));
       }
       o_vt[ip] = c_data_c.vt[i_c]; on_q -= c_data_c.vq[i_c];  i_c = i_c+1;
    } else {
       o_vt[ip] = c_data_o.vt[i_o];
       if c_data_o.vpi[i_o] == on_pi {
         on_q += c_data_o.vq[i_o];
       } else {
         on_q = c_data_o.vq[i_o]; on_pi = c_data_o.vpi[i_o];
       }
       i_o = i_o + 1;
    }
    let ipm1 = ip.saturating_sub(1);
    o_vq[ip] = on_q;  o_vpi[ip] = on_pi;
    if o_vt[ip] != o_vt[ipm1] {
      if o_vpi[ip] == o_vpi[ipm1] && o_vq[ip] == o_vq[ipm1] {
      } else {
        ip += 1;
      }
    } else {
      if o_vpi[ip] != o_vpi[ipm1] {
        ip += 1;
      }
    }
  }
  o_vt.resize(ip,0 as TTime); o_vpi.resize(ip,0 as TPi); o_vq.resize(ip,0 as TQ);
  Ok(InputSideStruct{bs:bs,u_p:u_p.clone(),vpi:o_vpi,vt:o_vt,vq:o_vq,vr:vec![0 as TNRi;0]})
}

pub fn concatenate_cds(bs:TBS, l_cds:&[InputSideStruct], u_p:Vec<TP>) -> Result<InputSideStruct,String>{
  if l_cds[0].vt.len() == 0 {
    println!("concatenaate_cds: the first Market InputSideStruct must be non zero in length!");
    return Err(String::from("concatante_cds: Error, first element of l_cds is zero length."));
  }
  if l_cds.len() == 1 { 
    return Ok(InputSideStruct{bs:l_cds[0].bs,u_p:u_p.clone(),vpi:l_cds[0].vpi.clone(),vt:l_cds[0].vt.clone(),vq:l_cds[0].vq.clone(),
                         vr:vec![0 as TNRi;l_cds[0].vpi.len()]});
  }
  if l_cds.len() > 254 {
    return Err(String::from("HEY, concatenate_cds, we do not permit more than u8 simulataneous Related parties!"));
  }
  let nl:usize = l_cds.len();
  let mut idxL: Vec<usize> = vec![0;nl];
  let mut nL: Vec<usize> = vec![0,nl]; let mut nn = 0;
  let mut mxt: TTime = l_cds[0].vt[l_cds[0].vt.len().saturating_sub(1)];
  for ii in 0..nl {
    nL[ii] = l_cds[ii].vpi.len(); nn += nL[ii];
    if nL[ii] > 0 { if mxt < l_cds[ii].vt[ nL[ii].saturating_sub(1)] { mxt = l_cds[ii].vt[nL[ii].saturating_sub(1)] } }
  }
  let mxt = mxt;  let mxtB = mxt.checked_add(1).expect("concatenate_cds: we wanted to add 1 to mxt");
  let nn = nn;
  let mut o_vpi = vec![0 as TPi;nn]; let mut o_vq = vec![0 as TQ; nn]; let mut o_vt = vec![0 as TTime;nn];
  let mut o_vr = vec![0 as TNRi;nn];
  let nr = l_cds.len() -1; let nru8:u8 = nr.try_into().expect("How could nr not become u8!");
  for ii in 0..nn {
    let mut on_r = 0; 
    let mut on_t = if nL[0] >= idxL[0] { l_cds[0].vt[idxL[0]] } else { mxtB };
    for ir in 1..nr { 
      if (idxL[ir] < nL[ir]) && (l_cds[ir].vt[idxL[ir]] < on_t) {  on_r = ir;  on_t = l_cds[ir].vt[idxL[ir]]; }
    }
    o_vt[ii] = on_t; 
    if idxL[on_r] >= nL[on_r] { 
      return Err(String::from(format!["concatenate_cds: idxL[on_r={}] = {}  but length is {}.", on_r, idxL[on_r], nL[on_r]]));
    }
    o_vq[ii] = l_cds[on_r].vq[idxL[on_r]]; 
    o_vpi[ii] = l_cds[on_r].vpi[idxL[on_r]];
    o_vr[ii] = if on_r == 0 { nru8 } else  { on_r.saturating_sub(1).try_into().expect("How could on_r not subtract to u8!") };
    idxL[on_r]+=1;
  }
  Ok(InputSideStruct{bs:bs,u_p:u_p.clone(),vpi:o_vpi,vt:o_vt,vq:o_vq,vr:o_vr}) 
}
