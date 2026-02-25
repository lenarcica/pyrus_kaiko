/// cumulate.rs
///   Cumulation Algorithm
///
///   Alan Lenarcic 2025
///
///   GPL_v2 licensed (you might be able to do similar accumulation with duckdb!)
///
///   Preperatory algorithm to Orderbook totalling implemented in Rust.  This solves the "need to
///   quickly accumulate the order table" problem we face is assembling order plots.
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
//use crate::ord_struct::{TQ};
pub type TQ = f64;
use arrow::array::{Int64Array, Float64Array, BinaryArray,TimestampNanosecondArray};
use arrow::datatypes::{DataType, TimeUnit};

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
     if i0 >= nn-1 { false 
     } else if v_t0[s0[i0+1]] != on_t { false
     } else if v_ven[s0[i0+1]] != on_ven { false
     } else if v_price[s0[i0+1]] != on_p { false
     } else if v_b0s1[s0[i0+1]] != on_s { false
     } else { true }
  }}
  macro_rules!next_t1_is_same {
     () => {
     if i1 >= nn-1 { false 
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
