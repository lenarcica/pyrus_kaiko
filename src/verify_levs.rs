///~/////////////////////////////////////////////
///  verify_levs.rs
///
///
///  Verify a NBBO sequence against arbitrary ordered series of resting price requests
///
///  Searches for situation where a break occurs and the resting order in question
///  either violateds NBBO sequence or comes within a pre-determined distance 'pdiff' of goal
///
use crate::sip_struct::{TT,TP};
use std::cmp::Ordering;

use std::thread;
use std::thread::scope;
//use itertools::Itertools;


///~///////////////////////////////////////////////////////////////////////
/// Safely subtract one from a usize
///
/// It can be annoying that we have to create safety checks a lot in terms of a try_from
/// Ideally a macro might save space here, though often we can safely subtract 1 from
/// this value
macro_rules!um1{
  ($u0:expr) => { match usize::try_from(if $u0 == 0 {0 as i64} else {($u0 as i64)-1}) {Ok(val)=>val,Err(_e)=>0}}
}

///~///////////////////////////////////////////////////////////
///  Print a vector
///
///  For convenience we'd like to be able to understand system state, and these are Rust-based
///    macros for printing vectors in a few cases (basic case one vector, then more
///    advanced cases where a new sort order is requested or we have other information to add
///
///  Note, you may need to sort in order according to another vector
///  Of course we get so complicated we wind up doing vectors in vectors or vectors of vectors
///  sorted by other vectors.
///
///  If we get a negative number we want a little more printed out in details
macro_rules!print_v {
  ($i0:expr, $v0:expr, $n0:expr) => {
      let v0len:usize = ($v0).len() as usize; let v0lenm1 = match usize::try_from( (v0len as i64)-1) {Ok(val)=>val,Err(_e)=>0};
      let nwant:usize = match usize::try_from(if $n0 < 0 { -1*($n0 as i64) } else { $n0 }) { Ok(val)=>val,Err(_e)=>0 };
      if nwant == 0 {
        println!("{}[len={}] zero lines print requested", $i0, v0len);
      } else if v0len <= 0 {
        println!("{}[len={}] -- Has zero lines.", $i0, v0len);
      } else if nwant < (v0len+1) {
        println!("{}[len={}]: [{},...,{}]", $i0, v0len, ($v0).iter().take(nwant).map(|x| x.to_string()).collect::<Vec<String>>().join(","),
               ($v0)[v0lenm1]);
      } else {
        println!("{}[len={}]: [{}]", $i0, $v0.len(), ($v0).iter().map(|x| x.to_string()).collect::<Vec<String>>().join(","));
      }
  };
  ($i0:expr,$v0:expr,$n0:expr,$srt0:expr,$type0:ty) => {
     print_v!($i0, (($srt0).iter().filter_map(|&idx| ($v0).get(idx).cloned()).collect::<Vec<$type0>>()),$n0);
  };
  ($i0:expr,$v0:expr,$n0:expr,$srt0:expr,$type0:ty,$othv0:expr) => {
     let nwant:usize = match usize::try_from(if $n0 < 0 { -1*($n0 as i64) } else { $n0 }) { Ok(val)=>val,Err(_e)=>0 };
     if $n0 >= 0 {
       print_v!($i0, (($srt0).iter().filter_map(|&idx| $v0.get($othv0[idx]).cloned()).collect::<Vec<$type0>>()),$n0);
     } else if $n0 < 0 {
       if nwant < (1+($srt0).len()) as usize {
         println!("{}[len={}]: [{},...,{}]", $i0, ($v0).len(), 
           ($srt0).iter().take(nwant).filter_map(
             |&idx| format!("[{}]{}", match $othv0.get(idx).cloned() { Some(val)=> val,None=>0 as usize}.to_string(), 
                                      match $v0.get($othv0[idx]).cloned() { Some(val)=>val,None=>0 as f64 }.to_string()
                         ).into()  ).collect::<Vec<String>>().join(","),
           format!("[{}]{}", $othv0[ $srt0[$srt0.len()-1] ], $v0[$othv0[ $srt0[$srt0.len()-1] ] ]   )); 
       } else {
         println!("{}[len={}]: [{}]", $i0, ($v0).len(), ($srt0).iter().filter_map(
           |&idx| format!("[{}]{}", $othv0[idx].to_string(), match $v0.get($othv0[idx]).cloned() { Some(val)=>val,None=>0 as f64 }.to_string()).into()  ).collect::<Vec<String>>().join(",")); 
       }
     }
  }
}

///~//////////////////////////////////////////////////////////////////////////////////////
///  Algorithm diagram. 
///
///
///    Consider a sequence of prices "NBBOmid[t]" and we look for its intersection
///     of a price level "P
///
///   |
///   |__                 ___               _________      NBBOmid[t]
///   |  |   __       ___|   |            _|         |___/
///   |  |__|  |     |       |           | 
/// -P|        X     X       XXX        XX
///   |        |__   |         |_      _| \
///   |           |__|           |   _|    "X" for all crossings"
///   |                          |__|
///   |
///   |_______________________________________________________
///                         time
///
///
///  In our algoirthm we caculate for a given price level the known crossing points of the target
///  price P
///  We check these finite number of points against all resting orders resting at P.
///  Assuming the NBBOmid[t] "breaks" through we know the order is hit, also if it "rests" on the
///  order, the order "may" have been hit, depending on whether we will accept "instantaneous"
///  touching or not (say an order closes exactly at the point the NBBOmid[t] level reaches it)
///
///  This algorithm should be efficient because "updates" to P+1 are easy.  We first would check
///  whether it is time to dump "existing cross events", or else add maybe a new "dip" that first
///  hits.  Finding the additional points to add can be easy if we have sorted all NBBO mid events
///
///  Note that the "NBBOmid[t]" moves as a right angled step function.  Prices change "immediately"
///  on an update.  But it is possible that NBBOmid can jump 2 or more levels ina  single update.
///   (This would be if NBBO goes from $100.02 from $100.00, skipping $100.01)
///

///~//////////////////////////////////////////////////////////////
/// fast_alter_cross_above (delicate fix assumption)
///
/// Here we assume that the ability to move up in price levels from lowest price level to highest
/// in v_p vector, which will cover increasing levels of prices in v_nop NBO price.
///
/// We assume that as we move a step up to the next level set, we assume that
/// we can 
///   1. only keep or remove existing crossers
///     Those are previous (+2) for fully under to fully over to (-2) fully over to fully under
///   2.  If we were previous below to equal (+1) we are definitely being removed
///   3.  If we were equal moving to below (-1) we are possibly becoming a (-2)
///   4.  If we were above to equal we are possibly becoming (-2)
///   5.  If we were equal to above we are possibly becoming (+2)
///
///   So some of the previous will go away or become stronger moves.
///   However, there can be none more from the set of v_nop[ik] for ik < previous ip_eq.
///   This is because we know that none of those other values jumped up to previous ptgt
///   And for the new ptgt the next answer can't be that big.
fn fast_alter_cross_above(ptgt: TP, st_n:&[usize], v_nop:&[TP], 
  out_x:&mut Vec<usize>, 
  ip_eq: &mut usize, ip_u: &mut usize) -> usize {
  //let st_n = *st_n; let v_nop = *v_nop;
  let n_n = st_n.len();
  //let mut t_out_x:Vec<usize> = Vec::with_capatcity(out_x.len());
  let mut write_i: usize = 0;
  for on_i in 0_usize..out_x.len() {
    let on_i:usize = on_i as usize;
    if out_x[on_i]+1 >= n_n {
      println!("fast_alter_cross_above: Error out_x[{}] ={} versus n_n={}", on_i, out_x[on_i], n_n);
    }
    if v_nop[out_x[on_i]] > ptgt {
      if (out_x[on_i] + 1 < n_n) && (v_nop[out_x[on_i]+1] <= ptgt) {
          out_x[write_i] = out_x[on_i]; write_i=write_i+1;
      } 
    } else if v_nop[out_x[on_i]] == ptgt {
      out_x[write_i] = out_x[on_i]; write_i=write_i+1; 
    } else {
      if (out_x[on_i] +1 < n_n) && (v_nop[out_x[on_i]+1] >= ptgt) {
        out_x[write_i] = out_x[on_i]; write_i=write_i+1; 
      }
    }
  }
  //if verbose >= 1 {
  //  println!("alter_cross_above: shrinking out_x from {} to {}. ", out_x.len(),write_i);
  //}
  out_x.truncate(write_i); // Everyone else is removed as they are no longer valid skips
  //if verbose >= 1 {
  //  println!("alter_cross_above: after that out_x is length {}.", out_x.len());
  //}
  // Next added memvers must be from 
  *ip_u = *ip_eq;
  if *ip_u > 0 { *ip_u = *ip_u-1; };
  while ( (*ip_u) < n_n) && (v_nop[st_n[ *ip_u ]] <= ptgt) {
    let on_i:usize = *ip_u as usize; let onp = st_n[on_i];
    if (onp as usize) >= (1 as usize) {
      let pp:usize = um1![onp];
      if v_nop[pp] >= ptgt {
        if (v_nop[onp] <= ptgt) && !(out_x.contains(&pp)) { 
          out_x.push(pp); write_i=write_i+1;
        }
      }
    }
    if v_nop[onp] < ptgt {
      if (onp+1 < n_n) && (v_nop[onp+1] >= ptgt) && !(out_x.contains(&onp)) {
        out_x.push(st_n[on_i]); write_i=write_i+1;
      } 
    } else if v_nop[onp] == ptgt {
      if !(out_x.contains(&onp)) {
        out_x.push(onp);  write_i = write_i+1;
      }
    }
    if (v_nop[st_n[*ip_u]] < ptgt) &&  (v_nop[st_n[*ip_u+1]] >= ptgt) { *ip_eq = *ip_u+1 }
    (*ip_u) = (*ip_u) + 1;
  } 
  let try0:usize = 0;
  if (v_nop[try0] <= ptgt) && !out_x.contains(&try0) {
    out_x.push(try0);
  }
  return out_x.len();
}


///~/////////////////////////////////
/// slow_calculate_cross_above:
///
///  This works to identify all distinct crossings of price "ptgt" using the fact
///  We know this data is sorted by price than time using price then time using st_n
///  But that it is sorted by time (then we assume no two prices can occur at same time)
///
///  We assume that the number of prices less than ptgt is a small subset of all prices.
///  However, we could always start at minimum price and jump to a higher price above the ptgt at
///  anytime.
///
///  Thus we have to walk all prices up to "ip_u"  (but this number is hopefully
///  Small.)  By update, this number is likely close to a previous estimate.
///
///  For every price point, we look at tick before and after.
///  If we can tell we just jumped from higher prices to inside, we push onto out_x the index
///  where we know next time is a cross.  This is a "-2" event (outside to inside).
///  If the previous was "equal" to ptgt, we push the event, but only push a -1 (on border to
///  inside).
///
///  If the 
///
fn slow_calculate_cross_above(ptgt: TP, st_n:&[usize], v_nop:&[TP], 
  out_x:&mut Vec<usize>, 
  ip_eq: &mut usize, ip_u: &mut usize) -> usize {
  let n_n = st_n.len();
  out_x.clear(); 
  while ( (*ip_eq) < n_n) && (v_nop[st_n[ *ip_eq ]] < ptgt) {
    (*ip_eq) = (*ip_eq) + 1;
  } 
  if (*ip_eq) >= n_n {  
    return 0;
  }
  *ip_u = *ip_eq;
  while (*ip_u < n_n) && (v_nop[st_n[*ip_u]] <= ptgt) {
    *ip_u = *ip_u + 1;
  }
  for on_i in 0..*ip_eq {
    let on_i:usize = on_i as usize;
    if ( st_n[on_i] as usize) >= (1 as usize) {
        let pp:usize = um1![ st_n[on_i] ];
        if v_nop[pp] > ptgt {
          out_x.push(pp); 
        } else if v_nop[pp] == ptgt {
          out_x.push(pp);  
        }
    }
    if (st_n[on_i] + 1) < n_n {
       let np:usize = match usize::try_from( st_n[on_i] + 1 ) {
          Ok(val) => val, Err(_e) => 0 };
       if v_nop[np] > ptgt {
         out_x.push(st_n[on_i]); 
       } else if v_nop[np] == ptgt {
         out_x.push(st_n[on_i]);  
       }
    }
  }
  for on_i in (*ip_eq)..(*ip_u) {
    if ( st_n[on_i]) >= (1 as usize) {
      let pp:usize = um1![ st_n[on_i]];
      if v_nop[pp] > ptgt {
        out_x.push(pp); 
      } 
    }
    if (st_n[on_i] + 1)< n_n {
       let nextp:usize = st_n[on_i] + 1; 
       if v_nop[nextp] > ptgt {
         out_x.push(st_n[on_i]);  
       } 
    }
  }
  let try0:usize = 0;
  if (v_nop[try0] <= ptgt) && !out_x.contains(&try0) {
    out_x.push(try0);
  }
  return out_x.len();
}


// Generates out_idx sort of v_s/v_p_vot by side/price/time
pub fn sort_v_pt(v_s:&[i8], v_p:&[TP], v_ot:&[TT]) -> Vec<usize> {
  let mut out_idx:Vec<usize> = (0..v_p.len()).collect();
  out_idx.sort_by(|i,j| if v_s[*i] < v_s[*j] { Ordering::Less 
                        } else if v_s[*i] > v_s[*j] { Ordering::Greater
                        } else if v_p[*i] < v_p[*j] { Ordering::Less
                        } else if v_p[*i] > v_p[*j] { Ordering::Greater
                        } else if v_ot[*i] < v_ot[*j] { Ordering::Less
                        } else if v_ot[*i] > v_ot[*j] { Ordering::Greater
                        } else { Ordering::Equal });
  //v_p[*i].cmp(&v_p[*j]).then_with(|| v_ot[*i].cmp(&v_ot[*j]))); 
  return out_idx;
}

// Sorts the NBBO price level by price/time
fn sort_nbbo_pt(v_nop:&[TP], v_not:&[TT]) -> Vec<usize> {
  let mut out_idx:Vec<usize> = (0_usize..(v_nop.len() as usize)).collect();
  out_idx.sort_by(|i,j| if v_nop[*i] < v_nop[*j] { Ordering::Less
                        } else if v_nop[*i] > v_nop[*j] { Ordering::Greater
                        } else if v_not[*i] < v_not[*j] { Ordering::Less
                        } else if v_not[*i] > v_not[*j] { Ordering::Greater
                        } else { Ordering::Equal });
  //out_idx.sort_by(|i,j| v_nop[*i].cmp(&v_nop[*j]).then_with(|| v_not[*i].cmp(&v_not[*j]))); 
  return out_idx;
}
// We do this sort because we sometimes don't calculate cross points in exact order they occur
fn resort_out_x(out_x: &Vec<usize>) -> Vec<usize> {
  let mut out_idx:Vec<usize> = (0_usize..(out_x.len() as usize)).collect();
  out_idx.sort_unstable_by(|i,j| out_x[*i].cmp(&out_x[*j]) );
  return out_idx;
}

///~//////////////////////////////////////////////////////////////////////
/// "fast_alter_cross_below" is the "fast update" companion to slow_calculate_cross_below()
///
/// One complete calculates the series fresh, the other assumes we are coming from a nearby price 
/// level and acts accordingly.
///
///
fn fast_alter_cross_below(ptgt: TP, st_n:&[usize], v_nop:&[TP], 
  out_x:&mut Vec<usize>, 
  ip_eq: &mut usize, ip_u: &mut usize) -> usize {
  let n_n = st_n.len();
  let mut write_i: usize = 0;

  // Step: Eliminate any crossing points from existing list that are no longer valid
  for on_i in 0_usize..out_x.len() {
    let on_i:usize = on_i as usize;
    //if out_x[on_i]+1 >= n_n {
    //  println!("fast_alter_cross_below: Error out_x[{}] ={} versus n_n={}", on_i, out_x[on_i], n_n);
    //}
    if v_nop[out_x[on_i]] > ptgt {
      if (out_x[on_i] + 1 < n_n) && (v_nop[out_x[on_i]+1] <= ptgt) {
        out_x[write_i] = out_x[on_i]; write_i=write_i+1;
      }
    } else if v_nop[out_x[on_i]] == ptgt {
      out_x[write_i] = out_x[on_i]; write_i=write_i+1; 
    } else {
      if (out_x[on_i] + 1 < n_n) && (v_nop[out_x[on_i]+1] >= ptgt) {
        out_x[write_i] = out_x[on_i]; write_i=write_i+1; 
      }
    }
  }
  //if verbose >= 2 {
  //  println!("{} -- We have written {} but had {} originally in out_x", vstr, write_i, out_x.len());
  //}
  out_x.truncate(write_i); // Everyone else is removed as they are no longer valid skips
  // Next add new members that are known to just "exactly" hit our price level of interest. 
  *ip_u = if (*ip_eq+1) < n_n { *ip_eq+1 } else { *ip_eq };
  while v_nop[st_n[ *ip_u ]] >= ptgt {
    let on_i:usize = *ip_u as usize; let onp = st_n[on_i];
    let pp:usize = um1![onp];
    //if verbose >= 3 {
    //  println!("{} --- on_i={}, st_n[{}]={}, v_nop[{}]={} versus ptgt={}, on_condition={}.",
    //    vstr, on_i, on_i, st_n[on_i], onp, v_nop[onp], ptgt,
    //    if v_nop[onp] == ptgt { "NOW_EQUAL" } else if v_nop[onp] < ptgt { "NOW_LESS" }
    //    else if v_nop[pp] < ptgt { "WAS_LESS"} else if v_nop[pp] == ptgt { "WAS_EQUAL" }
    //   else if (v_nop[onp] > ptgt) && (onp + 1< n_n) && (v_nop[onp+1] <= ptgt) { "GOING_LESS" }
    //    else { "SOMETHING_ELSE } "}
    //    );
    //}
    if (onp as usize) >= (1 as usize) {
      if (v_nop[pp] <= ptgt) && !out_x.contains(&pp) {
        if v_nop[onp] >= ptgt { out_x.push(pp); }
      }
    }
    if v_nop[onp] > ptgt {
      if (onp+1 < n_n) && (v_nop[onp+1] <= ptgt) && !out_x.contains(&onp) {
         out_x.push(onp); write_i = write_i+1;
      } 
    } else if (v_nop[onp] == ptgt) && !out_x.contains(&onp) { 
       out_x.push(onp); write_i = write_i + 1;
    } else if v_nop[onp] < ptgt {
      if (onp+1 < n_n) && (v_nop[onp+1] >= ptgt) && !out_x.contains(&onp) {
         out_x.push(onp); write_i = write_i + 1;
      }
    }
    let next_i:usize = match usize::try_from((*ip_u)-1) { Ok(val)=>val,Err(_e)=>0 };
    if (v_nop[st_n[*ip_u]] > ptgt) &&  (v_nop[st_n[next_i]] >= ptgt) { *ip_eq = next_i; }
    (*ip_u) = next_i; 
    if *ip_u == 0 { break; }
  } 

  let try0:usize = 0 as usize;
  if (v_nop[try0] >= ptgt) && !out_x.contains(&try0) {
    out_x.push(try0);
  }
  return out_x.len();
}


///~////////////////////////////////////////////////////////////////////
///  "slow_calculate_cross_below"
///   This is used to get crossing events absent any previous information.
///   This is good for checking quality of the fast_alter_cross_below calculation
///
fn slow_calculate_cross_below(ptgt: TP, st_n:&[usize], v_nop:&[TP], 
  out_x:&mut Vec<usize>, 
  ip_eq: &mut usize, ip_u: &mut usize) -> usize {
  //let st_n = *st_n; let v_nop = *v_nop;
  let n_n:usize  = st_n.len() as usize;
  out_x.clear(); 
  while ( (*ip_eq) > 0) && (v_nop[st_n[ match usize::try_from( (*ip_eq as i64)-1) { Ok(val)=>val,Err(_e)=>0} ]] > ptgt) {
    (*ip_eq) = match usize::try_from( (*ip_eq as i64) -1 ) { Ok(val)=>val, Err(_e)=>0 };
  } 
  if ((*ip_eq) == 0)  && (v_nop[st_n[ *ip_eq]] > ptgt) {  
    return 0;
  }
  *ip_u = *ip_eq;
  while (*ip_u > 0) && (v_nop[st_n[ match usize::try_from( (*ip_u as i64)-1) { Ok(val)=>val, Err(_e)=>0} ]] >= ptgt) {
    let n_val:i64 = (*ip_u as i64) - (1 as i64);
    (*ip_u) = match usize::try_from( n_val ) { Ok(val)=>val, Err(_e)=>0 };
    //(*ip_u) = <i64 as Into<usize>>::n_val.into().expect("slow_calculate_cross_below: This subtraction should never be less than zero"); 
  }
  for on_i in (((*ip_eq))..n_n).rev() {
    if ( st_n[on_i] as usize) >= (1 as usize) {
        let pp:usize = um1![st_n[on_i]];
        if v_nop[pp] < ptgt {
          out_x.push(pp); 
        } else if v_nop[pp] == ptgt {
          out_x.push(pp);  
        }
    }
    if ((st_n[on_i] + 1) as usize) < n_n {
       let nextp:usize = st_n[on_i] + 1;
       if v_nop[nextp] < ptgt {
         out_x.push(st_n[on_i]);  
       } else if v_nop[nextp] == ptgt {
         out_x.push(st_n[on_i]);  
       }
    }
  }
  if (*ip_eq) == (*ip_u) {
    return out_x.len();
  }
  for on_i in ((*ip_u)..(*ip_eq)).rev() {
    let on_i:usize = on_i as usize;
    if ( st_n[on_i]) >= (1 as usize) {
      let pp:usize = um1![st_n[on_i]];
      if v_nop[pp] < ptgt {
        out_x.push(pp); 
      } 
    }
    if ((st_n[on_i] + 1) as usize) < n_n {
       let nextp:usize = st_n[on_i]+1;
       if v_nop[nextp] < ptgt {
         out_x.push(st_n[on_i]);  
       } 
    }
  }

  // Odd case where the first point is technically a point we should check
  let try0:usize = 0;
  if (v_nop[try0] >= ptgt) && !out_x.contains(&try0) {
    out_x.push(try0);
  }
  return out_x.len();
}

///~//////////////////////////////////////////////////////////////////////////////////////////////
/// basic_hms
///   "Basic Hours Min Second".  This is used to print out time in more human readable units
///
///   Basic string generation of hours/min/second from i64 data ignoring the year/date/month
fn basic_hms(itime:i64) -> String {
   let frac:i64 =  itime %     1000000000;
   let secs:i64 = (itime %    60000000000) /  1000000000;
   let mins:i64 = (itime %  3600000000000) / 60000000000;
   let  hrs:i64 = (itime % (24 * 3600 * 1000000000)) / (3600 * 1000000000);
   let  hrs = if hrs.to_string().len()==0 {format!("00")} else if hrs.to_string().len() == 1 {format!("0{}",hrs.to_string()) } else {hrs.to_string()};
   let mins = if mins.to_string().len()==0 {format!("00")} else if mins.to_string().len() == 1 {format!("0{}",mins.to_string()) } else {mins.to_string()};
   let secs = if secs.to_string().len()==0 {format!("00")} else if secs.to_string().len() == 1 {format!("0{}",secs.to_string()) } else {secs.to_string()};
   let frac = if frac.to_string().trim_end_matches('0').len() == 0 { format!("") } else {format!(".{}",frac.to_string().trim_end_matches('0'))};
   return format!("{}:{}:{}{}", hrs, mins, secs, frac);
}

///~//////////////////////////////////////////////////////////////////////////////////////////////
///  full_check_jumps
///
///  "Lazy Complete Search algorithm"
///
///  Checks for jumps the slowest possible way, sweeping completely thru price levels.
///  The result of check should only capture the moment when we "jump" from one zone to another.
///
///  Note, there is an edge case, depending on side for whether the it=0th or first moment counts
///  as a checkable event.  In general, we are hoping never to have to reject on it=0 since
///  that is when time window opens and most order messages don't exist (and hopefully they don't
///  open in violation)
///
///  This is a slow and most complete check, walking the entire v_nop price time series for all
///  times this series crosses ptgt (or lands on it exactly)
fn full_check_jumps(onstr: &str, ptgt:TP, v_nop:&[TP],sort_cross_idx:&Vec<usize>,cross_idx:&Vec<usize>) -> u64 {
  let np = v_nop.len(); let mut nerr:u64=0;
  let npm1:usize = um1![np];
  let mut cm_idx:Vec<usize> = Vec::new();
  for it in 0..npm1 {
    let it:usize=it as usize;
    if v_nop[it] > ptgt {
      if v_nop[it+1] < ptgt {
        cm_idx.push(it); 
      } else if v_nop[it+1] == ptgt {
        cm_idx.push(it); 
      }
    } else  if v_nop[it] < ptgt {
      if v_nop[it+1] > ptgt {
        cm_idx.push(it); 
      } else if v_nop[it+1] == ptgt {
        cm_idx.push(it); 
      }
    } else if v_nop[it] == ptgt {
      if v_nop[it+1] > ptgt {
        cm_idx.push(it); 
      } else if v_nop[it+1] < ptgt {
        cm_idx.push(it); 
      } // We consider "same price" should be #1 impossible, #2 ignorable
    } 
  }
  if (cm_idx.len() == 0) && (cross_idx.len() == 0) {
    // Well, all clear
    return 0;
  } else  if (cm_idx.len() > 0) && (cross_idx.len() == 0) {
    println!("full_check_jump->{} -- ERROR Found crosses versus Zero cm_idx len={}, cross_idx len={}", onstr, cm_idx.len(), cross_idx.len());  nerr = nerr+1;
  }
  let mut d1:usize = 0;
  if cm_idx.len() == um1![cross_idx.len()]  && (cm_idx[0] != 0) && (cross_idx[sort_cross_idx[0]] == 0) {
    // Edge case where we stuck an early zero because it is an "edge case check" (but you need
    // side input to know this.
    d1 = 1 as usize;
  } else  if cm_idx.len() != cross_idx.len() {
      println!("full_check_jump->{} -- ERROR, cm_idx len={}, cross_idx len={}", onstr, cm_idx.len(), cross_idx.len());  nerr = nerr+1;
  }
  let nmin= if cm_idx.len() < cross_idx.len() { cm_idx.len() } else { cross_idx.len() };
  for ix in 0..nmin {
      let ix:usize = ix as usize;
      if cm_idx[ix] != cross_idx[sort_cross_idx[ix+d1]] {
        println!("full_check_jump->{} Error on ix={}/{}, {}, we had cm_idx[{}] = {}, but cross_idx[{}=st[{}]]={}.",
          onstr, ix, nmin, if d1 == 1 { "cross_idx[0]=0 so jump 1" } else {""},
          ix, cm_idx[ix], sort_cross_idx[ix+d1], ix+d1, cross_idx[sort_cross_idx[ix+d1]]);  nerr = nerr+1;
        println!(" --- According to our issues v_nop[cm_idx[{}]={}]={}, and v_nop[(cm_idx[{}]+1)={}]={}, note ptgt={}", 
          ix, cm_idx[ix], v_nop[cm_idx[ix]], ix, cm_idx[ix]+1, v_nop[if (cm_idx[ix]+1)<np { cm_idx[ix]+1 } else {cm_idx[ix]}], ptgt);
      } 
  }
  if nerr == 1 {
    println!("Weird error must have been on a neq.");
    println!("  Last cm_idx[{}] = {}, for v_nop[{}]={} and v_nop[{}]={}.",
      cm_idx.len()-1, cm_idx[cm_idx.len()-1], cm_idx.len()-1, v_nop[cm_idx[cm_idx.len()-1]], cm_idx[cm_idx.len()-1] + 1,
        v_nop[cm_idx[cm_idx.len()-1] + 1]); 
  }
  println!("full_check_jump->{} full_check_jump finished with {} errors.", onstr, nerr);
  return nerr;
}
// Note  for check_jumps you need to pick appropriate target
///~///////////////////////////////////////////////////////////////////////////////////////////////
///  check_jumps()
///
///  TEST UTILITY
///
///  Check our calculated "cross_idx" vector and verifies that the events it calculated really do
///  jump accross the price level ptgt.
fn check_jumps(onstr: &str, ptgt:TP, v_nop:&[TP],sort_cross_idx:&Vec<usize>,cross_idx:&Vec<usize>) -> u64 {
  let np = v_nop.len(); let nidx = sort_cross_idx.len();
  let mut nerr: u64 = 0;
  for ix in 0..nidx {
    let ix:usize = ix as usize;
    let oncross_idx:usize = cross_idx[sort_cross_idx[ix]];
    let next_idx:usize = oncross_idx + 1;
    let on_price = v_nop[cross_idx[sort_cross_idx[ix]]];
    if (ix == 0) && (cross_idx[sort_cross_idx[ix]] == 0) {
      // Edge case, yes, this is technically only good if the data is on "wrong side" at t=0, but
      // we don't really need to check this one.
    } else if next_idx >= np {
      println!("check_jumps: Really awkward error {}, next_idx = {} but np={}, we were ix={}, oncross_idx={}, ptgt={}",
         onstr, next_idx, np, ix, oncross_idx, ptgt); nerr = nerr + 1;
    } else {
      let next_price = v_nop[next_idx];
      if (on_price > ptgt) && (next_price < ptgt) {
      } else if (on_price > ptgt) && (next_price == ptgt) {
      } else if (on_price == ptgt) && (next_price < ptgt) {
      } else if (on_price == ptgt) && (next_price > ptgt) {
      } else if (on_price < ptgt) && (next_price > ptgt) {
      } else if (on_price < ptgt) && (next_price == ptgt) {
      } else {
           println!("check_jumps: ERROR {}, we ix={}, price location {}, jump from {} to {} but ptgt={}",
                  onstr, ix, oncross_idx, on_price, next_price, ptgt); nerr = nerr + 1;
      }

    }

  } 
  return nerr;
}
///~////////////////////////////////////////////////////
/// Our methodology to verify that every order at price v_p[i], with open and close times
/// is satisfactory against a price level + delta of v_nop-pdiff for every price in v_no[
///
/// What we do is figure out a vector of all known "Next stop is a switch" points.
/// We generate this with slow_calculate_cross_above()/fast_alter_cross_above() 
///
/// This gives us a logical list of all points, which, if the open/close times overlap one of
/// these switch points, then we can make a determination if the
///
/// Note if "eqstate==True" then any points "barely touching" are rejectable.
/// This makes it easy to say that if a touch point even exists, then we reject.
/// We have a little hard time if eqstate==False and we can tell we just start to sit beside
/// the price.  Even so, we can determine by cross_type whether this spot eventually jumps
/// back to clear (positive jump is next), or jumps to danger zone (negative jump comes next) and
/// if the time of that jump is also within the open_time/close_time interval
///
/// The resulting return vector "out_v" will be 0 if we never rejected the point, and 1 if we found
/// one or more reasons to reject a point (once we find at least one reason we are done).
///
/// Mathematically we need to calculate for every resting order v_p[i],
///    lasting from v_ot[i] to v_ct[i], is there a price event v_nop[j] occuring at v_not[j]
///    with v_ot[i] <= v_not[j] <= v_ct[i],
///    Such that v_nop[j] - pdiff  < v_p[i];
///
///  This is typically a n_v={Length v_p} times n_n {Length v_nop} search.
///     But this result would take too much time and energy for large vectors.
///
///  Inputs
///    pdiff: A price difference.  We want to know if there is a cross between v_nop[j]-pdiff and
///    v_p[], v_ot[], v_ct[]: Vectors for resting messages at price p, open ot, close ct.
///    v_nop, v_not: The continuous NBB or NBO price sequence we use as barrier price
///    verbose: Whether to print details
///    "eqstate":  We could use a "Kissing cutoff" (if eqstate==True) or a "breakthrough cuttoff"
///        So our condition of "<" or "<=" as a rejector.  Typically we will need a full
///        break through.  But because there is a potential for v_nop to "bounce and touch"
///        We have to consider tighter criteria
///
///
pub fn verify_buy(out_v:&mut Vec<u64>, st_v: &Vec<usize>, pdiff:f64, v_s:&[i8], v_p: &[TP], v_ot:&[TT], v_ct:&[TT],
              v_nop: &[TP], v_not: &[TT], eqstate:bool, verbose:i8, timetouch:bool) -> i8 {
  //let v_nop = *v_nop; let v_not = *v_not;
  let n_v:usize = v_p.len(); let n_n:usize = v_nop.len(); let n_n64:u64 = match u64::try_from(n_n) { Ok(val)=>val,Err(_e)=>0 };

  let vstr = format!("verify_buy({},n_v={},n_n={},tt={})",pdiff, n_v, n_n, if timetouch {1} else {0});
  if verbose >= 1 { println!("{}: We begin", vstr); }
  if n_v <= 0 {  return -1; } // No buys to match, return nothing. 
  if n_n <= 0 {  return -1; } // No NBBO cant do anything. 

  let n_buy = find_1(st_v, v_s);
  if n_buy <= 0 { return 0; }
  //let st_v = sort_v_pt(&v_s, &v_p, &v_ot);
  let st_n_v = sort_nbbo_pt(v_nop,v_not);
  let st_n = &st_n_v[..];
  if verbose >= 2 {
    println!("{} --- We are beginning with the vectors we have sorted. ", vstr);
    print_v!(format!("{}: st_v,", vstr), st_v,5); print_v!(format!("{}: st_n,", vstr), st_n,5);
  }
  //let mut out_v: Vec<u64> = vec![n_n as u64;n_v as usize];
  let mut on_plev = v_p[st_v[0]] + pdiff;
  let mut cross_idx:Vec<usize> = Vec::new();
  let mut on_ix = 0;
  
  let mut ip_eq:usize = 0; let mut ip_u: usize = 0;
  let mut nx = slow_calculate_cross_above(on_plev, st_n, v_nop, 
    &mut cross_idx, &mut ip_eq, &mut ip_u);
  let mut sort_cross_idx = resort_out_x(&cross_idx);
  if verbose >= 2 {
     println!("{} - for first price on_plev={}={}[{}=st_v[{}]]-{} we have ip_eq={}, ip_u={}",
       vstr, on_plev, v_p[st_v[0]], st_v[0], 0, pdiff, ip_eq, ip_u);
     print_v!(format!("{}: sort_cross_idx",vstr), sort_cross_idx, 5);
     print_v!(format!("{}: cross_idx",vstr), cross_idx, 5, sort_cross_idx, usize); 
     print_v!(format!("{}: v_nop[cidx]",vstr), v_nop,5,sort_cross_idx,TP,cross_idx); 
  }

  let n_nm1:usize = um1![n_n]; 
  for i_v in (0 as usize)..n_v {
    let i_v: usize = i_v as usize;  let v_on_p = v_p[st_v[i_v]]; let v_on_ct = v_ct[st_v[i_v]]; 
    if v_s[st_v[i_v]] == 1 { break; }
    let v_on_ot = v_ot[st_v[i_v]];
    if (v_on_p+pdiff) > on_plev {
      on_plev = v_on_p + pdiff;
      // Note "update_cross_above(...) will start from zero no matter price target.
      //   It will typically start at "minimum price" however, so it goes up.
      //
      // Fast Alter should be faster that just
      // alters existing solution
      nx = fast_alter_cross_above(on_plev, st_n, v_nop, 
        &mut cross_idx, &mut ip_eq, &mut ip_u);
      sort_cross_idx = resort_out_x(&cross_idx);  // We need to move forward through cross_idx in time.
      on_ix = 0;
      if verbose >= 2 {
        println!("{} - for renewed price on_plev={}={}[{}=st_v[{}]]-{} we have ip_eq={}, ip_u={}",
          vstr, on_plev, v_p[st_v[i_v]], st_v[i_v], i_v, pdiff, ip_eq, ip_u);

        print_v!(format!("{}: sort_cross_idx",vstr), sort_cross_idx, 5);
        print_v!(format!("{}: cross_idx",vstr), cross_idx, 5, sort_cross_idx, usize); 
        print_v!(format!("{}: v_nop[cidx]",vstr), v_nop,5,sort_cross_idx,TP,cross_idx); 
        let nerr = check_jumps(&format!("verify_buy[i_v={}]",i_v), on_plev, v_nop, &sort_cross_idx, &cross_idx);
        let nerr2 = full_check_jumps(&format!("verify_buy[i_v={}]",i_v), on_plev, v_nop, &sort_cross_idx, &cross_idx);
        println!("{} -- We checked jumps and got {} errors for {} finds, on_plev={} n_n={}, nerr2={}", 
          vstr, nerr, nx, on_plev, n_n, nerr2);  if nerr > 0 { println!("{} ERROR we know we aren't doing good. ip_eq={}, ip_u={}", vstr, ip_eq, ip_u); }
      }
    }
    if verbose >= 3 {
        println!("\n\n{}: on (i_v={}/{}=st_v={}), v_p[{}]={}, ot--ct=[{}--{}] ",
           vstr, i_v, st_v.len(), st_v[i_v], st_v[i_v], v_p[st_v[i_v]], basic_hms(v_ot[st_v[i_v]]), basic_hms(v_ct[st_v[i_v]])); 
        print_v!(format!("      : sort_cross_idx"), sort_cross_idx, 5);
        print_v!(format!("      : cross_idx"), cross_idx, 5, sort_cross_idx, usize); 
        print_v!(format!("      : v_nop[cidx]"), v_nop,5,sort_cross_idx,TP,cross_idx); 

    }
    if nx == 0 {
      if verbose >= 2 {
        println!("{} -- in this case nx={}, ip_eq={}/{}, on_plev={}, but local price={} for pdiff={}, v_nop ranges = [{},{}], what do we detect?",
          vstr, nx, ip_eq, n_n, on_plev, v_p[st_v[i_v]], pdiff, v_nop[st_n[0]], v_nop[st_n[n_nm1]]);
      }
      if ip_eq >= n_n {
        // In this case all prices are always less if there is any overlap.  Reject
        if (v_on_ot <= v_not[um1![n_n]]) && (v_on_ct >= v_not[0]) {
          out_v[st_v[i_v]] = 0 as u64;
        } else {
          // Crazy case everyone out of range
          out_v[st_v[i_v]] = n_n64;
        }
      } else if (v_nop[st_n[0]] - pdiff) > on_plev {
        // In this case all prices are north of issue.  We should keep this point.
        out_v[st_v[i_v]] = n_n64;
      } else {
        println!("{}  ISSUE nx = 0, for v_nop range [{},{}], ip_eq={}/{}, on_plev={} but pdiff={}",
          vstr,  v_nop[st_n[0]], v_nop[st_n[n_nm1]], ip_eq, n_n, on_plev, pdiff);
      }
    } else {
       // Note: our logic for determining a "full/partial"touch could be better organized, but
       // there is some complexity as to what counts as a "damning cross" versus a "near miss"
       for i_x in (on_ix as usize)..(nx as usize) {
         let i_x:usize = i_x as usize;
         let v_on_not = v_not[cross_idx[sort_cross_idx[i_x]]];
         let v_on_not_next = if i_x+1 < nx { v_not[cross_idx[sort_cross_idx[i_x+1]]] } else { v_on_not };
         let current_n = cross_idx[sort_cross_idx[i_x]];
         let next_n = if (current_n+1) < n_n { current_n+1 } else {current_n};
         if v_on_not < v_on_ot {
            if  (i_x +1 < nx) && (v_on_not_next < v_on_ot) {
              on_ix = on_ix + 1;
            } else if (timetouch==false) && (v_on_not_next == v_on_ot) {
            } else {
              // If we are at last break before v_ot we have to as questions.
              let lastn:usize = if (i_x+1) < sort_cross_idx.len() { cross_idx[sort_cross_idx[i_x+1]] } else { v_not.len() };
              for on_n in current_n..lastn {
                let on_n:usize = on_n as usize; let next_n:usize = if on_n+1 < v_not.len() { on_n+1 } else { on_n };
                if (v_not[on_n] > v_ct[st_v[i_v]]) || ((timetouch==false) && (v_not[on_n] >= v_on_ct)) { 
                  break;
                } else if (v_not[next_n] > v_ot[st_v[i_v]]) || ((timetouch==true) && (v_not[next_n] >= v_on_ot)) {
                  if (v_nop[on_n]-pdiff < v_p[st_v[i_v]]) || ((eqstate==true) && (v_nop[on_n]-pdiff <= v_p[st_v[i_v]])) {
                    out_v[st_v[i_v]] = on_n as u64; break;
                  }
                }
              }
            }
         } else if (v_on_not > v_ct[st_v[i_v]]) || ((timetouch==false) && v_on_not == v_ct[st_v[i_v]]) {
           break;
         } else {
           if (v_nop[current_n] - pdiff < v_on_p) || ((eqstate) && (v_nop[current_n]-pdiff <= v_on_p)) {
             out_v[st_v[i_v]] = current_n as u64;
           } else if (v_not[next_n] < v_on_ct) || ((timetouch) && (v_not[next_n] <= v_on_ct)) {
              if (v_nop[next_n]-pdiff < v_on_p) || ((eqstate==true) && (v_nop[next_n]-pdiff <= v_on_p)) {
                out_v[st_v[i_v]] = next_n as u64;
              }
           }
         }
         if out_v[st_v[i_v]] < n_n64 { break; }
       }
       if out_v[st_v[i_v]] >= n_n64 {
         out_v[st_v[i_v]] = n_n64;  // We have passed all challenges.
       }
    }
    if verbose >= 2 {
      println!("{} -- We reached i_v={}/{}, setting out_v[{}=st_v{}] = {}",
        vstr, i_v, n_v, st_v[i_v], i_v, out_v[st_v[i_v]]);
    }
  }
  return 1;
}

// Searches for first v_s=1 value, or the beginning of our sells
pub fn find_1(st_v: &Vec<usize>, v_s:&[i8]) -> usize {
  let mut i0 = 0; if v_s[st_v[i0]] == 1 { return i0; }
  let mut i1 = um1![st_v.len()]; if v_s[st_v[i1]] == 0 { return st_v.len(); }
  loop {
    if i0 == i1 { return i0; }
    let tp:usize  = match usize::try_from( (i0 + i1) / 2) { Ok(val)=>val,Err(_e)=>um1![i1] };
    if tp == i1 { if v_s[st_v[um1![tp]]] == 0 { return tp; } else { i1 = um1![i1] }
    }  else if tp == i0 { if v_s[st_v[tp+1]] == 1 { return tp+1; } else {i0 = tp+1; }
    }  else if v_s[st_v[tp]] == 1 { i1 = tp; if v_s[st_v[um1![tp]]] == 0 { return tp; }
    }  else if v_s[st_v[tp]] == 0 { i0 = tp;
    }
  }
}

/// parallel_verify_buy: Same logic as above, but splits data into to n threads.
///
/// Note that out_v could potentially be any point in big vector.  We might as well just give
/// chunks of this data, but we got lazy with mutexes.
pub fn parallel_verify_buy(out_v:&mut Vec<u64>, st_v: &Vec<usize>,  pdiff:f64, v_s:&[i8], v_p: &[TP], v_ot:&[TT], v_ct:&[TT],
              v_nop: &[TP], v_not: &[TT], eqstate:bool, verbose:i8, timetouch:bool) -> i8 {
  let n_v:usize = v_p.len(); let n_n:usize = v_nop.len(); let n_n64:u64 = match u64::try_from(n_n) { Ok(val)=>val,Err(_e)=>0 };
  let n_b_v:usize = find_1(st_v, v_s);
  if n_b_v <= 0 { return 0; } // no buys, don't do anything
  let vstr = format!("parallel_verify_buy({},n_b_v={},n_n={},tt={})",pdiff, n_b_v, n_n, if timetouch {1} else {0});
  if verbose >= 1 { println!("{}: We begin", vstr); }
  if n_v <= 0 {  return -1; } // No buys to match, return nothing. 
  if n_n <= 0 {  return -1; } // No NBBO cant do anything. 
  let st_n_v = sort_nbbo_pt(v_nop,v_not); // Note: NBBO reverse sort is unique to us being given nbb or nbo here
  let st_n = &st_n_v[..];
  if verbose >= 2 {
    println!("{} --- We are beginning with the vectors we have sorted. ", vstr);
    print_v!(format!("{}: st_v,", vstr), st_v,5); print_v!(format!("{}: st_n,", vstr), st_n,5);
  }
  
  // We are scoping this 
  //let mutex_out_v = Arc::new(Mutex::new(vec![n_n64 as usize;n_b_v]));
  let parallelism = thread::available_parallelism().unwrap().get();
  let chunk_size: usize = n_b_v / parallelism;
  let chunk_size: usize = if (chunk_size * parallelism) < n_b_v { chunk_size + 1 } else { chunk_size };
  let tgt_length = chunk_size * parallelism;
  let mut final_out_v:Vec<u64> = vec![n_n64 as u64; tgt_length];  // We know final_out_v is exactly divisible by chunk_size;
  println!("Estimated available parallelism: {}", parallelism); 
  scope(|s| {
    let mut chunks = final_out_v.chunks_exact_mut(chunk_size);
    for (p_i, chunk) in chunks.by_ref().enumerate() { 
      let start_i:usize = chunk_size * p_i;
      let end_i: usize = chunk_size * (p_i+1);
      let vstr_clone = vstr.clone();
      s.spawn(move || {
        let mut on_plev = v_p[st_v[start_i]] + pdiff;
        let mut cross_idx:Vec<usize> = Vec::new();
        let mut on_ix = 0;
        let mut ip_eq:usize = 0; let mut ip_u: usize = 0;
        let mut nx = slow_calculate_cross_above(on_plev, st_n, v_nop, &mut cross_idx, &mut ip_eq, &mut ip_u);
        let mut sort_cross_idx = resort_out_x(&cross_idx);
        if verbose >= 3 {
          println!("{}", format!("{} - thread {}/{} start with on_plev={}, nx={}, ip_eq={}: Prices at cross", vstr_clone, p_i, parallelism, on_plev, nx,ip_eq));
          if sort_cross_idx.len() > 0 {
            print_v!(format!("{} - thread {}/{} start with on_plev={}, nx={}, ip_eq={}: Prices at cross", vstr_clone, p_i, parallelism, on_plev, nx,ip_eq),
                   v_nop, -5 as i64 ,sort_cross_idx,TP,cross_idx);
          }
        }
        //let n_nm1:usize = um1![n_n]; 
        for i_c in 0..chunk.len() {
          let mut pushout: u64 = n_n64;
          let i_v: usize = start_i + i_c as usize;
          if i_v >= n_b_v { break; }
          
          if st_v[i_v] >= v_p.len() {
            println!("Error we have i_v={}, st_v={}, v_p_len={}", i_v, st_v[i_v], v_p.len())
          }
          let v_on_p = v_p[st_v[i_v]]; let v_on_ct = v_ct[st_v[i_v]]; 
          let v_on_ot = v_ot[st_v[i_v]];
          if (v_on_p+pdiff) > on_plev {
            on_plev = v_on_p + pdiff;
            nx = fast_alter_cross_above(on_plev, st_n, v_nop, &mut cross_idx, &mut ip_eq, &mut ip_u);
            sort_cross_idx = resort_out_x(&cross_idx);  // We need to move forward through cross_idx in time.
            on_ix = 0;
            if verbose >= 3 {
              print_v!(format!("{} - thread {}/{} resorted to on_plev={}, nx={}, ip_eq={}: Prices at cross", vstr_clone, p_i, parallelism, on_plev, nx,ip_eq),
                   v_nop, -5 as i64 ,sort_cross_idx,TP,cross_idx);
            }
          }
          if nx == 0 {
            if ip_eq >= n_n {
              if (v_on_ot <= v_not[um1![n_n]]) && (v_on_ct >= v_not[0]) {
                pushout = 0 as u64; } 
              } else if (v_nop[st_n[0]] - pdiff) > on_plev {  pushout  = n_n64; }
          } else {
            for i_x in (on_ix as usize)..(nx as usize) {
              let i_x:usize = i_x as usize;
              let v_on_not = v_not[cross_idx[sort_cross_idx[i_x]]];
              let v_on_not_next = if i_x+1 < nx { v_not[cross_idx[sort_cross_idx[i_x+1]]] } else { v_on_not };
              let current_n = cross_idx[sort_cross_idx[i_x]];
              let next_n = if (current_n+1) < n_n { current_n+1 } else {current_n};
              if v_on_not < v_on_ot {
                if  (i_x +1 < nx) && (v_on_not_next < v_on_ot) { on_ix = on_ix + 1;
                } else if (timetouch==false) && (v_on_not_next == v_on_ot) {
                } else {
                  let lastn:usize = if (i_x+1) < sort_cross_idx.len() { cross_idx[sort_cross_idx[i_x+1]] } else { v_not.len() };
                  for on_n in current_n..lastn {
                    let on_n:usize = on_n as usize; let next_n:usize = if on_n+1 < v_not.len() { on_n+1 } else { on_n };
                    if (v_not[on_n] > v_ct[st_v[i_v]]) || ((timetouch==false) && (v_not[on_n] >= v_on_ct)) {  break;
                    } else if (v_not[next_n] > v_ot[st_v[i_v]]) || ((timetouch==true) && (v_not[next_n] >= v_on_ot)) {
                      if (v_nop[on_n]-pdiff < v_p[st_v[i_v]]) || ((eqstate==true) && (v_nop[on_n]-pdiff <= v_p[st_v[i_v]])) {
                        pushout = on_n as u64; break;
                      }
                    }
                  }
                }
              } else if (v_on_not > v_ct[st_v[i_v]]) || ((timetouch==false) && v_on_not == v_ct[st_v[i_v]]) {
                break;
              } else {
                if (v_nop[current_n] - pdiff < v_on_p) || ((eqstate) && (v_nop[current_n]-pdiff <= v_on_p)) {
                  pushout = current_n as u64;
                } else if (v_not[next_n] < v_on_ct) || ((timetouch) && (v_not[next_n] <= v_on_ct)) {
                  if (v_nop[next_n]-pdiff < v_on_p) || ((eqstate==true) && (v_nop[next_n]-pdiff <= v_on_p)) {
                    pushout = next_n as u64;
                  }
                }
              }
              if pushout < n_n64 { break; }
            }
          }
          if pushout < n_n64 {
              chunk[i_c] = pushout;
          }
        }
      });
    }
  });
  //let final_out_v_data = mutext_out_v.lock().unwrap();
  for i_v in 0..n_b_v {
    out_v[st_v[i_v]] = final_out_v[i_v];
  }
  return 1;
}

// The Slow brute force algorithm
//
// We could potentially GPU this algorithm, though many of the price checks it makes are likely
// unnecessary
pub fn slow_verify_price(out_v: &mut Vec<u64>, pdiff:f64, v_s:&[i8], v_p: &[TP], v_ot:&[TT], v_ct:&[TT],
              v_nop: &[TP], v_not: &[TT], eqstate:bool, verbose:i8, timetouch:bool) -> i8 {
  //let v_nop = *v_nop; let v_not = *v_not;
  let n_v:usize = v_p.len(); let n_n:usize = v_nop.len();
  let n_n64:u64 = n_n as u64;
  let vstr = format!("slow_verify_price(pdiff={},n_v={},n_n={})", 
     pdiff, n_v, n_n);
  if verbose >= 1 {
    println!("{} -- Begin ", vstr);
  }
  //let n_nm1: usize = match usize::try_from( n_n-1 ) { Ok(val)=>val,Err(_e)=>0 };
  //let mut out_v:Vec<u64> =  vec![n_n64;n_v as usize];
  for iv in 0..n_v {
    let iv:usize = iv as usize;
    out_v[iv] = n_n64;
    if v_s[iv] == 0 {
      for ik in 0..n_n {
        let ik:usize = ik as usize;
        if (v_not[ik] > v_ct[iv]) || ((timetouch==false) && (v_not[ik]>=v_ct[iv])) {
          break;
        } else if (ik+1 < n_n) && ( (v_not[ik+1] < v_ot[iv]) || ((timetouch==false) && (v_not[ik+1] <= v_ot[iv]))) {
        } else if (v_p[iv] + pdiff > v_nop[ik]) || ((eqstate==true) && (v_p[iv]+pdiff >= v_nop[ik])) {
          out_v[iv] = ik as u64; break;
        }
      }
    } else {
      for ik in 0..n_n {
        let ik:usize = ik as usize;
        if (v_not[ik] > v_ct[iv]) || ((timetouch==false) && (v_not[ik]>=v_ct[iv])) {
          break;
        } else if (ik+1 < n_n) && ( (v_not[ik+1] < v_ot[iv]) || ((timetouch==false) && (v_not[ik+1] <= v_ot[iv]))) {
        } else if (v_p[iv] -pdiff < v_nop[ik]) || ((eqstate==true) && (v_p[iv]-pdiff <= v_nop[ik])) {
          out_v[iv] = ik as u64; break;
        }
      }
    }     
  }
  return 1;
}

///~///////////////////////////////////////////////////////////////////////////////////
/// verify_sell
///
/// Verifying each sell order does not cross nbb+pdiff (or nbo+pdiff if desired)
///
/// Note this is a reverse/inverse process of verifying buys and this can be confusing.
///
/// We will work in order of Highest price (earliest open time) to Lowest Price (latest open time) on the orders.
///
/// This means we are looking at all NBB crossings of the orders at
///   each fixed "v_p[st_v[i_v]]-pdiff = on_plev": price level
///
///
pub fn verify_sell(out_v: &mut Vec<u64>, st_v: &Vec<usize>, pdiff:f64, v_s:&[i8], v_p: &[TP], v_ot:&[TT], v_ct:&[TT],
          v_nop:&[TP], v_not:&[TT], eqstate:bool, verbose:i8, timetouch:bool) -> i8 {
  //let v_nop = *v_nop; let v_not=*v_not; 
  let n_v:usize = v_p.len(); let n_n:usize = v_nop.len();
  let n_n64:u64 = n_n as u64; //match u64::try_from(n_n) { Ok(val)=>val,Err(_e)=>0 };
  let n_nm1:usize = um1![n_n]; 
  let vstr = format!("verify_sell({},n_v={},n_n={},t_t={})",pdiff, n_v, n_n, if timetouch==false {0} else {1});
  if verbose >= 1 { println!("{} -- We begin.", vstr); }
  if n_n <= 0 {
    // Dumb Case, no NBBO data at all, reject nothing.
    return -1; 
  }

  if n_v <= 0 {
    // Dumb Case, no message data at all, reject nothing.
    return 0; 
  }

  // Problem with Selling, we sort of need to go forward with time, backwards with price.  Thus we
  // find ourselves descending this sort order weirdly.
  //
  // Of course, if I were doing again I guess I'd look for sort solution where ot, vp aren't
  // necessarily in order we see
  //let st_v = sort_v_reverse_pt(&v_s, &v_p, &v_ot); // We get v_p, v_
  let st_n_v = sort_nbbo_pt(v_nop,v_not);
  let st_n = &st_n_v[..];
  // Only Buys return
  if v_s[st_v[match usize::try_from(st_v.len() - 1) { Ok(val)=>val,Err(_e)=>0}]] == 0 { return 1; }
  //let n_vm1 = match usize::try_from( (n_v as i64)-1) { Ok(val)=>val,Err(_e)=>0 };
  let mut on_plev = v_p[st_v[match usize::try_from(st_v.len()-1) { Ok(val)=>val, Err(_e)=>0}]] - pdiff;
  let mut cross_idx:Vec<usize> = Vec::new();
  //let mut on_ix = 0;
   
  let mut ip_eq:usize = um1![n_n]; 

  let mut ip_u: usize = ip_eq;
  let mut nx = slow_calculate_cross_below(on_plev, st_n, v_nop, 
    &mut cross_idx, &mut ip_eq, &mut ip_u);
  let mut sort_cross_idx = resort_out_x(&cross_idx);

  if verbose >= 2 {
     println!("{} - for first price on_plev={}={}[{}=st_v[{}]]-{} we have ip_eq={}/{}, ip_u={}",
       vstr, on_plev, v_p[st_v[0]], st_v[0], 0, pdiff, ip_eq, n_n, ip_u);

     print_v!(format!("{}: sort_cross_idx",vstr), sort_cross_idx, 5);
     print_v!(format!("{}: cross_idx",vstr), cross_idx, 5, sort_cross_idx, usize); 
     print_v!(format!("{}: v_nop[cidx]",vstr), v_nop,5,sort_cross_idx,TP,cross_idx); 
  }
  #[cfg(debug_assertions)]
  {
    // Helpful Debug stuff
    //   Of course, this time to logic-time reverser is useless unless you know a target base time
    //   you want to convert to, rewrite tm_un based upon your content.
    const ORIGCONST:u64 = 1753272000000000000;
    const TMSEC:u64 = 1000000000; 
    fn tm_un(atime: i64)->i64 {
      return (atime-(ORIGCONST as i64))/( TMSEC as i64);
    }
    println!(" --- DEBUG Mode activated tm_un(v_ot[st_v[0]]) = {}", tm_un(v_ot[st_v[0]]));
  }
  for i_v in ((0 as usize)..n_v).rev() {
    let i_v:usize=i_v as usize;  let v_on_ot = v_ot[st_v[i_v]];  let v_on_ct = v_ct[st_v[i_v]]; let v_on_p = v_p[st_v[i_v]];
    if v_s[st_v[i_v]] == 0 { break; }
    #[cfg(debug_assertions)] 
    {
      if verbose >= 3 {
        println!("\n\n----------------------------------------------------------------------------------");
        println!("\n\n Start i_v=({}/{}) or st_v[{}] = {}, plev={}, ot/ct=[{}--{}] ip_eq={}/{}, ip_u={} for v_nop[st_npip_eq={}]]={}, v_nop[st_n[ip_u={}]]={}.", 
         i_v, n_v, i_v, st_v[i_v], v_p[st_v[i_v]], basic_hms(v_ot[st_v[i_v]]), basic_hms(v_ct[st_v[i_v]]),
         ip_eq, n_n, ip_u, 
        ip_eq, v_nop[st_n[ip_eq]], ip_u, v_nop[st_n[ip_u]]);
        println!("     previous on_plev ={}, but new on_plev will be {}. ", on_plev, v_p[st_v[i_v]]-pdiff);
      }
      // Debug Mode chekcks
      if (v_p[st_v[i_v]] == 103.0) && (v_ot[st_v[i_v]] == (175327200000000000 + 7*1000000000 as i64)) {
        println!("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
        println!("     here we are with v_on_p={}, v_ot on target., i_v={}/{} st_v[i_v] = {}", v_p[st_v[i_v]], i_v, st_v.len(), st_v[i_v] );
        println!("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
      } else if v_p[st_v[i_v]] == 103.0 {
         println!("YYY NOT here, we are 103, but v_ot[{}] = {} ", st_v[i_v], basic_hms(v_ot[st_v[i_v]]));
      }
    }
    
    if (v_on_p-pdiff) < on_plev {
      let old_ip_eq=ip_eq;  let old_plev = on_plev;
      on_plev = v_p[st_v[i_v]] - pdiff;  
      // Note you can use update_cross_below, but alter_cross_below should be faster
      nx = fast_alter_cross_below(on_plev, st_n, v_nop, 
        &mut cross_idx, &mut ip_eq, &mut ip_u); //, if i_v==3 {10} else {0});
      sort_cross_idx = resort_out_x(&cross_idx);  // We need to move forward through cross_idx in time.
      //on_ix = 0;
      
      #[cfg(debug_assertions)]
      {
       if verbose >= 2 {
        println!("{} - for renewed price on_plev={}={}[{}=st_v[{}]]-{} we have ip_eq={}/{}, ip_u={}/{}: range v_nop=[{},{}]",
          vstr, on_plev, v_p[st_v[i_v]], st_v[i_v], i_v, pdiff, ip_eq, n_n, ip_u, n_n, v_nop[st_n[0]], v_nop[st_n[n_nm1]]);
        print_v!(format!("{}: sort_cross_idx",vstr), sort_cross_idx, 5);
        print_v!(format!("{}: cross_idx",vstr), cross_idx, 5, sort_cross_idx, usize); 
        print_v!(format!("{}: v_nop[cidx]",vstr), v_nop,5,sort_cross_idx,TP,cross_idx); 
        let nerr = check_jumps(&format!("verify_sell[i_v={}]",i_v), on_plev, v_nop, &sort_cross_idx, &cross_idx);
        let nerr2 = full_check_jumps(&format!("verify_sell[i_v={}]",i_v), on_plev, v_nop, &sort_cross_idx, &cross_idx);
        println!("{} -- We checked jumps and got {} errors for {} finds, on_plev={} n_n={}, nerr2={}", 
          vstr, nerr, nx, on_plev, n_n, nerr2);  if nerr > 0 { println!("{} ERROR we know we aren't doing good.", vstr); }
        if nerr2 > 0 {
          println!("What is wrong with our errors?");
          println!("   -- on_plev={}, i_v={}/{} for st_v[{}]={} and v_p[{}]={}, v_ot[{}]={}, v_ct[{}]={} ",
            on_plev, i_v, n_v, i_v, st_v[i_v], st_v[i_v], v_p[st_v[i_v]], st_v[i_v], v_ot[st_v[i_v]], st_v[i_v], v_ct[st_v[i_v]]);
          println!("   Old ip_eq={} for old_ip_eq={} for v_nop[st_n[{}]={}]=={} and v_nop[st_n[{}]={}] = {}, on_plev={}, old_plev={} ",
             ip_eq, old_ip_eq, ip_eq, st_n[ip_eq], v_nop[st_n[ip_eq]], old_ip_eq, st_n[old_ip_eq], v_nop[st_n[old_ip_eq]], on_plev, old_plev);
          println!(" --- Printing the st_n and v_nop data");
          let rev_idxs = (0..n_n).rev().collect::<Vec<usize>>();
          print_v!(" --- To start rev_idxs is ", rev_idxs, 200);
          let st_n_above_eq = rev_idxs.iter().filter_map(|&idx| if idx > ip_eq { Some(st_n[idx]) } else {None}).rev().collect::<Vec<usize>>();
          print_v!(format!("  reverse order st_n[n_n={}..>ip_eq={}] is: ", n_n, ip_eq), st_n_above_eq, 26);
          print_v!(format!("  Prices[n_n..>ip_eq={}] is: ", ip_eq), v_nop, 26, st_n_above_eq, TP); 
          //print_v!(format!("  Prices[n_n..>ip_eq={}] is: ", ip_eq), (st_n_above_eq.iter().filter_map(|&idx| Some(v_nop[idx])).collect::<Vec<f64>>()),26);
          println!(" --- Note ip_eq ={}(st_n[ip_eq]={}]) and ip_u={}(st_n[ip_eq={}) for prices {} and {}", ip_eq, st_n[ip_eq], ip_u, st_n[ip_u],v_nop[st_n[ip_eq]], v_nop[st_n[ip_u]]);
          let st_n_eq_to_u = rev_idxs.iter().filter_map(|&idx| if (idx <= ip_eq) && (idx > ip_u) { Some(st_n[idx]) } else {None}).collect::<Vec<usize>>(); 
          print_v!(format!("  reverse order st_n[ip_eq={}..>ip_u={}] is: ", ip_eq, ip_u), st_n_eq_to_u, 26); 
          print_v!(format!("  Prices[ip_eq={}..>ip_u={}] is: ", ip_eq, ip_u), v_nop, 26, st_n_eq_to_u, TP);
          //print_v!(format!("  Prices[ip_eq={}..>ip_u={}] is: ", ip_eq, ip_u), (st_n_eq_to_u.iter().filter_map(|&idx| v_nop.get(idx).cloned()).collect::<Vec<f64>>()),26);
          println!(" ---- And then what happens after ip_u={} or Price={}", ip_u, v_nop[st_n[ip_u]]);
          let st_n_from_u = rev_idxs.iter().filter_map(|&idx| if idx <= ip_u { Some(st_n[idx]) } else {None}).collect::<Vec<usize>>(); 
          print_v!(format!("  reverse order st_n[ip_u={}..>=0] is: ", ip_u), st_n_from_u, 26);
          print_v!(format!("  Prices[ip_u={}..>=0] is: ", ip_u), v_nop, 26, st_n_from_u, TP); 
          //print_v!(format!("  Prices[ip_u={}..>=0] is: ", ip_u), (st_n_from_u.iter().filter_map(|&idx| v_nop.get(idx).cloned()).collect::<Vec<f64>>()),26);
        }
       }
      }
    }
    if (verbose >= 2) ||  ((verbose >= 1) && (i_v % 1000 == 0)) {
      println!("   i_v={}/{} for st_v[{}]={}. ", i_v, n_v, i_v, st_v[i_v]);
      print_v!(format!("{}: sort_cross_idx",vstr), sort_cross_idx, 5);
      print_v!(format!("{}: cross_idx",vstr), cross_idx, 5, sort_cross_idx, usize); 
      print_v!(format!("{}: v_nop[cidx]",vstr), v_nop,5,sort_cross_idx,TP,cross_idx); 
    }
    if nx == 0 {
      if ip_eq <= 0 {
        // In this case all prices are always more if there is any overlap.  Reject
        if (v_on_ot <= v_not[n_n-1]) && (v_on_ct >= v_not[0]) {
          out_v[st_v[i_v]] = 0 as u64;
        } else {
          // Crazy case everyone out of range
          out_v[st_v[i_v]] = n_n64;
        }
      } else if (v_nop[st_n[0]] + pdiff) > on_plev {
        // In this case all prices are north of issue.  We should keep this point.
        out_v[st_v[i_v]] = n_n64;
      }
    } else {
      for i_x in 0..(nx as usize) {
         let i_x:usize = i_x as usize;
         let v_on_not = v_not[cross_idx[sort_cross_idx[i_x]]];
         let current_n:usize = cross_idx[sort_cross_idx[i_x]]; let next_n:usize = if (current_n+1) < n_n { current_n+1} else {current_n}; 
         //let next_n:usize = if (on_n+1) &&
         let v_on_not_next = if i_x + 1 < nx { v_not[cross_idx[sort_cross_idx[i_x+1]]] } else { v_on_not };
         if v_on_not < v_on_ot {
            if (i_x+1 < nx) && (v_on_not_next < v_on_ot) {
              //on_ix = on_ix + 1;
            } else if (timetouch==false) && (v_on_not_next == v_on_ot) {
            } else {
              // If we are at last break before v_ot we have to as questions.
              let lastn:usize = if (i_x+1) < sort_cross_idx.len() { cross_idx[sort_cross_idx[i_x+1]] } else { v_not.len() };
              for on_n in current_n..lastn {
                let on_n:usize = on_n as usize; let next_n:usize = if on_n+1 < v_not.len() { on_n+1 } else { on_n };
                if (v_not[on_n] > v_ct[st_v[i_v]]) || ((timetouch==false) && (v_not[on_n] >= v_on_ct)) { 
                  break;
                } else if (v_not[next_n] > v_ot[st_v[i_v]]) || ((timetouch==true) && (v_not[next_n] >= v_on_ot)) {
                  if (v_nop[on_n] > on_plev) || ((eqstate==true) && (v_nop[on_n] >= on_plev)) {
                    out_v[st_v[i_v]] = on_n as u64; break;
                  }
                }
              }
            }
         } else if (v_on_not > v_ct[st_v[i_v]]) || ((timetouch==false) && v_on_not == v_ct[st_v[i_v]]) {
           //println!("   Ending i_v={}/{} [{}] on i_x={}/{} with out_v[st_v[{}]] = {}, n_n64={}.",
           //  i_v, st_v.len(), st_v[i_v], i_x,nx, st_v[i_v], out_v[st_v[i_v]], n_n64);
           break;
         } else {
           if (v_nop[current_n] > on_plev) || ((eqstate==true) && (v_nop[current_n] >= on_plev)) {
             out_v[st_v[i_v]] = current_n as u64;
           } else if (v_not[next_n] < v_on_ct) || ((timetouch) && (v_not[next_n] <= v_on_ct)) {
             if (v_nop[next_n] > on_plev) || ((eqstate == true) && (v_nop[next_n] >= on_plev)) {
               out_v[st_v[i_v]] = next_n as u64;
             }
           }
         }
         if out_v[st_v[i_v]] < n_n64 { break; }
       }
       if out_v[st_v[i_v]] >= n_n64 {
         out_v[st_v[i_v]] = n_n64;  // We have passed all challenges.
       }
    }
  }

  return 1;
}

///~//// parralel version of sell verification algorithms
pub fn parallel_verify_sell(out_v: &mut Vec<u64>, st_v: &Vec<usize>, pdiff:f64, v_s:&[i8], v_p: &[TP], v_ot:&[TT], v_ct:&[TT],
          v_nop:&[TP], v_not:&[TT], eqstate:bool, verbose:i8, timetouch:bool) -> i8 {
  //let v_nop = *v_nop; let v_not=*v_not; 
  let n_v:usize = v_p.len(); let n_n:usize = v_nop.len();
  let n_n64:u64 = n_n as u64; //match u64::try_from(n_n) { Ok(val)=>val,Err(_e)=>0 };
  //let n_nm1:usize = um1![n_n]; 
  let n_b_v:usize = find_1(st_v, v_s);
  if n_v - n_b_v <= 0 { return(0); }
  let vstr = format!("parallel_verify_sell({},n_v={},n_n={},t_t={})",pdiff, n_v, n_n, if timetouch==false {0} else {1});
  if verbose >= 1 { println!("{} -- We begin.", vstr); }
  if n_n <= 0 {
    // Dumb Case, no NBBO data at all, reject nothing.
    return -1; 
  }

  if n_v <= 0 {
    // Dumb Case, no message data at all, reject nothing.
    return 0; 
  }

  // Problem with Selling, we sort of need to go forward with time, backwards with price.  Thus we
  // find ourselves descending this sort order weirdly.
  //
  // Of course, if I were doing again I guess I'd look for sort solution where ot, vp aren't
  // necessarily in order we see
  //let st_v = sort_v_reverse_pt(&v_s, &v_p, &v_ot); // We get v_p, v_
  let st_n_v = sort_nbbo_pt(v_nop,v_not);
  let st_n = &st_n_v[..];
  // Only Buys return
  if v_s[st_v[match usize::try_from(st_v.len() - 1) { Ok(val)=>val,Err(_e)=>0}]] == 0 { return 1; }

  // We are scoping this 
  //let mutex_out_v = Arc::new(Mutex::new(vec![n_n64 as usize;n_b_v]));
  let parallelism = thread::available_parallelism().unwrap().get();
  let n_s_v:usize = match usize::try_from ( (n_v as i64) - (n_b_v as i64)) { Ok(val)=>val,Err(_e)=>0 };
  let chunk_size: usize = (n_s_v)  / parallelism;
  let chunk_size: usize = if (chunk_size * parallelism) < n_s_v { chunk_size + 1 } else { chunk_size };
  let tgt_length = chunk_size * parallelism;
  let mut final_out_v:Vec<u64> = vec![n_n64 as u64; tgt_length];  // We know final_out_v is exactly divisible by chunk_size;
  println!("Estimated available parallelism: {}", parallelism); 
  scope(|s| {
    let mut chunks = final_out_v.chunks_exact_mut(chunk_size);
    for (p_i, chunk) in chunks.by_ref().enumerate() { 
      let start_i:usize = n_b_v + chunk_size * p_i;  //match usize::try_from((n_v as i64) - (chunk_size * p_i as i64)) { Ok(val)=>val,Err(_e)=>0 };
      let end_i: usize =  if start_i + chunk_size  < n_v { start_i + chunk_size } else {n_v};
      let n_c_i: usize = match usize::try_from((end_i as i64)-(start_i as i64)) { Ok(val)=>val,Err(_e)=>0 };
      let vstr_clone = vstr.clone();
      s.spawn(move || {
        if end_i > n_v { return; }
        let mut on_plev = v_p[st_v[um1![end_i]]] + pdiff;
        let mut cross_idx:Vec<usize> = Vec::new();
        let mut ip_eq:usize = um1![n_n]; let mut ip_u: usize = ip_eq;
        let mut nx = slow_calculate_cross_below(on_plev, st_n, v_nop, &mut cross_idx, &mut ip_eq, &mut ip_u);
        let mut sort_cross_idx = resort_out_x(&cross_idx);
        if verbose >= 2 {
          print_v!(format!("{} - thread={}/{} start={} for first price on_plev={}={}[{}=st_v[{}]]-{} we have ip_eq={}/{}, ip_u={}",
            vstr_clone, p_i, parallelism, start_i, on_plev, v_p[st_v[0]], st_v[0], 0, pdiff, ip_eq, n_n, ip_u), v_nop, -5 as i64, sort_cross_idx, TP, cross_idx);
        }
        for i_c in (0..n_c_i).rev() {
          let i_v:usize=i_c + start_i as usize;  let v_on_ot = v_ot[st_v[i_v]];  let v_on_ct = v_ct[st_v[i_v]]; let v_on_p = v_p[st_v[i_v]];
          let mut pushback: u64 = n_n64;
          if v_s[st_v[i_v]] == 0 { break; }
          if (v_on_p-pdiff) < on_plev {
            on_plev = v_p[st_v[i_v]] - pdiff;  
            nx = fast_alter_cross_below(on_plev, st_n, v_nop, &mut cross_idx, &mut ip_eq, &mut ip_u); //, if i_v==3 {10} else {0});
            sort_cross_idx = resort_out_x(&cross_idx);  // We need to move forward through cross_idx in time.
            if verbose >= 2 {
              print_v!(format!("{} - thread={}/{} start= {}: for alteration to  price on_plev={}={}[{}=st_v[{}]]-{} we have ip_eq={}/{}, ip_u={}",
              vstr_clone, p_i, parallelism, start_i, on_plev, v_p[st_v[0]], st_v[0], 0, pdiff, ip_eq, n_n, ip_u), v_nop, -5 as i64, sort_cross_idx, TP, cross_idx);
            }
          }
      
          if nx == 0 {
            if (ip_eq <= 0) && (v_on_ot <= v_not[n_n-1]) && (v_on_ct >= v_not[0]) { pushback = 0 as u64;
            } 
          } else {
            for i_x in 0..(nx as usize) {
              let i_x:usize = i_x as usize;
              let v_on_not = v_not[cross_idx[sort_cross_idx[i_x]]];
              let current_n:usize = cross_idx[sort_cross_idx[i_x]]; let next_n:usize = if (current_n+1) < n_n { current_n+1} else {current_n}; 
              //let next_n:usize = if (on_n+1) &&
              let v_on_not_next = if i_x + 1 < nx { v_not[cross_idx[sort_cross_idx[i_x+1]]] } else { v_on_not };
              if v_on_not < v_on_ot {
                if (i_x+1 < nx) && (v_on_not_next < v_on_ot) {
                } else if (timetouch==false) && (v_on_not_next == v_on_ot) {
                } else {
                  // If we are at last break before v_ot we have to as questions.
                  let lastn:usize = if (i_x+1) < sort_cross_idx.len() { cross_idx[sort_cross_idx[i_x+1]] } else { v_not.len() };
                  for on_n in current_n..lastn {
                    let on_n:usize = on_n as usize; let next_n:usize = if on_n+1 < v_not.len() { on_n+1 } else { on_n };
                    if (v_not[on_n] > v_ct[st_v[i_v]]) || ((timetouch==false) && (v_not[on_n] >= v_on_ct)) { 
                      break;
                    } else if (v_not[next_n] > v_ot[st_v[i_v]]) || ((timetouch==true) && (v_not[next_n] >= v_on_ot)) {
                      if (v_nop[on_n] > on_plev) || ((eqstate==true) && (v_nop[on_n] >= on_plev)) {
                        pushback = on_n as u64; break;
                      }
                    }
                  }
                }
              } else if (v_on_not > v_ct[st_v[i_v]]) || ((timetouch==false) && v_on_not == v_ct[st_v[i_v]]) {
                break;
              } else {
                if (v_nop[current_n] > on_plev) || ((eqstate==true) && (v_nop[current_n] >= on_plev)) {
                  pushback = current_n as u64;
                } else if (v_not[next_n] < v_on_ct) || ((timetouch) && (v_not[next_n] <= v_on_ct)) {
                  if (v_nop[next_n] > on_plev) || ((eqstate == true) && (v_nop[next_n] >= on_plev)) {
                    pushback = next_n as u64;
                  }
                } 
              }
              if pushback < n_n64 { break; }
            }
          }
          if pushback <= n_n64 {
            chunk[i_c] = pushback;
          }
        }
      });  
    }
  });  
  for i_z in 0..n_s_v {
     out_v[st_v[n_b_v + i_z]] = final_out_v[i_z];
  }
  return 1;
}
