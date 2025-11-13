///~//////////////////////////////////////////////////////////////////////
/// sip_algo.rs
///
///  Alan Lenarcic, November 2024
///
/// Algorithms for calculating the National Best Bid/Offer from a running
///  feed of uqdf/cqs System Information Processors (SIP) feed are ubiquitous
///  and a basic challenge.
///
/// Here we adapt a typical direct C-based strategy using Rust idioms.
///
/// We keep a running vector of the best price/best offer of the n_venues
/// 
/// The current state is stored in a TSip object ("tsip" in this case).
/// The state is updated with the next record (record index "i_r" in 0..n_r
/// Each record updates the best bid and/or best ask from one venue.
///
/// We merely need to identify if this venue is joining the national best
///   bid or offer, or if it sets a new better price, if the venue is 
///   leaving the nbbo or, if the change is merely at later venues.
///
/// We can store a simple binary structure which records exactly which venues
///   stand at the maximum in each state.  (this is tsip.b.b2v or tsip.s.b2v)
///   That way, we are not spending much time
///
/// We record each step in RecSip/VPublishSip structures. We only need one
///   but it is useful to conser column-based and row-based published records
/// This affects whether we can export the array in an Arrow based structure
///   or a numpy--pandas based structure
///
/// Numpy (rec arrays that map to pandas) can be a more common tabular 
///   structure returned by algorithms.  However, there are not easy programatic
///   safe ways to generate an arbitrary shaped numpy table in rust.
///   (This is hard to do with C-Api too, its easier to declare these in python)
///
/// Arrow is a more unsual tabular structure to emit results in, but these results
///   can be saved immediately to disk as .parquet, so they are potentially
///   platform agnostic (python/R/lua all accept these. 
///
/// It should be possible to very quickly evaluate multiple millions of these records
///   in a second.
///
/// The only extra challenge to this algorithm is that it is difficult to predict
///   how many "useful" changes occure.  Any time there a is a record that doesn't
///   change the NBBO meaningfully (it depends what the user thinks is meaningful)
///   or if two records share the same timestamp and happen "simultaneously" and
///   we don't want to record an artificial intermediary state.  
use crate::sip_struct::{TSip, RecSip, VPublishSip};
//use crate::sip_struct::{PublishSipI};
use crate::sip_struct::{TNV, TT, TP, TQ, ILine, dt64_string, u64_string};

const DEBUG_ALGO:u8 = 0;

////////////////////////////////////////////////////////////////////////////////
// Compute_SIP_NBBO(...)
//
//    A Rust demonstration of common SIP NBBO algo
//    For two Venue sequences, for bid (buy) and ask (sell), 
//     track running maximum and which venues sit in the place.
//
//  Requires two Memory structures:
//    1. a TSIP (contains two pSIPS, one for bid and ask each)
//    2. a recSIP (records the running changes [often less than current changes]
//
//  Inputs:
//    nVenues: Number of Venues (will call it "nv" in TSIP memory structure)
//    vT: a vector of Time, sorted in increasing sequence
//    vBP, vBQ, vBV: the Bid vectors of Price, Quantity, Venue from SIP records
//    vSP, vSQ, vSV: the Ask vectors for Price, Quantity, Venue from SIP records
//
//  Output: Ideally a recSIP, that can be exported later to R/Arrow/Python with
//     the dangling memory removed. It is difficult to determine how many
//     records are actually required until algorithm is finished.
pub fn compute_sip_nbbo(n_venues: u8, 
  vt_in:Vec<TT>, vb_q_in:Vec<TQ>, vs_q_in:Vec<TQ>, 
  vb_p_in:Vec<TP>, vs_p_in: Vec<TP>,
  vb_v_in:Vec<TNV>, vs_v_in:Vec<TNV>, verbose: u8) -> Option<VPublishSip> { 
  let n_r: u32 = vt_in.len().try_into().unwrap();
  let n_ur = n_r as usize;
  if (n_ur != vb_q_in.len()) || (n_ur != vb_p_in.len()) || (n_ur != vb_v_in.len()) ||
     (n_ur != vs_q_in.len()) || (n_ur != vs_p_in.len()) || (n_ur != vs_v_in.len()) {
    println!("sip_algo.rs->ComputeSIP -- lens:[vt={},vb_p={},\
      vb_q={},vb_v={},vs_p={},vs_q={},vs_v={}]. -- ERROR should return",
      vt_in.len(), vb_p_in.len(), vb_q_in.len(), vb_v_in.len(), 
      vs_p_in.len(), vs_q_in.len(), vs_v_in.len()); 
    return None ;
  }
  if verbose > 0 {
    println!("sip_algo.rs->Compute_SIP_NBBO(VB={}, n_venues={}, n_r={}) -- initiate() -> ",
      verbose, n_venues, n_r);
  } 
  let mut tsip: TSip = TSip::new(n_venues as TNV);
  let mut rec_sip: RecSip = RecSip::new(n_r, n_venues);
  let mut vpublish_sip: VPublishSip = VPublishSip::new(n_r);
  if verbose >= 1 {
    println!("sip_algo.rs->ComputeSIP(verbose={},nv={},input_len={})",
      verbose, n_venues, vt_in.len());
  }
  tsip.b.b2v.zero_all();  tsip.s.b2v.zero_all();
  let mut on_r = 0;
  for i_r in 0..n_r {
    if (verbose >= 2) && (i_r % 1000) == 0 {
      println!("sip_algo.rs({}/{}) -- {}: moving with B(q,p,v)=({},{},{}),S(q,p,v)=({},{},{}) ",
        i_r, n_r, dt64_string(vt_in[i_r as usize]), vb_q_in[i_r as usize], 
		vb_p_in[i_r as usize], vb_v_in[i_r as usize], 
        vs_q_in[i_r as usize], vs_p_in[i_r as usize], vs_v_in[i_r as usize]);
    }
    // Pretty Simple, update Bid.  Update Ask
    //  Then Verify quality, and record.
    let _last_b_p = tsip.b.p;  let _last_s_p = tsip.s.p;
    tsip.tt = vt_in[i_r as usize];  tsip.i_r = i_r as usize;
    tsip.b.update(vb_v_in[i_r as usize], vb_p_in[i_r as usize], vb_q_in[i_r as usize]);
    tsip.s.update(vs_v_in[i_r as usize], vs_p_in[i_r as usize], vs_q_in[i_r as usize]);
    if  DEBUG_ALGO > 0 {
      if tsip.verify() > 0 {
        println!("Fail on an update i_r = {}/{}", i_r, n_r);
        println!("We had vb_vpq_in[{}=[{},{},{}]", i_r, vb_v_in[i_r as usize], vb_p_in[i_r as usize], vb_q_in[i_r as usize]);
        println!("We had vs_vpq_in[{}=[{},{},{}]", i_r, vs_v_in[i_r as usize], vs_p_in[i_r as usize], vs_q_in[i_r as usize]);
        println!(" Last b/s prices were {},{}.", _last_b_p, _last_s_p);
        return None;
      }
    }
    // Note record_state is vector record, we don't really need but here for practice
    rec_sip.record_state(&tsip);
    if (i_r >= n_r-1) || (vt_in[i_r as usize] != vt_in[(i_r+1) as usize]) {
      // Deciding when to print, we want to only print on new time stamps
      // We also only want to print if NBBO has noticably changed somehow
      if on_r == 0 {
        // Print because we never printed.
        vpublish_sip.publish(vt_in[i_r as usize], i_r as ILine, &tsip); on_r +=1;
      } else if (tsip.b.p != vpublish_sip.v[vpublish_sip.on_r-1].b_p) || 
         (tsip.s.p != vpublish_sip.v[vpublish_sip.on_r-1].s_p) ||
         (tsip.b.q != vpublish_sip.v[vpublish_sip.on_r-1].b_q) || 
         (tsip.s.q != vpublish_sip.v[vpublish_sip.on_r-1].s_q) ||
         (tsip.b.b2v.out64() != vpublish_sip.v[vpublish_sip.on_r-1].b_vu) || 
         (tsip.s.b2v.out64() != vpublish_sip.v[vpublish_sip.on_r-1].s_vu) {
        // Only need to print noticed differences
        vpublish_sip.publish(vt_in[i_r as usize], i_r as ILine, &tsip); on_r +=1;
      }
    }
    if verbose >= 2 && ((i_r % 1000 == 0) || (i_r < 10)) {
      println!("i_r={}/{} - [t={},b[p,q,n,b2]=[{},{},{},{}],s[p,q,n,b2]=[{},{},{},{}]]",
         i_r, n_r, dt64_string(tsip.tt), tsip.b.p,tsip.b.q,
         tsip.b.nwv,tsip.b.b2v.p_string(n_venues as u16),
         tsip.s.p,tsip.s.q,tsip.s.nwv,tsip.s.b2v.p_string(n_venues as u16));
      if on_r > 0 {
        let o_r:usize = (on_r -1).try_into().unwrap();
          println!("PUBLISH -- [t={},b[p,q,n,b2]=[{},{},{},{}={}]", 
            dt64_string(vpublish_sip.v[o_r].tt), vpublish_sip.v[o_r].b_p, 
            vpublish_sip.v[o_r].b_q, vpublish_sip.v[o_r].b_nwv,
            vpublish_sip.v[o_r].b_vu,
            u64_string(vpublish_sip.v[o_r].b_vu,n_venues as u16));
          println!("  --------- s[p,q,n,b2]=[{},{},{},{}={}]]",
            vpublish_sip.v[o_r].s_p, vpublish_sip.v[o_r].s_q, 
            vpublish_sip.v[o_r].s_nwv,vpublish_sip.v[o_r].s_vu,
            u64_string(vpublish_sip.v[o_r].s_vu, n_venues as u16));
      }
    }
  }
  return Some(vpublish_sip);
}

// Same code, now it returns the rec_sip instead, which is better column supported	
pub fn compute_sip_nbbo_rec_sip(n_venues: u8, 
  vt_in:Vec<TT>, vb_q_in:Vec<TQ>, vs_q_in:Vec<TQ>, 
  vb_p_in:Vec<TP>, vs_p_in: Vec<TP>,
  vb_v_in:Vec<TNV>, vs_v_in:Vec<TNV>, verbose: u8) -> Option<RecSip> { 
  let n_r: u32 = vt_in.len().try_into().unwrap();
  let n_ur = n_r as usize;
  if (n_ur != vb_q_in.len()) || (n_ur != vb_p_in.len()) || (n_ur != vb_v_in.len()) ||
     (n_ur != vs_q_in.len()) || (n_ur != vs_p_in.len()) || (n_ur != vs_v_in.len()) {
    println!("sip_algo.rs->compute_sip_nbbo_rec_sip -- lens:[vt={},vb_p={},\
      vb_q={},vb_v={},vs_p={},vs_q={},vs_v={}]. -- ERROR should return",
      vt_in.len(), vb_p_in.len(), vb_q_in.len(), vb_v_in.len(), 
      vs_p_in.len(), vs_q_in.len(), vs_v_in.len()); 
    return None ;
  }
  if verbose > 0 {
    println!("sip_algo.rs->compute_sip_nbbo_rec_sip(VB={}, n_venues={}, n_r={}) -- initiate() -> ",
      verbose, n_venues, n_r);
  } 
  let mut tsip: TSip = TSip::new(n_venues as TNV);
  let mut rec_sip: RecSip = RecSip::new(n_r, n_venues);
  let mut vpublish_sip: VPublishSip = VPublishSip::new(n_r);
  if verbose >= 1 {
    println!("sip_algo.rs->ComputeSIP(verbose={},nv={},input_len={})",
      verbose, n_venues, vt_in.len());
  }
  tsip.b.b2v.zero_all();  tsip.s.b2v.zero_all();
  let mut on_r = 0;
  for i_r in 0..n_r {
    if (verbose >= 2) && (i_r % 1000) == 0 {
      println!("sip_algo.rs({}/{}) -- {}: moving with B(q,p,v)=({},{},{}),S(q,p,v)=({},{},{}) ",
        i_r, n_r, dt64_string(vt_in[i_r as usize]), vb_q_in[i_r as usize], 
		vb_p_in[i_r as usize], vb_v_in[i_r as usize], 
        vs_q_in[i_r as usize], vs_p_in[i_r as usize], vs_v_in[i_r as usize]);
    }
    // Pretty Simple, update Bid.  Update Ask
    //  Then Verify quality, and record.
    tsip.tt = vt_in[i_r as usize];  tsip.i_r = i_r as usize;
    tsip.b.update(vb_v_in[i_r as usize], vb_p_in[i_r as usize], vb_q_in[i_r as usize]);
    tsip.s.update(vs_v_in[i_r as usize], vs_p_in[i_r as usize], vs_q_in[i_r as usize]);
    if  DEBUG_ALGO > 0 {
      tsip.verify();
    }
    // Note record_state is vector record, which seems easier
    // to work with for returning a Arrow based table
    rec_sip.record_state(&tsip);
    if (i_r >= n_r-1) || (vt_in[i_r as usize] != vt_in[(i_r+1) as usize]) {
      // Deciding when to print, we want to only print on new time stamps
      // We also only want to print if NBBO has noticably changed somehow
      if on_r == 0 {
        // Print because we never printed.
        vpublish_sip.publish(vt_in[i_r as usize], i_r as ILine, &tsip); on_r +=1;
      } else if (tsip.b.p != vpublish_sip.v[vpublish_sip.on_r-1].b_p) || 
         (tsip.s.p != vpublish_sip.v[vpublish_sip.on_r-1].s_p) ||
         (tsip.b.q != vpublish_sip.v[vpublish_sip.on_r-1].b_q) || 
         (tsip.s.q != vpublish_sip.v[vpublish_sip.on_r-1].s_q) ||
         (tsip.b.b2v.out64() != vpublish_sip.v[vpublish_sip.on_r-1].b_vu) || 
         (tsip.s.b2v.out64() != vpublish_sip.v[vpublish_sip.on_r-1].s_vu) {
        // Only need to print noticed differences
        vpublish_sip.publish(vt_in[i_r as usize], i_r as ILine, &tsip); on_r +=1;
      }
    }
    if verbose >= 2 && ((i_r % 1000 == 0) || (i_r < 10)) {
      println!("i_r={}/{} - [t={},b[p,q,n,b2]=[{},{},{},{}],s[p,q,n,b2]=[{},{},{},{}]]",
         i_r, n_r, dt64_string(tsip.tt), tsip.b.p,tsip.b.q,
         tsip.b.nwv,tsip.b.b2v.p_string(n_venues as u16),
         tsip.s.p,tsip.s.q,tsip.s.nwv,tsip.s.b2v.p_string(n_venues as u16));
      if on_r > 0 {
        let o_r:usize = (on_r -1).try_into().unwrap();
          println!("PUBLISH -- [t={},b[p,q,n,b2]=[{},{},{},{}={}]", 
            dt64_string(vpublish_sip.v[o_r].tt), vpublish_sip.v[o_r].b_p, 
            vpublish_sip.v[o_r].b_q, vpublish_sip.v[o_r].b_nwv,
            vpublish_sip.v[o_r].b_vu,
            u64_string(vpublish_sip.v[o_r].b_vu,n_venues as u16));
          println!("  --------- s[p,q,n,b2]=[{},{},{},{}={}]]",
            vpublish_sip.v[o_r].s_p, vpublish_sip.v[o_r].s_q, 
            vpublish_sip.v[o_r].s_nwv,vpublish_sip.v[o_r].s_vu,
            u64_string(vpublish_sip.v[o_r].s_vu, n_venues as u16));
      }
    }
  }
  if verbose >= 2 {
    println!("sip_algo.rs->compute_sip_rec_sip -- Concluded with on_r={},n_r={}",
      on_r, n_r);
  }
  return Some(rec_sip);
}
