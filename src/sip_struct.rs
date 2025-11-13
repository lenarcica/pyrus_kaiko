//////////////////////////////////////////////////////////////////
// Structures for reproducing market sip
//
//  Alan B Lenarcic 2024/11/09
//
//  Learning some efficiencies in SIP structures in a 
//    RUST environment.  
//
//  We will borrow some finance KDB+ style to illustrate using simple
//    short as possible names for common elements in a type that
//    features often in code will be helpful
//
//    Trying to use condensed KDB style term names

use crate::b2v_struct::{B2Vs};
//use create::b2v_struct::{B2V, NMAXVENS, NMAXBITVENS};
use numpy::pyo3::{Python};
use numpy::pyo3::{Bound,types::{PyAnyMethods} };
use numpy::PyArrayDescr;
use numpy::Element;
use numpy::get_array_module;
use pyo3::types::IntoPyDict;
use std::fmt;
use std::vec::Vec;
use sprintf::sprintf;
//extern crate derive_more;

// Attempt to not have to write (*self) over again
// Alas, Macros in Rust do not permit taking over this convention.
// Instead Rust has a sort of quasi "(*self)."=="self." convention
// It is called, "automatic referencing and dereferencing"
//  But frankly it's properties are confusing.
//macro_rules!ME{()=>{(*self)};}


pub type TNV = u8; //Note this number will be different if NMAXVENS goes more than 128 (probably 20)
pub type TP = f64;
pub type TQ = f64;
//pub type TPi = u64;
pub type TT = i64;
pub type ILine = u64;

#[derive(Debug,  PartialEq, Eq)]
pub enum Emside {
  Buy, 
  Sell
}

// Represent datetime64 integer as string.
pub fn dt64_string(x:TT) -> String  {
  //let ss = String::new();
  const ML: [i64;12] = [31,28,31,30,31,30,31,31,30,31,30,31];
  if x > 0 {
    let ns: i64 = x % 1000000000;
    let sec: i64 = (x / 1000000000) % 60;
    let min: i64 = (x / 60000000000)  % 60;
    let hr: i64 = (x / 3600000000000) % 24;
    let days: i64 = x / (24*3600000000000);
    let mut yr = 1970; let mut tdays = days; let mut found=0; 
    while found == 0 {
      let pump = if (yr%4)==0 {366} else {365 };
	 if tdays < (pump) {
        found = 1;
	 } else {
        tdays = tdays -pump; yr = yr + 1;
	 }
    }
    found = 0;  let mut mth = 0;
    while found == 0 {
      let pump = if mth==1 { if yr%4==0 { 29 } else { 28 } } else { ML[mth] };
      if tdays < pump {
        found = 1;
      } else {
        tdays = tdays - pump; mth = mth +1;
      }
    } 
    return sprintf!("%04ld-%02ld-%02ld %02ld:%02ld:%02ld.%09ld", yr, mth, tdays, hr, min, sec, ns).unwrap();
  }
  let ns: i64 = 1000000000 + (x % 1000000000);
  let sec: i64 = 60 + (x / 1000000000) %60;
  let min: i64 = 60 + (x / 60000000000) % 60;
  let hr: i64 = 24 + (x / 3600000000000) % 24;
  let days: i64 = x / (24*3600000000000);
  let mut yr = 1969; let mut tdays = days; let mut found = 0; 
  while found == 0 {
    let pump = if yr%4==0 { 366 } else { 365 };
    if tdays < -1 * pump {
      found = 1;
    } else {
      tdays = tdays + pump; yr = yr -1;
    }	   
  }
  found = 0; let mut mth = 0;
  while found == 0 {
    let pump = if mth==1 {if yr%4==0 { 29 } else { 28 }} else { ML[mth] };  
    if tdays < -1*pump {
      found = 1;
      tdays = pump + tdays;
    } else {
      tdays = tdays + pump; mth = mth-1;
    }
  }
  return sprintf!("%04ld-%02ld-%02ld %02ld:%02ld:%02ld.%09ld", yr, mth, tdays, hr, min, sec, ns).unwrap();
}

  
impl fmt::Display for Emside {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    match self {
      Emside::Buy => write!(f, "buy"),
      Emside::Sell => write!(f, "sell"),
    }
  }
}
// PSip: PART of SIP state
//   will only be Bid or Ask part of SIP 
// bs: buy or sell
// nv: number of Venues printing SIP
// b2V: Bit Vector representing Venue on state
// q: Quantity at best price level
// p: Price of best price level
// vq: Vector of resting quantities for venues (0 if venue off)
// vp: Vector of resting prices of best bid/ask at venues
// nW: number of venues at Winning price
pub struct PSip {
  pub bs:Emside,
  pub nv: TNV,
  pub b2v: B2Vs,
  pub q: TQ,
  pub p: TP,
  pub vq: Vec<TQ>,
  pub vp: Vec<TP>,
  pub nwv: TNV 
}

impl TSip {
  pub fn verify(&self) -> u32 {
    (*self).b.verify() + (*self).s.verify()
  }
}

// Really only price and quantity come into NBBO
//  There is no price integer vector we will use in real orderbook
impl PSip {
  pub fn new(bs:Emside, n_venues: u8) -> PSip {
    return PSip { bs:bs, nv: n_venues as TNV,
      b2v: B2Vs::new(),
        q: 0 as TQ, p: 0.0 as TP,
        vq: vec![0 as TQ;n_venues as usize], vp: vec![0 as TQ;n_venues as usize], nwv: 0 as TNV 
      };
  }
  pub fn verify(&self) -> u32 {
    let mut n_err:u32 = 0 as u32; let mut sumq:TQ = 0 as TQ;
    let mut tgtp: TP = 0.0 as TP;
    let mut tot_w:TNV = 0.0 as TNV;
    if self.bs == Emside::Buy {
      for ii in 0..self.nv {
        if self.vq[ii as usize] > (0 as TQ) {
          if (tgtp <= 0.0 as TP) || (self.vp[ii as usize] > tgtp) {
            tgtp = self.vp[ii as usize];
          }			   
        }
      }
    } else {
      for ii in 0..self.nv {
        if self.vq[ii as usize] > (0 as TQ) {
          if (tgtp <= 0.0 as TP) || (self.vp[ii as usize] < tgtp) {
            tgtp = self.vp[ii as usize];
          }			   
        }
      }
    }
    if (self.p <= 0.0) && (tgtp <= 0.0) {
    } else if self.p != tgtp {
      println!("PSip verify on bs side {}, self.p={}, maxP={}, nVenues={}.",
        self.bs, self.p, tgtp, self.nv); n_err+=1;		 
    }
    for iv in 0..self.nv {
      if (self.vp[iv as usize] == tgtp) && (tgtp > 0.0) {
        sumq += self.vq[iv as usize];  tot_w+=1;
        if self.b2v.is_on_i(iv as TNV) == 0 {
          println!("PSip verify on {} but on iv={}/{}, its price was {}, best is {} but is_on_i({}) is {}",
            if self.bs == Emside::Buy { "BUY" } else { "SELL" }, 
            iv, self.nv, self.vp[iv as usize], tgtp, iv, self.b2v.is_on_i(iv)); 
          n_err+=1;			 
        }
      } else {
         if self.b2v.is_on_i(iv as TNV) > 0 {
           println!("PSip verfiy on {} on iv={}/{}, its price was {}, best is {} but is_on_i({}) is {}",
            if self.bs == Emside::Buy { "BUY" } else { "SELL" },  iv, self.nv, self.vp[iv as usize], tgtp, iv, self.b2v.is_on_i(iv as TNV));
           n_err+=1;
         }
      }
    }
    if sumq != self.q {
      println!("PSip -- Verify on {} -- fail to have Sum of winning Quantity was {} but self.q={}",
        if (*self).bs==Emside::Buy { "BUY" } else { "SELL" },
          sumq, self.q); 
        n_err+=1;
    }
    if tot_w != self.nwv {
      println!("PSip -- Verify on {} -- fail to have tot Winning, totW = {} but self.nvw={}",
        if self.bs == Emside::Buy { "BUY" } else { "SELL" },
        tot_w, self.nwv); 
      n_err+=1;
    }
    if (n_err as u32) > (0 as u32) {
      let out_st = if self.bs == Emside::Buy { "maxP" } else {"minP" }; 
      println!("PSip verify, there were n_err={} failures on side={} with {}={} ",
        n_err, self.bs, out_st, tgtp);
    }
    return n_err;
  }
  pub fn update(&mut self, newv: TNV, newp: TP, newq: TQ) -> u8 {
     let newvi: usize = newv as usize;
     let newp = if newq > 0 as TQ { newp } else { 0.0 };
     if self.p == newp && self.p <= 0.0 {
     } else if (self.p == 0.0) && (newp > 0.0) {
       self.b2v.zero_all(); self.nwv = 1;
       self.b2v.one_i(newv); self.p = newp; self.q = newq;
     } else if (newp > 0.0) &&
               ((self.bs == Emside::Buy && self.p < newp) ||
                (self.bs == Emside::Sell && self.p > newp)) { 
       self.b2v.zero_all(); self.nwv = 1;
       self.b2v.one_i(newv); self.p = newp; self.q = newq;
     } else if self.p == newp {
        if self.b2v.is_on_i(newv) > 0 {
          self.q += newq - self.vq[newvi];
          self.vq[newvi] = newq; 
        } else {
          self.nwv += 1; 
          self.q += newq;  self.b2v.one_i(newv);
        }
     } else if (((*self).bs == Emside::Buy) && ((*self).p > newp)) ||
               (((*self).bs == Emside::Sell) && ((newp <= 0.0) || ((*self).p < newp))) {
        if (*self).b2v.is_on_i(newv) > 0 {
          if (*self).nwv > 1 {
            (*self).q -= self.vq[newvi]; self.nwv -= 1;
            (*self).b2v.zero_i(newv);
          } else {
            (*self).b2v.zero_all(); (*self).nwv = 0;
            (*self).vp[newvi] = newp;  (*self).vq[newvi] = newq;
            (*self).q = 0 as TQ;
            let mut tgtp: TP = 0.0 as TP;
            for ii in 0..self.nv {
              if (tgtp <= 0.0) && ((*self).vp[ii as usize] > 0.0) {
                tgtp = (*self).vp[ii as usize];
              } else if (*self).bs == Emside::Buy && (*self).vp[ii as usize] > tgtp {
                tgtp = (*self).vp[ii as usize];
              } else if (*self).bs == Emside::Sell && ((*self).vp[ii as usize] > 0.0) && (*self).vp[ii as usize] < tgtp {
                tgtp = (*self).vp[ii as usize];
              }
            }
            if tgtp > 0.0 {
            for ii in 0..(*self).nv {
              if (*self).vp[ii as usize] == tgtp {
                (*self).b2v.one_i(ii); (*self).q += (*self).vq[ii as usize];
                (*self).nwv += 1;
              }
            }
            (*self).p = tgtp;
            } else {
              (*self).p = 0.0; (*self).q = 0 as TQ; (*self).nwv = 0;
            }
          }
        }
     } else {
          // Pass, nothing to do?
          println!("Why are we here in edge case?: newv = {}, newp={}, self.p={}", newv, newp, self.p);
     }
     (*self).vp[newvi] = newp; (*self).vq[newvi] = newq;
     return newv;
  }
}

// TSip, total SIP, combines, 
//   b and s sides to have NBB and NBO calculator
//
//   i_r which record in input data we are on
//   tR which Time in input data we are on
pub struct TSip {
  pub b: PSip,
  pub s: PSip,
  pub i_r: usize,
  pub tt: TT
}
impl TSip {
  pub fn new(n_venues: u8) -> TSip {
    TSip{b:PSip::new(Emside::Buy,n_venues as TNV),
	    s:PSip::new(Emside::Sell,n_venues as TNV),
         i_r: 0 as usize,tt: 0 as TT
    }
  }	  	
}

/////////////////////////////////////////////
// RecSip  -- Column based recording.
//
// RecSip structure records the changes of state
//   NBB is National Best Bid (highest bid price
//   NBO is National Best Offer (lowest sell price)
//
//   nNr: number of total records struct can hold
//   ir: Recording index of current position to write
//   nv: number of venues in sip
//   vT: time vector
pub struct RecSip{ 
   pub n_r: u32,
   pub i_r: usize,
   pub nv: TNV,
   pub vt:Vec<TT>,
   pub vb_nwv: Vec<TNV>, pub vs_nwv:Vec<TNV>,
   pub vb_b2v:Vec<B2Vs>,pub vs_b2v:Vec<B2Vs>,
   pub vb_u: Vec<u64>, pub vs_u:Vec<u64>,
   pub vb_q: Vec<TQ>,pub vb_p: Vec<TP>,
   pub vs_q: Vec<TQ>,pub vs_p: Vec<TP>,
   pub vi_line:Vec<usize>
}

// PublishSIP -- row based recording.
// b_... {buy/Bids/NBB}, s_... {sell/Asks/NBO}
// _q ... Total Quantity at NBB/NBO
// _p ... Price of NBB/NBO
// _nwv ... Number of venues with prices at NBB/NBO
// _b2v ... Bit vector identifying which venues at NBB/NBO
pub struct PublishSip{
  pub tt:TT, 
  pub b_p:TP, pub s_p:TP,
  pub b_q:TQ, pub s_q:TQ,
  pub b_vu:u64, pub s_vu:u64, // We will cast as u64 for Numpy publish
  pub i_line:ILine,
  pub b_nwv:TNV, pub s_nwv: TNV
  //pub b_b2v:B2Vs, pub s_b2v:B2Vs,
}

pub fn u64_string(x:u64, n_rv:u16) -> String  {
    let mut ss = String::new();
    ss.reserve(n_rv.into());
    for ii in (0..n_rv).rev() {
      if (x & ( ( 1 as u64) << ii)) > 0 {
        ss.push('1');
      } else {
        ss.push('0');
      }
    }
    return ss
}
  
impl PublishSip {
  pub fn new() -> Self {
    PublishSip{tt:0 as TT,
      b_p:0 as TP, s_p:0 as TP,
      b_q:0 as TQ, s_q: 0 as TQ, 
	 b_vu:0 as u64, s_vu:0 as u64,
      i_line:0 as ILine,
      b_nwv:0 as TNV,s_nwv: 0 as TNV} 
      //b_b2v: B2Vs::new(), s_b2v: B2Vs::new(),
  }
}
pub struct VPublishSip {
  pub v: Vec<PublishSip>,
  pub on_r:usize
}
impl VPublishSip {
  pub fn new(_n:u32) -> VPublishSip {
    // VPublishSip{v: vec![PublishSip::new();n as usize],on_r:0 as usize}
    //let mut av:VPublishSip = VPublishSip{v: Vec::new(),on_r:0 as usize};
    //av.v.reserve(n.try_into().unwrap());
    VPublishSip{v: Vec::new(),on_r:0 as usize}
  }
}
impl VPublishSip {
  // Note: typically printLine < iLine
  pub fn publish(&mut self, ttime: TT, i_line: ILine, tsip: &TSip) -> () {
    //(*self).v[(*self).on_r].publish(ttime, i_line, tsip); (*self).on_r+=1;
      (*self).v.push(PublishSip{tt:ttime, 
            b_p:(*tsip).b.p, s_p:(*tsip).s.p,
            b_q:(*tsip).b.q, s_q:(*tsip).s.q,
            b_vu:(*tsip).b.b2v.out64(), s_vu:(*tsip).s.b2v.out64(),
            i_line:i_line,
            b_nwv:(*tsip).b.nwv, s_nwv: (*tsip).s.nwv});
      (*self).on_r += 1;
  }
}
/***
// No way to implement "from" idiom in Rust because existing types can't
//  inherit traits? 
impl TryFrom<Vec<usize>> for Vec<u64> {
  type Error = ();
  fn try_from(v_usize:Vec<usize>) -> Result<Self, Self::Error> {
    let mut v_64 = vec![0 u64; v_usize.len()];
    for ii in 0..v_usize.len() {
      v_64[ii] = v_usize[ii] as u64;
    }
    Ok(v_64)
  }
}
***/
// Advantage of this: No need for Result Error handling
pub fn to_v64(v_usize:Vec<usize>) -> Vec<u64> {
  let mut v_64 = vec![0 as u64; v_usize.len()];
  for ii in 0..v_usize.len() {
    v_64[ii] = v_usize[ii] as u64;
  }
  return v_64; 
}

impl Clone for PublishSip {
  fn clone(&self) -> Self {
    PublishSip{tt:self.tt, b_p:self.b_p, s_p:self.s_p,
               b_q:self.b_q, s_q:self.s_q,
			b_vu: self.b_vu, s_vu: self.s_vu,
               i_line:self.i_line,
               b_nwv:self.b_nwv, s_nwv:self.s_nwv,
               //b_b2v:self.b_b2v.clone(), s_b2v:self.s_b2v.clone(),
    }
  }
}
impl PublishSip {
  pub fn publish(&mut self, ttime: TT, i_line: ILine, t_sip: &TSip) -> () {
    (*self).i_line = i_line; (*self).tt = ttime;
    (*self).b_p = (*t_sip).b.p; (*self).s_p = (*t_sip).s.p;
    (*self).b_q = (*t_sip).b.q; (*self).s_q = (*t_sip).s.q;
    (*self).b_nwv = (*t_sip).b.nwv; (*self).s_nwv = (*t_sip).s.nwv;
    //(*self).b_b2v = (*t_sip).b.b2v.clone(); (*self).s_b2v = (*t_sip).s.b2v.clone();
    (*self).b_vu = (*t_sip).b.b2v.out64();  (*self).s_vu = (*t_sip).s.b2v.out64();
  }
}

/****
 Old Numpy based export.  Does not appear to still run in updated Arrow code
// For SIP, so long as we can figure out way to publish activity
unsafe impl Element for PublishSip {
  const IS_COPY: bool = false;  
  fn get_dtype_bound(py: Python<'_>) -> Bound<'_, PyArrayDescr> {
    let locals = [("np", get_array_module(py).unwrap())].into_py_dict_bound(py);
      //let i_str = "import numpy as np; np.dtype([('time','datetime64[ns]'),('price',np.float64),('side','|S5')])";
      let result = py
        .eval_bound("np.dtype([('time','datetime64[ns]'),('nbb','float64'),('nbo','float64'), \
           ('nbbq','float64'),('nboq','float64'), \
           ('II_b','uint64'),('II_s','uint64'),('iline','uint64'),
           ('n_bid','uint8'),('n_ask','uint8')])", Some(&locals), None);
      return result.unwrap().downcast_into::<PyArrayDescr>().unwrap();
    }
    // Weird: because _py lifetime is unimportant, it asks that we call "_py" by underscore
    fn clone_ref(&self, _py: Python<'_>) -> Self {
      return PublishSip{tt:self.tt, 
               b_p:self.b_p, s_p:self.s_p,
               b_q:self.b_q, s_q:self.s_q,
               b_vu: self.b_vu, s_vu: self.s_vu,
               i_line:self.i_line,
               b_nwv:self.b_nwv as TNV, s_nwv:self.s_nwv as TNV,
               //b_b2v:self.b_b2v.clone(), s_b2v:self.s_b2v.clone(), 
               };
    }
    fn vec_from_slice(_py: Python<'_>, slc: &[Self]) -> Vec<Self> { 
      return slc.to_vec();
    }
}  
****/
impl RecSip {
  // nv is number of venues in this case
  pub fn new(n_r:u32, nv: TNV) -> RecSip {
    let n_ur: usize = n_r.try_into().unwrap();
    return RecSip{n_r:n_r, i_r: 0,
       nv:nv,
	  vb_nwv:vec![0 as TNV;n_ur],   vs_nwv:vec![0 as TNV;n_ur],
	  vb_b2v: vec![B2Vs::new();n_ur], vs_b2v: vec![B2Vs::new();n_ur],
          vb_u:vec![0 as u64; n_ur], vs_u: vec![0 as u64;n_ur],
	  vb_q: vec![0 as TQ;n_ur], vb_p: vec![0 as TP;n_ur],
	  vs_q: vec![0 as TQ;n_ur], vs_p: vec![0 as TP;n_ur],
       vt: vec![0 as TT;n_ur], vi_line: vec![0;n_ur]
    }
  }
  // Little helper to determine if no change.
  pub fn state_match(&self, ptsip:&TSip) -> u8 {
    if (*self).nv < ((*ptsip).b.nwv) {
      println!("Error statematch, self nv={}, but ptsip b nwv is {}.",
         (*self).nv, (*ptsip).b.nwv);  return 7;
    }
    if (*self).nv < ((*ptsip).s.nwv) {
      println!("Error statematch, self nv={}, but ptsip s nwv is {}.",
         (*self).nv, (*ptsip).s.nwv);  return 7;
    }
    if ((*ptsip).b.p == self.vb_p[self.i_r]) &&
       ((*ptsip).b.q == (*self).vb_q[(*self).i_r]) &&
       ((*ptsip).b.nwv == (*self).vb_nwv[(*self).i_r]) &&
       ((*ptsip).b.b2v == (*self).vb_b2v[(*self).i_r]) &&
       ((*ptsip).s.p == (*self).vs_p[(*self).i_r]) &&
       ((*ptsip).s.q == (*self).vs_q[(*self).i_r]) &&
       ((*ptsip).s.nwv == (*self).vs_nwv[(*self).i_r]) &&
       ((*ptsip).s.b2v == (*self).vs_b2v[(*self).i_r]) {
      return 1;
    } else {
      return 0;
    }
    //return 128; // UNREACHABLE?
  }
  pub fn record_state(& mut self, ptsip: &TSip) -> u32 {
    if self.state_match(ptsip) > 0 {
      return 0;
    }
    if (*self).i_r >= (*self).n_r.try_into().unwrap() {
      println!("Error record_state:  i_r={} but n_r={}. ",
        (*self).i_r, (*self).n_r);
      return 0;
    }
    (*self).vb_p[(*self).i_r] = (*ptsip).b.p;(*self).vb_q[(*self).i_r] = (*ptsip).b.q;
    (*self).vb_b2v[(*self).i_r] = (*ptsip).b.b2v.clone();(*self).vb_nwv[(*self).i_r] = (*ptsip).b.nwv;
    (*self).vs_p[(*self).i_r] = (*ptsip).s.p;(*self).vs_q[(*self).i_r] = (*ptsip).s.q;
    (*self).vs_b2v[(*self).i_r] = (*ptsip).s.b2v.clone();(*self).vs_nwv[(*self).i_r] = (*ptsip).s.nwv;
    (*self).vb_u[(*self).i_r]= (*ptsip).b.b2v.out64(); (*self).vs_u[(*self).i_r] = (*ptsip).s.b2v.out64();
    (*self).vi_line[(*self).i_r] = (*ptsip).i_r;  // note tsip.i_r >= self.i_r
    (*self).vt[(*self).i_r] = (*ptsip).tt;
    (*self).i_r += 1;
    return (*self).i_r as u32;
  }
}

