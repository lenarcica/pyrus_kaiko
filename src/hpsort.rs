use std::vec::Vec;

macro_rules!triomax01 {
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

pub fn idx_sort_trio(v_s: Vec<char>, v_p: Vec<f64>, v_v: Vec<i64>, v_t:Vec<i64>) -> Vec<usize> {
  println!("idx_sort_trio: we start with v_s.len() = {}, v_t.len() = {}.",
    v_s.len(), v_t.len());
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
  println!("Hey who is max 0 or 3? triomax02(0,3)={}",
    (triomax02!(0,2)));
  let n: usize = v_s.len() as usize;
  println!("idx_sort_trio: n = {}",n);
  let mut i: usize = 0; let mut ir: usize = 0;
  let mut j: usize = 0; let mut l: usize = (n >> 1) + 1;
  println!("idx_sort_trio: start with l = {}, n={}",l,n);
  let mut rra_ix = 0;
  ir = n-1;
  let mut out_ix:Vec<usize> = vec![0 as usize;n];
  for ii in 0..n { out_ix[ii as usize] = ii as usize; } 

  print!(" -- v_S=[");
  if n > 1 {
    for ii in 0..(n-1) { print!("{},", v_s[ii]); }
  }
  print!("{}]", v_s[n-1]);

  if n < 2 { return out_ix; }
  println!(" -- Start Algo, n = {}, l = {}, ir={}", n, l,ir);
  let mut nloop = 0;
  while true {
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
      if (j < ir) && (1==(triomax02!( (out_ix[j]), (out_ix[j+1]) ))) {
        j=j+1;
      }
      if 1 == (triomax02!( (rra_ix), (out_ix[j]))) {
        out_ix[i] = out_ix[j];
        i=j;
        j = ((j+1) << 1)-1;
      } else {
        break;
      }
      println!(" -- inner loop nloop={}, l={},j={}, i={}, ir={}, rra_ix={} ", nloop, l, j, j, ir, rra_ix);
      nloop = nloop + 1;
      if nloop > n*n {
        break;
      }
    }
    out_ix[i] = rra_ix;
    println!("nloop={}, l={}, ir={}, i={}, j={}, rra_ix={}", nloop, l, ir, i, j, rra_ix);
    if j < n {
      println!("out_ix[j={}] is {}, rra_ix < out_ix[j] = {}", 
      j, out_ix[j],
      triomax02!( (rra_ix), (out_ix[j]) ));
    } else {
      println!("out_ix[ir={}] is {},  rra_ix < out_ix[ir] = {}",
        ir, out_ix[ir], triomax02!( (rra_ix), (out_ix[ir])));
    }
    nloop = nloop + 1;
    if nloop > n*n {
      break;
    }
  }
  return out_ix;
}



pub fn idx_sort_sidepricetime(v_s: Vec<char>, v_p: Vec<f64>, v_t:Vec<i64>) -> Vec<usize> {
  println!("idx_sort_trio: we start with v_s.len() = {}, v_t.len() = {}.",
    v_s.len(), v_t.len());
macro_rules!m_sptmax {
  ($i0:expr, $i1:expr) => {

    (      if v_s[ $i0 ] < v_s[ $i1 ] { 1 
    } else if v_s[ $i0 ] > v_s[ $i1 ] { 0
    } else if v_p[ $i0 ] < v_p[ $i1 ] { 1
    } else if v_p[ $i0 ] > v_p[ $i1 ] { 0
    } else if v_p[ $i0 ] < v_p[ $i1 ] { 1
    } else if v_t[ $i0 ] < v_t[ $i1 ] { 1 
    } else {0})
  
  }
}
  println!("Hey who is max 0 or 3? m_sidemax(0,3)={}",
    (m_sptmax!(0,2)));
  let n: usize = v_s.len() as usize;
  println!("idx_sort_sidepricetime: n = {}",n);
  let mut i: usize = 0; let mut ir: usize = 0;
  let mut j: usize = 0; let mut l: usize = (n >> 1) + 1;
  println!("idx_sort_sidepricetime: start with l = {}, n={}",l,n);
  let mut rra_ix = 0;
  ir = n-1;
  let mut out_ix:Vec<usize> = vec![0 as usize;n];
  for ii in 0..n { out_ix[ii as usize] = ii as usize; } 

  print!(" -- v_S=[");
  if n > 1 {
    for ii in 0..(n-1) { print!("{},", v_s[ii]); }
  }
  print!("{}]", v_s[n-1]);

  if n < 2 { return out_ix; }
  println!(" -- Start Algo, n = {}, l = {}, ir={}", n, l,ir);
  let mut nloop = 0;
  while true {
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
      if (j < ir) && (1==(m_sptmax!( (out_ix[j]), (out_ix[j+1]) ))) {
        j=j+1;
      }
      if 1 == (m_sptmax!( (rra_ix), (out_ix[j]))) {
        out_ix[i] = out_ix[j];
        i=j;
        j = ((j+1) << 1)-1;
      } else {
        break;
      }
      println!(" -- inner loop nloop={}, l={},j={}, i={}, ir={}, rra_ix={} ", nloop, l, j, j, ir, rra_ix);
      nloop = nloop + 1;
      if nloop > n*n {
        break;
      }
    }
    out_ix[i] = rra_ix;
    println!("nloop={}, l={}, ir={}, i={}, j={}, rra_ix={}", nloop, l, ir, i, j, rra_ix);
    if j < n {
      println!("out_ix[j={}] is {}, rra_ix < out_ix[j] = {}", 
      j, out_ix[j],
      m_sptmax!( (rra_ix), (out_ix[j]) ));
    } else {
      println!("out_ix[ir={}] is {},  rra_ix < out_ix[ir] = {}",
        ir, out_ix[ir], m_sptmax!( (rra_ix), (out_ix[ir])));
    }
    nloop = nloop + 1;
    if nloop > n*n {
      break;
    }
  }
  return out_ix;
}


