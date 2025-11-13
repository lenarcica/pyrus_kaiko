/////////////////////////////////////////////////////////////////////
// b2v_struct.rs
//
//  A Module defining and exploring a Bit mask to represent venue participation
//
//
//use macro_const::{macro_const};
pub const NMAXVENS:u8 = 64;

//macro_const! {
pub const NMAXBITVENS:u8=  (NMAXVENS/8) + (if NMAXVENS % 8 > 0 { 1  } else { 0 }) as u8;
//}
pub type TUV = u8;
pub type B2V = [u8;NMAXBITVENS as usize];

pub struct B2Vs{pub s:B2V}
impl Clone for B2Vs {
  fn clone(&self) -> Self {
    B2Vs{s:(*self).s}
  }
}
impl B2Vs {
  pub fn new() -> Self {
    B2Vs{s:[0;NMAXBITVENS as usize]}
  }
}

// Kernighan's algoritm should determine 1s in bite quickly
fn kern_count(bvp:u8) -> u8 {
  let mut ntt:u8 = 0;
  let mut bvpc = bvp;
  while bvpc > 0 {
    ntt = ntt + 1;
    bvpc = bvpc & (bvpc-1);
  }	  
  return ntt;
}
impl B2Vs {
  pub fn out64(&self) -> u64 {
    let mut o1:u64 = 0; let mut m8:u64=1;
    for ii in 0..NMAXBITVENS {
      o1 += m8 * ( (*self).s[ii as usize] as u64);
      if ii < NMAXBITVENS-1 { m8 = m8 * 256; }
    }
    o1
  }
}

impl B2Vs {
  pub fn zero_all(&mut self) -> () {
    for ii in 0..NMAXBITVENS {
  	 (*self).s[ii as usize] = 0;	
    }
  }
}
impl B2Vs {
  pub fn is_on_i(&self, i_c: TUV) -> u8 {
    if i_c > NMAXVENS.into() {
 	 return 0 as TUV;  
    }
    let div_l:u8 = (i_c / 8).try_into().unwrap(); 
    let rem_l:u8 = (i_c % 8).try_into().unwrap();
    if ( ((*self).s[div_l as usize]) & ( ( 1 as u8) << rem_l)) > 0 {
 	 return 1 as u8;  
    }
    return 0 as u8;
  }
}

impl B2Vs {
  pub fn zero_i(&mut self, i_c: TUV) -> ()  {
    let div_l: u8 = (i_c / 8).try_into().unwrap();
    let rem_l: u8 = (i_c % 8).try_into().unwrap();
    if rem_l == 0 {
      (*self).s[div_l as usize] = (*self).s[div_l as usize]- (1 as u8);
    } else {
      (*self).s[div_l as usize] = (*self).s[div_l as usize] & !(1 << ((rem_l)));
    }	
  }
}

impl B2Vs{
  pub fn one_i(&mut self, i_c: TUV) -> () {
    let div_l: u8 = (i_c / 8).try_into().unwrap();
    let rem_l: u8 = (i_c % 8).try_into().unwrap();
    (*self).s[div_l as usize] = (*self).s[div_l as usize] | ( (1 as u8) << rem_l);  
  }
}
impl PartialEq for B2Vs {
    // Required method
    fn eq(&self, other: &Self) -> bool {
      for ii in 0..NMAXBITVENS {
        if (*self).s[ii as usize] != (*other).s[ii as usize] {
          return false;
        }
      }
      return true;
    }

}

impl B2Vs{
  // Calculate total number of venues at max
  pub fn n_on_total(&self, n_rv: u16) -> u32 {
    let mut cnt:u32 = 0;
    let n_wrv = n_rv / 8 + (if (n_rv % 8) > 0 { 1 } else { 0 } );
    for ii in 0..n_wrv {
	 cnt += kern_count((*self).s[ii as usize]) as u32;
    }
    cnt
  }

  pub fn p_string(&self, n_rv:u16) -> String  {
    let mut ss = String::new();
    ss.reserve(n_rv.into());
    for ii in (0..n_rv).rev() {
      let val_e = (*self).is_on_i(ii.try_into().unwrap());
 	 ss.push(if val_e > 0 { '1' } else { '0'});
    }
    //unsafe{ *(ss + nRV) = '\0' };
    return ss
  }

  pub fn bit_tester( &mut self, n_rv: u16) -> ()  {
  print!("BitTester ------  EXECUTE -----------------------------------------------------------------------\n");	
  print!("---  Intitiate  nRV={}. NMAXVENS={}.\n", n_rv, NMAXVENS);
  print!("--- Number Of Bits On are {}. \n", (*self).n_on_total(n_rv));
  print!("--- We Zero ourselves out. \n");
  self.zero_all();
  print!("--- Now after Zero, number of bits are {}. \n", (*self).n_on_total(n_rv));
  print!("Now we will turn on bits 2 and 3. \n");
  (*self).one_i(2); (*self).one_i(3);
  print!(" --- After that number of bits are {}, which is string {}. \n", (*self).n_on_total(n_rv), (*self).p_string(n_rv));
  print!(" --- Turn on 9, turn on 10, turn of 2 \n");
  (*self).one_i(9); (*self).one_i(10); (*self).zero_i(2);
  print!(" --- After that number of bits are {}, which is string {}. \n", (*self).n_on_total(n_rv), (*self).p_string(n_rv));
  print!(" --- turn 2 on, turn 11 on, turn 10 off. \n");  
  (*self).one_i(11); (*self).zero_i(10); (*self).one_i(2);
  print!(" --- After that number of bits are {}, which is string {}. \n", (*self).n_on_total(n_rv), (*self).p_string(n_rv));
  }
}
