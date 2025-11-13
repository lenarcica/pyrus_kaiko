Default_Fake_OB = """
nMidPoints=5000;nGen=4000;nVen=20;verbose = 1;
stTime = '2024-10-03 10:30:32'; deltaSec= 1.53; startMid = 100;
dTas = 4; sd_Midpoint=4; spread_Orders=50; AddMidPrice=True; v_bhit = False; fast_check="FASTCHECK";
from pyrus_marketbook import pyrus_marketbook; 
"""
def Fake_OB(nGen:int, nMidPoints:int=5000,nVen:int=10,verbose:VerbosityLevel = 0, 
  stTime:str='2024-10-03 10:30:32', deltaSec:np.float64=1.53, 
  startMid:np.float64=100, sd_Midpoint:np.float64=4, spread_Orders:np.float64=50, v_bhit:bool = False, 
  AddMidPrice:bool = False, fast_check:str="FASTCHECK", fixed_safe_n:bool=False) :
  """
  Fake_OB
  
    A basic way to simulate market participant quotes that reflect a moving market.

    This is an artificial model, designed to be a simple way to create a diverse set of long-lasting
     quotes that conform to a volatile market.  Instead of generating customer quotes from stochastic agent
     process, we simulate first the midpoint price from a stationary model, and then simulate quotes at random
    times and distances away from the midpoint price.

   An issue though, is that the quotes may eventually be hit by midpoint motion.  Here we can use the verify_price
     algorithm to identify the first midpoint time where the price is crossed.

  ##############################################################################
  ## Simulating event times for customers giving multiple sides of trades.
  ##   
  ## Start times randomly uniform in time window
  ## Then simulate side for the limit order, a random quantity
  ## Then place the price of the limit order, starting somewhere are the right side of the midpoint.
  ## Then come up with a criterion to cancel the limit order
  ## 1. NBBO comes up to within v_b of the order
  ##     If "v_b=0" this is like a cancel just when order should have traded, others assume order cancels somewhat before that.
  ## 2. NBBO falls away a distance of v_c for the order.
  ##     This is a "give up hope" cancel.  One might reasonably believe a real actor would adjust their orders price to increase chance of fill.
  ## 
  ##  We add v_b, v_c, v_p to a measured criterion.  We will use a Verify Price algorithm to identify when first hit of v_b, v_c occurs
  v_t0 = np.sort(np.random.uniform(np.min(vTime).astype(np.int64), np.max(vTime).astype(np.int64), nMake)).astype('datetime64[ns]');
  """
  import pandas as pd; import numpy as np; import copy;
  from pyrus_marketbook.ob import BasicSimulateMidPrice, TimeGeneration;
  from pyrus_marketbook.ob import asjoin_py, unsorted_asjoin_py;
  #from pyrus_marketbook.ob import * asjoin_py
  StT = "Fake_OB(Vb=" + str(verbose) + ",nG=" + str(nGen) + "): ";
  if verbose >= 1 :
    print(StT + " -- Begin");
  # If we wanted fully random simulated data we might need to simulate more
  #  This is because we want to simulate a healthy amount of orders that started before v_t0 so that on window open
  #  many orders are already on the book.
  # However, some of those orders will be closed before the window starts too.  And it will be challenging to simulate
  #  right number of orders that survive window time into market open.  Hence, we simulate about 50% more orders
  nMake = int(np.floor(nGen * 1.5).astype(np.int64)); # simulate some orders starting before and after intended start time.
  vTmin, vTmax, vTime = TimeGeneration(stTime=stTime, nMidPoints=nMidPoints, deltaSec=deltaSec);
  MidPrice = BasicSimulateMidPrice(nMidPoints, sd_Midpoint, startMid);
  MidTable = pd.DataFrame({'time':vTime, 'nbm':MidPrice});
  ## Simulate Orders: #1 Simulate initial time, side, quantity
  v_t0 = np.sort(np.random.uniform(np.min(vTime).astype(np.int64), np.max(vTime).astype(np.int64), nMake)).astype('datetime64[ns]');
  v_s = np.random.choice([-1,1], size=nMake).astype(np.int8); ## v_s, vector of sides
  v_q = np.round(np.random.uniform(0,1000,size=nMake));       ## v_s vector of quantity
  # venue/participant indicator
  v_ven = np.random.choice(np.arange(nVen), size=nMake).astype(np.int64);
  # Midpoint price as of order simulated arrival time: (gained by as_of_join)
  n_mid_t =  pyrus_marketbook.asjoin_py(v_t0.astype(np.int64),vTime.astype(np.int64))
  mid_of_v = MidPrice[n_mid_t];                               ## v_s vector of our data
  if verbose >= 1 :
    print(StT + " -- Fake_OB Generated mid_of_v, now makign add points"); 
  # Simulate a deviation from midpoint "add_p", which can be zero/0.01
  add_p = (1 + np.round(np.abs(np.random.normal(0,spread_Orders,nMake))));
  add_p[add_p < 0] = 0.0;
  v_p = add_p*(2*(v_s < 0) - 1.0) + mid_of_v;  # simulate order price
  # Forced Cancellation:
  #  Trivially, v_b sets a level for which, if midpoint comes within v_b of the order, cancels the order.  Typically this happens when Midpoint=v_p
  #             v_c is more of a "I give up" situation, so if midpoint reduces to be more than v_c away from rest price, cancel the order.
  if v_bhit == True :
    v_b = 0.0 * add_p; # in a "v_bhit" simulation we trivially set cancel times to be midpoint reacging between 0, or departing by more than 5 add_p
    v_c = 5 * add_p;
  else :
    v_b = np.round(np.random.uniform(0,add_p,len(add_p)));
    v_c = add_p*2 + np.round(np.abs(np.random.normal(0,spread_Orders,nMake)));
    v_c[v_c < add_p] = add_p[v_c < add_p]  # if v_c=add_p we would be cancelling order at moment it is listed (which happens in real life but is not useful for simulation)
  # simulate initial survival times uniform as length of the window
  survival_time = np.random.uniform(np.min(vTime).astype(np.int64), np.max(vTime).astype(np.int64), size = nMake).astype('datetime64[ns]') - np.min(vTime);
  v_t1 = v_t0 + survival_time.astype('timedelta64[ns]'); v_t1[v_t1 > vTmax] = vTmax; ## Everything out of window just ends at widnow end
  v_t1c = copy.copy(v_t1);  # keep copy of original simulation time
  market_price_v_t1c = MidPrice[unsorted_asjoin_py(v_t1c.astype(np.int64), vTime.astype(np.int64))]
  ## Forcing numpy to type a variable at i8 can sometimes be a challenge
  ilhverbose = np.array([verbose - 3]).astype(np.int8)[0]
  if verbose >= 1 :
    print(StT + ": testing v_t1lc length " + str(len(v_t1c)) + " for data against midpoints at " + str(len(MidPrice)) + " points.");
  ## "pyrus_verify_price" algorithm used to determine whether a v_b/v_c level event happens before the simulated close.
  FASTCHECK_ANSWERS = ["FASTCHECK", "FAST_CHECK", "fCheck", "FCheck", "FCHECK", "fast_check", "fastcheck", "AUTOCHECK", "autocheck", "Autocheck"];
  PFAST_CHECK_ANSWERS =  ["PARALLELFASTCHECK","PARRALEL","PARALLEL","PARALLELFASTCHECK", "PARRALLEL", "PARRALEL"];
  SLOW_CHECKANSWERS = ["SLOWCHECK", "slow_price_verify", "SLOW_PRICE_VERIFY"];
  SLOWEST_LIMIT_HIT_ANSWERS = ["LIMIT_HIT", "limithit", "limit_hit", "lhit", "LHIT", "limhit", "Slowest", "SLOWEST"];
  if fast_check in FASTCHECK_ANSWERS+PFAST_CHECK_ANSWERS+SLOW_CHECKANSWERS : 
  if fast_check in FASTCHECK_ANSWERS+PFAST_CHECK_ANSWERS+SLOW_CHECKANSWERS : 
    import pyarrow as pa; pdiff = 0;
    tabOrd = pd.DataFrame({'bs':np.asarray(['b' if x == 1 else 's' for x in v_s], dtype="S1"), "price":v_p, "qty": np.zeros(len(v_p)) + 2,
                           'open': v_t0.astype('datetime64[ns]'), 'close': v_t1.astype('datetime64[ns]'), 'origidx':np.arange(len(v_s))});
    tabOrd['price'] = v_p + v_s * v_b;
    nbboBig = MidTable.copy(); nbboBig.columns = ['time','nbb'];  nbboBig['nbo'] = nbboBig['nbb']; 
    DoFast = 1 if fast_check in FASTCHECK_ANSWERS else (2 if fast_check in PFAST_CHECK_ANSWERS else 0);
    if verbose >= 1 :
      print(StT + " -- Trying to create minimum usage times using pyrus_verify_prices algorithms.");
    ## Verfiy a hit of v_p + v_s*v_b
    v_t1_b = pyrus_marketbook.pyrus_verify_prices(pdiff, pa.Table.from_pandas(tabOrd), pa.Table.from_pandas(nbboBig),
                    False, False, ilhverbose,
                    DoFast, False).to_numpy();
    tabOrd['price'] = v_p + v_s * v_c;  ## Note v_c is a point far away (on opposite side of NBBO) from v_p where we give up order
    tabOrd['bs'] = np.asarray(['s' if x == 1 else 'b' for x in v_s], dtype='S1');
    ## Verify a walk away up to v_p + v_s * v_c;
    v_t1_c = pyrus_marketbook.pyrus_verify_prices(pdiff, pa.Table.from_pandas(tabOrd), pa.Table.from_pandas(nbboBig),
                    False, False, 0,
                    DoFast, False).to_numpy();
    ## Admittedly we are still in trouble since it is hard to merge data when null matches exist.
    v_t1_both = copy.copy(v_t1_b);  v_t1_both[v_t1_b > v_t1_c] = v_t1_c[v_t1_b > v_t1_c]
    v_t1o_Both = pd.Series(np.repeat(-1, len(v_t1_both)).astype('datetime64[ns]'))
    v_t1o_Both[v_t1_both < len(vTime)] = vTime[v_t1_both[v_t1_both < len(vTime)]]; # v_t1o_new = copy.copy(v_t1o)
    v_t1o_Both[v_t1_both >= len(vTime)] = pd.NaT;
    v_t1o = v_t1o_Both
  elif fast_check in SLOW_CHECKANSWERS :
    v_t1o_Old = np.asarray(pyrus_marketbook.limit_hit(v_t0.astype(np.int64), v_t1c.astype(np.int64), 
      v_s.astype(np.int8), v_p.astype(np.float64), v_b.astype(np.float64), v_c.astype(np.float64), 
      vTime.astype(np.int64), MidPrice.astype(np.float64), ilhverbose)).astype('datetime64[ns]');
    v_t1o = v_t1o_Old
  else :
    print("fast_check = " + fast_check + " is invalid, will trigger an error.");
  if verbose >= 1 :
    print(StT + " we calculated MiPrices, now we will join our v_t1o hits to MidPrice table to identify points.");
  ##
  aHit = np.asarray((v_t1o < v_t1c) & (v_t1o >= v_t0)).astype(np.int8);  # How many orders get hit by NBBO
  if verbose >= 1:
    print(StT + " after that we hit " + str(aHit) + "/" + str(len(v_p)) + " prices, next we perform a shrinking");
  if np.sum(aHit) > 0 :
    ## Update the time to be a uniform somewhere between the initiation time v_t0 and the hit time v_t1o;
    v_t1[aHit == 1] = np.random.uniform(0,1.0,np.sum(aHit)) * (v_t1o[aHit==1]-v_t0[aHit==1]) + v_t0[aHit==1]
  ## New Midpoint price calculagted at unsorted_asjoin to v_t1
  n_mid_t1 =  unsorted_asjoin_py(v_t1.astype(np.int64),vTime.astype(np.int64))
  mid_of_v1 = MidPrice[n_mid_t1];
  dictS = { 1:'b', -1:'s' }
  order_DB = pd.DataFrame({
    "bs": np.asarray([dictS[a] for a in v_s],dtype="|S1"),
    "price": v_p, "mr":v_ven, "qty": v_q, "t0": v_t0, "t1": v_t1, "hit":aHit, "sim_t1": v_t1c, "v_c":v_c, "v_b":v_b, "mid_of_v": mid_of_v, "n_mid_t":n_mid_t,
    "v_s":v_s, "mid_of_v1":mid_of_v1,"mid_of_vt1c": market_price_v_t1c
  });
  ## Eliminate the orders that did not survive to occur in window.
  order_DB = order_DB[(order_DB.t0 <= vTmax) & (order_DB.t1>= vTmin)].reset_index(inplace=False, drop=True);
  if fixed_safe_n == True :
    ## This logic limits that if nGen was target we have up to nGen/2 buys and sells,
    ## Typically this is only a thing to force if we are worried about too many "safe" simulations.
    ## Of course, the "hits" are imaginary (we rederived v_t1o to guarantee pre-hit times are in v_t1, but v_t1c is a close time that "hits"
    ## this helps us test our algo for its ability to find a mix of data that indeed got hit by v_b/v_c criterion
    numSafe = len(order_DB) - np.sum(order_DB['hit']);
    numSafeSells = len(order_DB[(order_DB['v_s'] < 0) & (order_DB['hit']==0)]);
    numSafeBuys = len(order_DB[(order_DB['v_s'] > 0) & (order_DB['hit']==0)]);
    nNeed = round(nGen/2);
    if (numSafeSells > nNeed) and (numSafeBuys > nNeed) :
      nGetS = nNeed if numSafeSells > nNeed else numSafeSells; 
      nGetB = nNeed if numSafeBuys > nNeed else numSafeBuys;
      order_DB = pd.concat([order_DB[(order_DB['hit']==0) & (order_DB['v_s'] < 0)].reset_index(
                              inplace=False, drop=True)[0:nGetS].reset_index(inplace=False, drop=True),
                            order_DB[(order_DB['hit']==0) & (order_DB['v_s'] > 0)].reset_index(
                              inplace=False, drop=True)[0:nGetB].reset_index(inplace=False, drop=True),
                            order_DB[order_DB['hit']==1].reset_index(inplace=False, drop=True)]);
      
  ##  order_DB = order_DB[0:nGen].reset_index(inplace=False, drop=True);
  vt0 = order_DB.t0.to_numpy();
  vt0[vt0 < vTmin] = vTmin-1000;
  vt1 = order_DB.t1.to_numpy();
  vt1[vt1 > vTmax] = vTmax+1000;
  order_DB.t0 = vt0;  order_DB.t1 = vt1;
  if AddMidPrice == True :
    return(order_DB, MidTable);
  return(order_DB);
