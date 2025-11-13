# kaiko.py
#
# Alan Lenarcic, 2025/08/06
#
#  Algorithms around faking Kaiko Data, as well as opening Kaiko "fob" (full orderbook) files in gz format

import pandas as pd; import numpy as np;
import pyarrow as pa;
from pyrus_kaiko import pyrus_kaiko;
try :
  from pyrus_marketbook import pyrus_marketbook
  from pyrus_marketbook import ob;
except :
  print("pyrus_kaiko->kaiko.py --- Error, pyrus_marketbook did not load in.  Need to also include.");

# VL = Verbosity Level
import typing
from typing import TypeAlias;

VerbosityLevel: TypeAlias = int;


## fake_for_algo(nMidpoints=400,nGen=1000,nVen=20,verbose = 1,
##  stTime = '2024-10-03 10:30:32', deltaSec= 1.53, startMid = 100,
##  dTas = 4)

def Gen_Fake_Mktbook(nMidpoints=4000, nGen=10000, nVen=10, verbose=1,
  stTime = '2025-07-25 10:00', deltaSec=1, startMid = 1000) :
  ## Hard part will be to get certain price updates to occur at same time
  from pyrus_marketbook import pyrus_marketbook;
  from pyrus_marketbook import ob;
  import numpy as np; import pandas as pd;
  import pyarrow as pa;
  kalgoint1 = np.array([1], dtype=np.uint8)[0];
  mktbook, nr, v_d_fake, bprices, sprices, \
    write_sum_q, write_sum_pq, write_avp, write_worstpi, \
    write_nlevels, write_atq, \
    v_w, verbose, time_minmax = \
    ob.fake_for_algo(nMidpoints=4000,nGen=10000,nVen=10,verbose=1,
    stTime = '2025-07-25 10:00:00', deltaSec=1, startMid=100);
  kalgoint1 = np.array([1], dtype=np.uint8)[0];

  mktbook['time'] = np.floor(10000*np.floor(mktbook['time']/10000)).astype(np.int64).astype('datetime64[ns]').astype(np.int64)
  ctime = mktbook['time'].shift(-1,fill_value=pd.NaT).to_numpy();
  ctime[(mktbook['price'] != mktbook['price'].shift(-1,fill_value=-1)) | (mktbook['side'] != mktbook['side'].shift(-1,fill_value=b'n'))] = pd.NaT;
  mktbook['ctime'] = ctime;
  return(mktbook);

def genASnapShot(mktbook, snaptime) :
  """
  SnapShots

   By themselves, Snapshots can be queried from a pandas dataframe relatively quickly using basic Pandas Filters.
   That said, it is likely polars/duckdb filters could be faster or scale better.  It would be however needlessly
   complex to encode snapshot generation by hand.

   We filter for all snapshots where the opentime and close time straddle our snaptime.
   We also only need non-zero quantity in the snapshot (excepting that the zero print occurs exactly at our snaptime)

   The real work, after filtering is to create the somewhat complext "[[100.0,30],...." strings for the single row
   Returned by a Snapshot.  These strings however are relatively easy to join with built in python, and snapshots 
   Typically will include hundreds of price levels and be unwieldingly large strings.
   Converting prow to a String type appears to be a challenge
  """
  ## Note we want to include even "turn offs" at a snapshot
  allPs = mktbook[(mktbook['time'] <= snaptime) & (mktbook['ctime'] >= snaptime) & ((mktbook['time'] == snaptime) | (mktbook['cm_qty'] > 0))].reset_index(inplace=False, drop=True)
  askT = allPs[(allPs['side'] == b's') | (allPs['side'] == 's') ].reset_index(inplace=False, drop=True)
  askQ = "[" + ",".join(["[" + str(np.round(askT['price'][i]/1000,3)) + "," + str(np.round(askT['cm_qty'][i]/10,1)) + "]" for i in range(len(askT))]) + "]"
  bidT = allPs[(allPs['side'] == b'b') | (allPs['side'] == 'b') ].reset_index(inplace=False, drop=True)
  bidQ = "[" + ",".join(["[" + str(np.round(bidT['price'][i]/1000,3)) + "," + str(np.round(bidT['cm_qty'][i]/10,1)) + "]" for i in range(len(bidT))]) + "]"
  prow = pd.DataFrame({
    'timestamp': np.asarray([snaptime]),
    'type': np.asarray([b's']),
    'asks': np.asarray([askQ]),
    'bids': np.asarray([bidQ])})
  prow['type'] = prow['type'].astype(pd.StringDtype());
  return(prow);


def genAupdate(mktbook, uptime) :
  """
   Note genAupdate is rather slow, only generates one line summary of all times associated with uptime (occur at uptime)
   Comparison to snapshot is that this must only give updates localized at the report time
  """
  allPs = mktbook[(mktbook['time'] == uptime) ].reset_index(inplace=False, drop=True)
  askT = allPs[(allPs['side'] == b's') | (allPs['side'] == 's') ].reset_index(inplace=False, drop=True)
  askQ = "[" + ",".join(["[" + str(np.round(askT['price'][i]/1000,3)) + "," + str(np.round(askT['cm_qty'][i]/10,1)) + "]" for i in range(len(askT))]) + "]"
  bidT = allPs[(allPs['side'] == b'b') | (allPs['side'] == 'b') ].reset_index(inplace=False, drop=True)
  bidQ = "[" + ",".join(["[" + str(np.round(bidT['price'][i]/1000,3)) + "," + str(np.round(bidT['cm_qty'][i]/10,1)) + "]" for i in range(len(bidT))]) + "]"
  prow = pd.DataFrame({
    'timestamp': np.asarray([uptime]),
    'type': np.asarray([b'u']),
    'asks': np.asarray([askQ]),
    'bids': np.asarray([bidQ])})
  return(prow)

def Kaiko_FOB_format_mktbook(mktbook, nSnaps=20, verbose:int =1, 
  printEvery:int=1000,  dec_price:int=0, dec_qty:int=0) :
  ## Note we need to show updates that "turn off" quantity as much as we need non turn offs.  Zeros are good

  from pyrus_kaiko import pyrus_kaiko
  if 'time' in mktbook.columns :
    cntime = 'time';
  elif 'timestamp' in mktbook.columns :
    cntime = 'timestamp';
  times = np.unique(mktbook[cntime])

  ## We will take snapshots at some mkt book events but also ignore at others
  snaptimes = np.random.choice(times, replace=False, size=nSnaps)
  nonsnaptimes = np.asarray([x for x in times if x not in snaptimes])

  mkt2 = mktbook[mktbook[cntime].isin(nonsnaptimes)].reset_index(inplace=False, drop=True).sort_values('time').reset_index(inplace=False, drop=True)
  idtime = mkt2[cntime].to_numpy().astype(np.int64);
  bs01 = np.repeat(0, len(mkt2)).astype(np.uint8)
  bs01[np.asarray(mkt2['side'].isin(['s',b's', b'S'])).astype('bool')] = 1;
  idprice = mkt2['price'].to_numpy().astype(np.uint64);
  idqty = mkt2['cm_qty'].to_numpy().astype(np.uint64);

  ## Generate large amount of "update" lines in rust
  UFob = pyrus_kaiko.kaiko_make_u_fob(3,1,verbose,printEvery, 
    idtime, bs01, idprice, idqty) 
  NF = UFob.to_pandas().copy();
  
  for k in range(len(snaptimes)) :
    t1 = snaptimes[k];
    NF = pd.concat([NF, genASnapShot(mktbook,t1)]);
    k = 0;
  NF['timestamp'] = NF['timestamp'].astype(pd.StringDtype());
  return(NF.sort_values('timestamp').reset_index(inplace=False, drop=True));

def Kaiko_FOB_save_fobbook(fobbook, saveLoc = "") :
  if (saveLoc is not None) and (isinstance(saveLoc, str)) and (saveLoc != "") :
    if verbose >= 1:
      print("Kaiko_FOB_save_fobbook -- Saving fobbook length " + str(len(fobbook)) + " to " + saveLoc);
    try :
      fobbook.to_csv(filepath, sep=";", header=True, index=False, compression='gzip')
    except Exception as e :
      print("KaikoFOBformatOutT -- Error occured on to_csv");
      print(e);
      print(" ---- Recover error for bad format?");
  return(1);

def KaikoReadFOB(gzLoc: str="", verbose:int = 0) :
  """
  KaikoReadFOB(gzLoc, verbose)  -- 

    Uses Kaiko Rust algorithm to read in a gzLoc algorithm which generates a PyArrow Table

    We then use Polars to convert the Arrow Table into a usable table.
    This will include char-string formatting that is difficult to transfer from Rust to Python

    Unfortunately, we see how difficult here it is to convert Polars data (especially Decimal format qty/price)
    into something that can be studied for Diff Time/Diff Quantity that we need to do time series algorithms

    It is possible DuckDB may be better at formatting this data, but it might lose more type characteristics.
  """
  from pyrus_kaiko import pyrus_kaiko; import polars as pl; import numpy as np;
  import pandas as pd; import pyarrow as pa;
  mypat = pyrus_kaiko.kaiko_fob(gzLoc, 0);
  mypl = pl.from_arrow(mypat);
  mypl = mypl.sort(['bs01','price','time'], descending=False)
  mclose = pd.Series(mypl['time'].shift(-1).fill_null(np.int64(-1).astype('datetime64[ns]')).to_numpy().astype('datetime64[ns]'));
  mclose[mclose == np.int64(-1).astype('datetime64[ns]')] = pd.NaT;
  mclose[((mypl['bs01'].to_numpy() != mypl['bs01'].shift(-1).fill_null(2).to_numpy()) | 
          (mypl['price'].to_numpy() != mypl['price'].shift(-1).fill_null(-1).to_numpy())).astype('bool')] = pd.NaT;
  mside = np.repeat('b', len(mypl));
  mside[mypl['bs01'] == 1] = 's';  mside = pd.Series(mside).astype(pd.StringDtype());
  mtype = np.repeat('u', len(mypl));
  mtype[(mypl['us01'] == 1)] = 's';
  mynq = mypl['qty'].shift(1).fill_null(0.0).to_numpy().astype(np.float64);
  mynq[((mypl['price'].to_numpy() != mypl['price'].shift(1).fill_null(-1.0).to_numpy()) |
        (mypl['bs01'].to_numpy() != mypl['bs01'].shift(1).fill_null(2).to_numpy())).astype('bool')] = 0.0;
  my_onq = mypl['qty'].cast(pl.Float64).to_numpy() - mynq;
  mypl = mypl.with_columns( [pl.Series(mclose).alias('close'),
                             pl.Series(mypl['time'].to_pandas().astype('datetime64[ns]')).alias('open'),
                             pl.Series(mside).cast(pl.String).alias('side'),
                             pl.Series(mtype).cast(pl.String).alias('type'),
                             pl.Series(my_onq).alias('onqty')]);
  mypl = mypl.filter(pl.col('onqty').ne(pl.lit(0.0)))
  mynq = mypl['qty'].shift(1).fill_null(0.0).to_numpy().astype(np.float64);
  mynq[((mypl['price'].to_numpy() != mypl['price'].shift(1).fill_null(-1.0).to_numpy()) |
        (mypl['bs01'].to_numpy() != mypl['bs01'].shift(1).fill_null(2).to_numpy())).astype('bool')] = 0.0;
  my_onq = mypl['qty'].cast(pl.Float64).to_numpy() - mynq;
  mypl = mypl.with_columns( [pl.Series(my_onq).alias('onqty')])
  
  mypl = mypl.select(["type","side","price","open","close","qty","onqty"]);
  return(mypl);
                             

def genASnapShot(mktbook, snaptime) :
  ## Note we want to include even "turn offs" at a snapshot
  cntime = 'time';
  cltime = 'ctime';
  if ('timestamp' in mktbook.columns) and ('time' not in mktbook.columns) :
    cntime = 'timestamp';
  elif ('open' in mktbook.columns) and ('time' not in mktbook.columns) :
    cntime = 'open';
  if ('ctimestamp' in mktbook.columns) and ('ctime' not in mktbook.columns) :
    cltime = 'ctimestamp';
  elif ('close' in mktbook.columns) and ('ctime' not in mktbook.columns) :
    cltime = 'close';
  cqty = 'cm_qty';
  if ('mqty' in mktbook.columns) and ('cm_qty' not in mktbook.columns) :
    cqty = 'mqty';
  elif ('qty' in mktbook.columns) and ('cm_qty' not in mktbook.columns) :
    cqty = 'qty';

  allPs = mktbook[(mktbook[cntime] <= snaptime) & (mktbook[cltime] >= snaptime) & \
                  ((mktbook[cqty] < 0) | (mktbook[cntime] == snaptime))].reset_index(inplace=False, drop=True)
  askT = allPs[(allPs['side'] == b's') | (allPs['side'] == 's') ].reset_index(inplace=False, drop=True)
  askQ = "[" + ",".join(["[" + str(np.round(askT['price'][i]/1000,3)) + "," + str(np.round(askT['cm_qty'][i]/10,1)) + "]" for i in range(len(askT))]) + "]"
  bidT = allPs[(allPs['side'] == b'b') | (allPs['side'] == 'b') ].reset_index(inplace=False, drop=True)
  bidQ = "[" + ",".join(["[" + str(np.round(bidT['price'][i]/1000,3)) + "," + str(np.round(bidT['cm_qty'][i]/10,1)) + "]" for i in range(len(bidT))]) + "]"
  prow = pd.DataFrame({
    'timestamp': np.asarray([snaptime]),
    'type': np.asarray([b's']),
    'asks': np.asarray([askQ]),
    'bids': np.asarray([bidQ])})
  return(prow);




