import pandas as pd; import pyarrow as pa;
import numpy as np;
try :
  print(" pyrus_kaiko first try to from import");
  from pyrus_kaiko import pyrus_kaiko;
  pk = pyrus_kaiko;
except :
  print("--- Init -- unsuccessful import pyrus_kaiko from pyrus_kaiko");
  import pyrus_kaiko;
  pk = pyrus_kaiko;

if "asjoin_py" in pyrus_kaiko :
  asjoin_py = pyrus_kaiko.asjoin_py;
else :
  try :
    import pyrus_marketbook;
    asjoin_py = pyrus_marketbook.asjoin_py;
  except :
    print("pyrus_kaiko..__init__.py could not receive asjoin.py");

# VL = Verbosity Level
import typing
from typing import TypeAlias;

VerbosityLevel: TypeAlias = int;

print("------------------- pyrus_kaiko::__init__.py  running ---------- ---");
print("-- dir(pyrus_kaiko) is [" + ",".join(dir(pyrus_kaiko)) + "]");
print("--- Thank you for running init.");
print("-- dir(pk) is [" + ",".join(dir(pk)) + "]");
print("-------------------------------------------------------------");


if __name__ == "__main__" :
  print("------------------------------------------------------------------------------------");
  import pandas as pd; import pyarrow as pa; import numpy as np;
  print("initiate pyrus_kaiko inside main");
  import pyrus_kaiko;  
  try :
    from pyrus_kaiko import pk;
  except :
    import pyrus_kaiko as pk;
  print("--- dir pyrus_kaiko in main: = [" + ",".join(dir(pyrus_kaiko)) + "]");
  print("--- dir pk in main: = [" + ",".join(dir(pk)) + "]");
  print("-- wut?");
