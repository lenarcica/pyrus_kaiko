## pyrus_kaiko
##
## A collection of utilities for ssimulating or reading from Kaiko-style continuous time orderbook
##
##  Alan Lenarcic, 2025/10


# LICENSE:  GNU General Public License, version 2
  -- The following code is a prototype and demonstration of Rust and python capabilities with simulated data.
For use in commercial products, others are recommended to rewrite for your own purposes.

## Package features:
  1. Simulate Continuous Orderbook according to a stochasitic midpoint + clipped simulated order process
  2. Read typicla format Kaiko FOB files
  3. Use the code for futher analyses such as in nbbo calculations.

## Testing Python/Pyo3/Rust ability to unzip Kaiko Gzip files
### A supplemental file reader in RUST.
   Here we simulate gzip files and demonstrate Kaiko has features that can read these files.
   
   Due to the "grouped update" format of FOB files, this text reader algorithm will read the output of a FOB
  file and split into all price/time/side categories to be tracked separately.


#  KAIKO "fob" files (Full-Order-Book)
exist in a strange format

|time  | type ["U"PDATE/"S"NAP]  | Buys Orders     | Sells Orders 
| 0:00 | U                       | [[10.0,1000]]   | [11.0, 2000]
| 0:01 | S                       | [[9.0,,400]...] | [[10.0,1000], ....]
| 0:03 | U                       | [[10.1,900]]    | [[10.1,0],[10.2,140]]

We wish to break down into a side, price, time format: as an example:
| side  | price   | time           | qty 
|    B  |  0.01   | 09:40:00.302   | 100
|    B  |  0.01   | 09:40:03.532   | 200
|    B  |  0.01   | 10:20:04:323   |   0
|    B  |  0.02   | ....
| ...  ... ....
|    S  | 9999.99 | 23:29:59.999   |   0

This format separates all price levels into individual records and notes each change individually.  This is used for
 more performant orderbook depth algorithms.

### Project components:
  1. Python source files (pyrus_kaiko/pyrus_kaiko directoriy.  Includes
    -- __init__.py: main script file, tests to ensure that both Rust functions and internal script modules are loaded.
  2.


A "lib.rs" file, located in src directory, is the Rust file that will be compiled to start this library.

## To Compile with Maturin
   
  Installation of Rust is achieved using Anaconda python environment, python 3.12 or sooner, which would include Rust and pyo3 related packages.

  In windows, llvm must be installed Visual Studio tools.
 
### Compiling this package: 
To install package install pyo3, pyo3 tools, and "maturin" from Python PIP
```

```

Then initiate maturin
```
maturin init

```

And 
```
maturin develop
```

Appears to build the tool.  Note that access to CARGO and many pyo3 and arrow related packages will be required.

### Note to DROPBOX users
  This code is best compiled outside a Dropbox/OneDrive or other synced directory.

  Configure a ".cargo/config.toml" file where you will create a separate compilation location so that the downloaded cargo packages and the 
compiled code does not interfere with normal cloud syncronization.  

# Package contents
  - src: RUST source code
    a. lib.rs: Primary Matruin access code, demonstrating Arrow Parquet and other interaction formats between Rust/Python
      External Facing functions:
       kaiko_fob[]: main function called to extract a GZ zipped FOB file and return ARROW format table to python for further processing.
       kaiko_make_u_fob[]: Take simulated orderbook data and create a fake GZ zipped FOB file
       asjoin_py[]: simple demonstration of As_of_join in rust for pre-sorted columns of timestamp data
       unsorted_asjoin_py[]: Harder dmeonstration for pre-unsorted data
       pyrus_verify_prices[]: Used in simulating fake orderbook data to guarantee that simulations are "clipped" when they are hit by midpoint.
       sip_arrow[]: SIP NBBO from separate BBO streat demonstration
    b. sip_algo.rs/sip_struct.rs: Demonstration of a basic NBBO calculation on SIP format NBBO data (such as Kaiko TOB files that would need to be integrated)
    c. b2v_struct.rs:  Bit-based structures for keeping tabs on each Exchange's ability to keep at top of book (see sip_algo.rs material)
    c. verify_price.rs algo: Main algorithm for "clipping" simulated orders against the simulated midprice to create a large set of simulated orderbook data
