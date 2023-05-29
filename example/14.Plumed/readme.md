# Example folder for using plumed
To use plumed, you'll need to have an installed plumed version, and you need to patch it with upside and re-compile upside

## Set Up Environment
1. install plumed following its [documentation](https://www.plumed.org/doc-v2.8/user-doc/html/_installation.html)
2. patch plumed with upside: 
    1. go to the `$UPSIDE_HOME/src` directory and type `plumed patch --new upside` to create an empty path file for upside
    2. patch upside with the empty patch`plumed patch --patch --engine upside`
    3. save the modification `plumed patch --save`
    These steps will create three sym-linked files in `$UPSIDE_HOME/src` folder: `Plumed.h`, `Plumed.cmake` and `Plumed.inc`. 
    The `CMakeList.txt` will automatically include `Plumed.cmake` and `Plumed.cpp` if `Plumed.cmake` is detected.
3. re-compile upside: go to `$UPSIDE_HOME/obj` folder and type `make`. (A weird bug related to types will be thrown, unless we do `export PLUMED_TYPESAFE_IGNORE=yes`. A full bug report can be seen here: )
4. run simulations: check `0.run.py` to see usages using `advanced_config.py` to config plumed. You'll also need to write a plumed input file, see the example in `./plumed_input`
5. compare results between upside and plumed: run `compare_plumed_upside.ipynb`


