As the first example, the purpose of this example is to help you check whether Upside is installed correctly, and to show how Upside works.

You can use linux shell or python to run Upside respectively. In this folder (example/01.GettingStarted), just run:

./0.run.sh

or:

source ../../source.sh
python 0.run.py


The script 1.ana.sh shows how to convert the format of the Upside trajectory to vtf format, which can be read by visualization software such as VMD or pymol. Using the "mdtraj_upside" library, the Upside trajectory can also be read as mdtraj format.
just run the command:

./1.ana.sh

The results are stored in the "results" folder.
