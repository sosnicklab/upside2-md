As the first example, the purpose of this example is to help you check whether Upside is installed correctly, and to show how Upside works.

You can use linux shell or python to run Upside respectively. In this folder (example/1.GettingStarted), just run:

./run.sh

or:

source ../../source.sh
python run.py


The script ana.sh shows how to convert the format of the Upside trajectory to vtf format, which can be read by visualization software such as VMD or pymol. Using the "mdtraj_upside" library, the Upside trajectory can also be read as mdtraj format.
just run the command:

./ana.sh

The results are stored in the "results" folder.
