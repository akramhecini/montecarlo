

Author: 

akram.hecini@gmail.com


Create a program using the method described in an article. This method makes it possible to calculate the folding of a protein by means of a Monte Carlo algorithm . An arbitrary protein sequence can be submitted to the program in order to calculate its folding by this model.

This work is based on : A replica exchange Monte Carlo algorithm for protein folding in the HP model by Thachuk C, Shmygelska A, Hoos HH. 


Environment :

Python 3.7.3

Python modules used :

    numpy
    random
    sys
    argparse
    
    
Directories :

data : the sequence used in the implimented algorithm 

src : python code

docs : the report 


Program :

Only one script is meant to be excuted, it needs 2 arguments : 

1st argument is the fle of the amino acid sequence
2nd argument is the nsteps 

eg :  python3 montecarlo2.py seq.txt 1000
