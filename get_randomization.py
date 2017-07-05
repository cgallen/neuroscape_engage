#!/usr/bin/env python

# Script to return randomization scheme for each block of 8 participants
# written by Courtney Gallen, July 2017

#-----------------------------------------------------------------------------
# IMPORTS
#-----------------------------------------------------------------------------
import os, sys
from os.path import join as pjoin
from random import shuffle
import pandas as pd
from datetime import datetime

#-----------------------------------------------------------------------------
# FUNCTIONS
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Main Code
#-----------------------------------------------------------------------------
def main(argv = sys.argv):

    #-----
    # user inputs
    #-----
    n_engage = int(argv[1]) # number of engage ppts needed
    n_control = int(argv[2]) # number of control ppts needed
    #-----
    # EXAMPLE:
    # n_engage = 4
    # n_control = 4
    # NB: n_engage + n_control must equal 8
    #-----

    #-----
    # GLOBAL variables
    #-----
    N_TOTAL = 8

    # check that total number of participants equals 8
    if n_engage + n_control != N_TOTAL:
        sys.exit('total number of ppts input equals %d, not %d' %((n_engage+n_control), N_TOTAL))

    # set up lists of engage and control ppts
    list_engage = ['engage'] * n_engage
    list_control = ['control'] * n_control

    # combine the two lists
    list_comb = list_engage + list_control
    # shuffle the combined list
    shuffle(list_comb)

    # set up directory to save randomized list
    out_dir = pjoin('./randomization')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # add list to pandas dataframe
    df = pd.DataFrame({'group': list_comb})
    # add current date and time info to df
    df['datetime'] = str(datetime.now())

    # save df as a csv in out_dir
    df.to_csv(pjoin(out_dir, 'engage_randomization.csv'), index=False)
    1/0
                 
if __name__ == '__main__':
    main()
