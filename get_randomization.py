#!/usr/bin/env python

# Script to return randomization scheme for each block of 8 participants
# written by Courtney Gallen, July 2017
# updated by C Gallen, December 2017

#-----------------------------------------------------------------------------
# IMPORTS
#-----------------------------------------------------------------------------
import os, sys
from os.path import join as pjoin
from random import shuffle
import pandas as pd
import numpy as np
from datetime import datetime

#-----------------------------------------------------------------------------
# FUNCTIONS
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Main Code
#-----------------------------------------------------------------------------
def main(argv = sys.argv):

    #-----
    # EXAMPLE 1 (add blank rows for future randomization):
    # rand_type = extend_file
    # n_engage = 4
    # n_control = 4
    # NB: n_engage + n_control must equal 8
    #--
    # EXAMPLE 2 (randomize ppt id to the next blank row):
    # rand_type = rand_ppt
    # ppt_id = 600
    #-----
    
    #-----
    # user input to dictate rand_type
    #-----
    rand_type = argv[1] # extend_file or rand_ppt
    
    #-----
    # GLOBAL variables
    #-----
    N_TOTAL = 8
    ID_NUMS = '600'
    OUT_DIR = pjoin('./randomization')
    FNAME = pjoin(OUT_DIR, 'engage_randomization.csv')
    GROUP_HDR = 'group'
    PPT_HDR = 'ppt_id'
    RAND_HDR = 'rand_type'
    
    #-----
    # Randomize ppt to existing rows ('rand_ppt')
    #-----
    if rand_type == 'rand_ppt': # randomize ppt id to the next open row
        
        ppt_id = int(argv[2])
        # make sure id matches format
        if str(ppt_id)[0] != ID_NUMS[0]:
            print '---------------------------------------------'
            print 'ERROR: ppt id must begin with %s' %ID_NUMS[0]
            sys.exit()
        if len(str(ppt_id)) != len(ID_NUMS):
            print '---------------------------------------------'
            print 'ERROR: ppt id must be %s digits' %len(ID_NUMS)
            sys.exit()

        # if the rand file doesn't exist, error
        if not os.path.exists(FNAME):
            print '---------------------------------------------'
            print 'ERROR: randomization file doesnt exist'
            print 'run script with extend_file option first'
            sys.exit()

        # otherwise, load the file
        else:
            df = pd.read_csv(FNAME)

            # if ppt id is already in this file, error
            if ppt_id in np.array(df[PPT_HDR]):
                print '---------------------------------------------'
                print 'ERROR: %d already randomized' %ppt_id
                sys.exit()
           
            # if no blank lines left, error
            if not df[PPT_HDR].isnull().any():
                print '---------------------------------------------'
                print 'ERROR: no empty rows to randomize'
                print 'run script with extend_file option'
                sys.exit()

            # otherwise, find next blank line
            else:
                next_idx = df[PPT_HDR].index[df[PPT_HDR].apply(np.isnan)][0]
                # get group for this index
                group = df.loc[next_idx, GROUP_HDR]
                # add ppt id and rand_type to this index
                df.loc[next_idx, PPT_HDR] = ppt_id
                df.loc[next_idx, RAND_HDR] = rand_type
                
    elif rand_type == 'extend_file': # add more blank rows for future randomization
        
        n_engage = int(argv[2]) # number of engage ppts needed
        n_control = int(argv[3]) # number of control ppts needed
    
        # make directory if needed
        if not os.path.exists(OUT_DIR):
            os.makedirs(OUT_DIR)

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

        # add this info to a new df
        df_new = pd.DataFrame()
        # add list to pandas dataframe
        df_new[GROUP_HDR] = list_comb
        # add current date and time info to df
        df_new['datetime'] = str(datetime.now())
        # add field for ppt_id for later randomization
        df_new[PPT_HDR] = np.nan
        df_new[RAND_HDR] = np.nan
        
        # if a rand file already exists, load it
        if os.path.exists(FNAME):
            df_old = pd.read_csv(FNAME)
            
            # if there are still NANs in ppt ID, print a warning
            if df_old[PPT_HDR].isnull().any():
                print '---------------------------------------------'
                print 'WARNING: there are still empty rows not randomized'
                print 'only proceed if you need to manually randomize this ppt'
            
            # merge old and new
            df = pd.concat([df_old, df_new], ignore_index=True)
            
        # otherwise, set df_new to df (saved file)
        else:     
            df = df_new    
    
    # save df as a csv in out_dir
    df.to_csv(FNAME, index=False)
    # if rand_ppt, print result
    if rand_type == 'rand_ppt':
         print '---------------------------------------------'
         print 'ppt = %d, group = %s' %(ppt_id, group)
    1/0
                 
if __name__ == '__main__':
    main()
