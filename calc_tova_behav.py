#!/usr/bin/env python

# Script to calculate behavioral variables from TOVA .log Presentation files
# written by Courtney Gallen, July 2017
# usage: calc_tova_behav.py data_dir data_fname rt_outliers
# NB:
# participant ID should be an integer
# assumes sustained is run in first half, followed by impulsive

#-----------------------------------------------------------------------------
# IMPORTS
#-----------------------------------------------------------------------------
import os,string,sys
from os.path import join as pjoin
import pandas as pd
import re
import numpy as np
from scipy.stats import norm
Z = norm.ppf
import utils
reload(utils)

#-----------------------------------------------------------------------------
# FUNCTIONS
#-----------------------------------------------------------------------------
def get_behav_values(df, var, trial_type):
    '''
    Calculate behavioral variables from the pandas data frame
    '''

    # setup data and hdr depending on trial_type
    hdr = '%s_%s' %(var, trial_type)
    if trial_type == 'total':
        data = df # if total, just select all data
    elif trial_type == 'sustained' or trial_type == 'impulsive':
        data = df[df.loc[:, 'Type'] == trial_type] # otherwise, select the sub-data
    else:
        raise NotImplementedError('%s not a valid trial_type' %trial_type)

    # calculate some features of the data
    total_targets = float(np.sum(data['Code'] == '1')) # N targets
    total_nontargets = float(np.sum(data['Code'] == '2')) # N targets
    n_corr_rt = float(data['CorrectRT'].count())   # N targets correctly responded to
    n_commissions = np.sum(data['CommissionsRT'] > 0) # N non-targets responded to
    n_omissions = np.sum(data['Omissions']) # N omissions
    n_ant_targ = np.sum(data['Anticipatory_Target']) # N anticpatory responses to targets
    n_ant_nontarg = np.sum(data['Anticipatory_NonTarget']) # N anticpatory responses to non-targets

    # calculate rate info used in multiple variables
    omis_rate = (n_omissions) / (total_targets - n_ant_targ) * 100
    hit_rate = 1 - (omis_rate)/100
    comis_rate = (n_commissions) / (total_nontargets - n_ant_nontarg) * 100
    fa_rate = comis_rate/100
    
    # calculate behavior
    if var == 'RTMean':
        behav = data['CorrectRT'].mean()/10

    elif var == 'RTStdev':
        behav = data['CorrectRT'].std()/10

    elif var == 'CommissionRate':
        behav = comis_rate

    elif var == 'HitRate':
        behav = hit_rate
        
    elif var == 'FAR':
        behav = (n_commissions - n_ant_nontarg) / (total_nontargets - n_ant_nontarg)

    elif var == 'OmissionRate':
        behav = omis_rate

    elif var == 'Dprime':
        behav = utils.calc_dprime(hit_rate, fa_rate)

    elif var == 'PostCommissionRT':
        1/0
        
    elif var == 'AnticipatoryRate':
        behav = (n_ant_targ + n_ant_nontarg) / (total_targets + total_nontargets) * 100
        
    elif var == 'MultipleResponseRate':
        1/0

    return behav, hdr


def get_demographics(age, sex, norm_table):
    '''Function to get demographics within ranges of norm tables
    '''
    
    # age
    # first see if age is in table as is
    if age in norm_table['age'].unique():
        age_norm = age
    # if not (20-80s), round down to that decade
    elif age not in norm_table['age'].unique() and len(str(age)) == 2:
        age_norm =  int(str(age)[0] + '0')
        # make sure this age is in the table
        if age_norm not in norm_table['age']:
            sys.exit('no demographic data for age %s' %age)
    else:
         raise NotImplementedError('no demographic data for age %s' %age)

    # sex
    sex_norm = []
    # loop through the unique sex options in the table to find a match
    for sex_str in norm_table['sex'].unique():
        # see if first letter of normative sex matches participant sex
        if sex_str[0] == sex:
            sex_norm.append(sex_str)
    # make sure something wonky hasn't happened and you only have one value in sex_norm
    if len(sex_norm) != 1:
        sys.exit('too many matches for ppt sex %s' %(sex))

    return age_norm, sex_norm[0]


def acs_zscore(value, mean, std):
    '''Function to calculate zscore based on input mean and std
    '''
    return (value - mean) / std
            

#-----------------------------------------------------------------------------
# Main Code
#-----------------------------------------------------------------------------
def main(argv = sys.argv):

    #-----
    # user inputs
    #-----
    data_dir = argv[1] # top level directory where .log lives
    data_fname = argv[2] # name of .log file to load (NB: this must live in {data_dir})
    rt_outliers = argv[3] # clean_outliers (remove trials > 2SD RTs) -or- keep_outliers (do nothing)
    # OPTIONAL: input age and sex for ACS calculation
    if len(argv) > 4:
        age = int(argv[4]) # age of participant
        sex = argv[5] # sex as M or F
    #-----
    # EXAMPLE:
    # data_dir = /Users/courtneygallen/Dropbox/Gazzaley/Data/engage/TOVAprocessing/
    # data_fname = 2573_Pre-tova-v.log
    # rt_outliers = keep_outliers
    # age = 30
    # sex = F
    #-----

    #-----
    # GLOBAL variables
    #-----
    # all variables but ACS, which is calculated after these
    VARS = sorted(['RTMean', 'RTStdev', 'CommissionRate', 'OmissionRate', 'HitRate', 'FAR', 'Dprime',
                   'PostCommissionRT', 'AnticipatoryRate', 'MultipleResponseRate'])
    TRIAL_TYPES = ['total', 'sustained', 'impulsive']
    OUTLIER_THRESH = '2SD'
    DEFAULT_N_TRIALS = 508
    
    # load the data into a pandas dataframe, ignoring the first 2 rows
    df = pd.read_csv(pjoin(data_dir, data_fname), skiprows = 2, delimiter = '\t', dtype={'Code': str})
    sub, sess = utils.get_sub_info(df, 'tova')
    print '------------------------------------------------'
    print 'working on SUB: %s, SESS: %s' %(sub, sess)
    
    # get trial number info
    n_trials = df.loc[:, 'Trial'].max() # total number of trials
    n_subtrials = n_trials/2 # number of trials per sustained/impulsive
    n_qtrials = n_subtrials/2 # nuumber of trials per quarter
    # quit if n_trials not what expected
    #if n_trials != DEFAULT_N_TRIALS:
    #    sys.exit('%d trials different from %d expected' %(n_trials, DEFAULT_N_TRIALS))
        
    
    ### Add columns to df for behavioral calculations ###
    #-----
    # Trial: this is the trial type (sustained vs impulsive)
    #-----
    # get trial numbers corresponding to sustained/impulsive
    sus_trials = np.array((range(1, n_subtrials + 1))) # first half is sustained, 1-indexed
    imp_trials = np.array((range(1, n_subtrials + 1) + n_subtrials)) # second half is impulsive, 1-indexed
    
    # get rows related to sustained/impulsive (pictures, responses, fixation)
    sus_rows = utils.return_indices_of_a(df.loc[:, 'Trial'], sus_trials)
    imp_rows = utils.return_indices_of_a(df.loc[:, 'Trial'], imp_trials)
    
    # label the rows sustained/impulsive in column 'Type'
    df.loc[sus_rows, 'Type'] = 'sustained' 
    df.loc[imp_rows, 'Type'] = 'impulsive'

    #-----
    # Get rows/indices/rt info for code = 1 and code = 2 stimuli
    #-----
    # picture 1 (target)
    code1_rows = df.loc[:, 'Code'] == '1' # rows of code = 1
    code1_idx = df[code1_rows].index # indices of code1 rows
    code1_rtidx = [x+1 for x in code1_idx.tolist()] # RTs are in row after each picture
    code1_rts = df.loc[code1_rtidx, 'TTime']
    # picture 2 (non-target)
    code2_rows = df.loc[:, 'Code'] == '2' # rows of code = 2
    code2_idx = df[code2_rows].index # indices of code2 rows
    code2_rtidx = [x+1 for x in code2_idx.tolist()] # RTs are in row after each picture
    code2_rts = df.loc[code2_rtidx, 'TTime']
    1/0
    #-----
    # CorrectRT: these are response times for Pictures with Code = 1 stimulus
    # (analogous to Ttime in Erwin's calculations)
    #-----
    # add rts to df, in row where picture was presented
    df.loc[code1_idx, 'CorrectRT'] = np.array(code1_rts)

    #-----
    # CommissionsRT: these are RTs of responses for Pictures with Code = 2 stimulus
    #-----
    # loop through impulsive/sustained trial types
    # this is so the start of the impulsive trials aren't influenced by end of sustained
    # where the next row is an RT for the fixation cross between trials
    for tidx, ttype in enumerate(df['Type'].unique()):

        # get the inds just for this trial type
        ttype_inds = df.loc[:, 'Type'] == ttype
        
        # select data for this trial type and code 2
        code2_rows_pertrial = df.loc[ttype_inds, 'Code'] == '2'
        code2_idx_pertrial = df.loc[ttype_inds & code2_rows_pertrial].index # indices of code2 rows

        # get rt data
        code2_rtidx_pertrial = [x+1 for x in code2_idx_pertrial.tolist()] # RTs are in row after each picture
        code2_rts_pertrial = np.array(df.loc[code2_rtidx_pertrial, 'TTime'])
        
        # remove rt from start of next trial type if present
        if np.array(df[ttype_inds].index)[-1] != np.array(code2_rtidx_pertrial)[-1]:
            code2_rts_pertrial[-1] = 0
    
        df.loc[ttype_inds & code2_rows_pertrial, 'CommissionsRT'] = code2_rts_pertrial
    
    #-----
    # Omissions: these are booleans of yes/no RT for Pictures with Code = 1 stimulus
    #-----
    # find where RTs for code 1 are 0
    code1_omit = code1_rts == 0
    # add this to df in row where picture was presented
    df.loc[code1_idx, 'Omissions'] = np.array(code1_omit)

    #-----
    # Anticipatory_Target, Anticipatory_NonTarget: these are booleans of RTs < 150ms
    # NB: RTs in this file are *10ms
    #-----
    # code 1 (targets)
    # find where RTs for code 1 are les than 1500
    # add to row where pic was presented
    code1_ant = (code1_rts > 0) & (code1_rts < 1500)
    df.loc[code1_idx, 'Anticipatory_Target'] = np.array(code1_ant)
    # code 2 (targets)
    # find where RTs for code 2 are les than 1500
    # add to row where pic was presented
    code2_ant = (code2_rts > 0) & (code2_rts < 1500)
    df.loc[code2_idx, 'Anticipatory_NonTarget'] = np.array(code2_ant)

    
    ### Clean up CorrectRT ###
    
    # set RTs < 150ms to nan across all trials
    # cleans up: anticipatory RTs (RTs between 0 and 1500)
    # cleans up: omission error RTs (RTs == 0)
    all_ant = df.loc[:, 'CorrectRT'] < 1500
    # set these values to nans
    df.loc[all_ant, 'CorrectRT'] = np.nan
    
    # outliers: if clean_outliers, remove RTs > 2SD (separately for each trial type)
    if rt_outliers == 'clean_outliers' and OUTLIER_THRESH != '2SD':
        raise NotImplementedError('%s criteria for outliers not implemented')
    
    elif rt_outliers == 'clean_outliers' and OUTLIER_THRESH == '2SD':
        
        # loop through impulsive/sustained trial types
        for tidx, ttype in enumerate(df['Type'].unique()):

            # get the inds just for this trial type
            ttype_inds = df.loc[:, 'Type'] == ttype
            
            # get the upper limit for 2SD
            mean = df.loc[ttype_inds, 'CorrectRT'].mean()
            std = df.loc[ttype_inds, 'CorrectRT'].std()
            out_thresh = mean + (2*std)

            # set values above threshold to nan
            outliers = df.loc[ttype_inds, 'CorrectRT'] > out_thresh
            df.loc[ttype_inds & outliers, 'CorrectRT'] = np.nan
            
    
    ### Calculate performance variables ###
    
    # setup csv to add behavioral variables
    df_out = pd.DataFrame({'subject': sub, 'session': sess}, index = np.arange(1))
    
    # loop through variables
    for vidx, var in enumerate(VARS):

        # loop through trial types
        for tidx, ttype in enumerate(TRIAL_TYPES):

            # get behavioral variable and hdr
            val, hdr = get_behav_values(df, var, trial_type = ttype)
            # add to df_out
            df_out[hdr] = val

    
    # now calculate ACS variable using the already calculated variables
    # NB: this only runs if age and sex were entered as inputs
    try:
    
        # load normative data
        # NB: these need to live in the same directory as this script
        rtnorm_all = pd.read_csv('rt_normative_data.csv') # mean RT
        dpnorm_all = pd.read_csv('dprime_normative_data.csv') # dprime
        varnorm_all = pd.read_csv('variability_normative_data.csv') # RT variability (stdev)
        
        # get demographic info that matches table (age, sex)
        age_norm, sex_norm = get_demographics(age, sex, rtnorm_all)

        # get norm data for this demographic
        rtnorm_idx = (rtnorm_all.loc[:, 'age'] == age_norm) & (rtnorm_all.loc[:, 'sex'] == sex_norm)
        dpnorm_idx = (dpnorm_all.loc[:, 'age'] == age_norm) & (dpnorm_all.loc[:, 'sex'] == sex_norm)
        varnorm_idx = (varnorm_all.loc[:, 'age'] == age_norm) & (varnorm_all.loc[:, 'sex'] == sex_norm)
        
        # calculate z-scores
        # 1: rt - half 1 (sustained)
        rt_z = acs_zscore(float(df_out['RTMean_sustained']), float(rtnorm_all['half1_mean'].loc[rtnorm_idx]),
                          float(rtnorm_all['half1_std'].loc[rtnorm_idx]))
        # 2: dprime - half 2 (impulsive)
        dp_z = acs_zscore(float(df_out['Dprime_impulsive']), float(dpnorm_all['half2_mean'].loc[dpnorm_idx]),
                          float(dpnorm_all['half2_std'].loc[dpnorm_idx]))
        # 3: rt std (total)
        var_z = acs_zscore(float(df_out['RTStdev_total']), float(varnorm_all['total_mean'].loc[varnorm_idx]),
                           float(varnorm_all['total_std'].loc[varnorm_idx]))

        # calculate acs score
        acs = (rt_z*-1) + dp_z + (var_z*-1) + 1.80

        # add acs and demographic info to df_out
        df_out['ACS'] = acs
        df_out['Age'] = age
        df_out['Sex'] = sex

    except NameError:

        print '------------------------------------------------'
        print 'age & sex not entered, NOT calculating ACS value'


    # save the file
    # set up names for output directory and csv file
    if rt_outliers == 'keep_outliers':
        out_dir_name = 'outliers-kept'
        fname = '%s%s_%s.csv' %(sub, sess, out_dir_name)
    elif rt_outliers == 'clean_outliers':
        out_dir_name = '%soutliers-removed' %OUTLIER_THRESH
        fname = '%s%s_%s.csv' %(sub, sess, out_dir_name)

    # create output dir if it doesn't exist
    out_dir = pjoin(data_dir, out_dir_name)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # save csv
    df_out.to_csv(pjoin(out_dir, fname), index=False)
    1/0
        
    

if __name__ == '__main__':
    main()
