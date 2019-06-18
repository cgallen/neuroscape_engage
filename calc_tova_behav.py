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
import ExGUtils.pyexg as exg
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
    n_multiple_resp = np.sum(data['MultipleResponse']) # N multiple target responses (correct trials)
    
    # calculate rate info used in multiple variables
    omis_rate = (n_omissions) / (total_targets - n_ant_targ) * 100
    hit_rate = 1 - (omis_rate)/100
    comis_rate = (n_commissions) / (total_nontargets - n_ant_nontarg) * 100
    fa_rate = comis_rate/100

    # calculate behavior
    if var == 'rt_mean':
        behav = data['CorrectRT'].mean()/10

    elif var == 'rt_stdev':
        behav = data['CorrectRT'].std(ddof=0)/10

    elif var == 'commission_rate':
        behav = comis_rate

    elif var == 'hit_rate':
        behav = hit_rate
        
    elif var == 'far':
        behav = fa_rate

    elif var == 'omission_rate':
        behav = omis_rate

    elif var == 'dprime':
        behav = utils.calc_dprime(hit_rate, fa_rate)

    elif var == 'postcommission_rt':
        behav = data['PostCommissionsRT'].mean()/10
        
    elif var == 'anticipatory_rate':
        behav = (n_ant_targ + n_ant_nontarg) / (total_targets + total_nontargets) * 100
        
    elif var == 'multipleresponse_rate':
        behav = n_multiple_resp / total_targets * 100
        
    elif var == 'exgauss':
        # get exgauss variables on valid RT data
        #try:
        mu, sigma, tau = exg.maxLKHD(np.array(data['CorrectRT'].dropna()/10))
        exg_func = 'max_lkhd'
        #except ZeroDivisionError:
        #    [x, y] = exg.histogram(np.array(data['CorrectRT'].dropna()/10))
        #    mu, sigma, tau = exg.minSQR(x,y)
        #    exg_func = 'min_sqr'
        #    1/0
        behav = [mu, sigma, tau]
        hdr = ['mu_%s' %trial_type, 'sigma_%s' %trial_type, 'tau_%s' %trial_type]

        
    return behav, hdr


def get_code_rtidx(df, code_idx):
    '''Function to get associated RT rows after a code presentation, where Code = 100
    and the Trial number matches associated RT row
    '''
    code_rtidx = []
    omit_rtidx = []
    #loop through code_idx, see if there was a response
    for row in code_idx:
        #only run this analysis if it's not the last row in df
        last_row = len(df) - 1
        if row != last_row:
            # get trial number for picture
            trial_num = df.loc[row, 'Trial']
            # get trial num for next trial (where RTs are logged)
            next_row = row + 1
            # see if there was a Response in next_row, with same trial_num
            next_trial = df.loc[next_row, 'Trial']
            next_resp = df.loc[next_row, 'Code']
            next_event = df.loc[next_row, 'Event Type']

            # if there was a response, add it
            if next_resp == '100':
                code_rtidx.append(next_row)
                if next_trial != trial_num:
                    # add to omit_rtidx, this is a response for a separate trial
                    # usually the RT to start the 2nd block
                    omit_rtidx.append(next_row)

            # if there was no response, add it to be counted as a miss later
            elif next_event == 'Picture':
                code_rtidx.append(next_row)
                #1/0

        #if it is the last row, add this one to an omission
        #elif row == last_row:
        #    code
        #    1/0
    
    return code_rtidx, omit_rtidx

    
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
    # exgauss is a catch all for all three variables: mu, sigma, tau
    VARS = ['rt_mean', 'rt_stdev', 'commission_rate', 'omission_rate', 'hit_rate', 'far', 'dprime',
            'postcommission_rt', 'anticipatory_rate', 'multipleresponse_rate', 'exgauss']
    TRIAL_TYPES = ['sustained', 'impulsive', 'total']
    OUTLIER_THRESH = '2SD'
    DEFAULT_N_TRIALS = 500
    FIXATION_TRIALS = [1, 500*(1./4)+2, 500*(2./4)+3, 500*(3./4)+4] # ignore these RTs
    
    # load the data into a pandas dataframe, ignoring the first 2 rows
    df = pd.read_csv(pjoin(data_dir, data_fname), skiprows = 2, delimiter = '\t', dtype={'Code': str})
    sub, sess = utils.get_sub_info(df, 'tova')
    print '------------------------------------------------'
    print 'working on SUB: %s, SESS: %s' %(sub, sess)
    
    # get trial number info
    code1_rows = df.loc[:, 'Code'] == '1' # rows of code = 1 (target)
    code2_rows = df.loc[:, 'Code'] == '2' # rows of code = 2 (non-target)
    n_trials = sum(code1_rows) + sum(code2_rows) # total number of trials
    n_subtrials = n_trials/2 # number of trials per sustained/impulsive
    n_qtrials = n_subtrials/2 # nuumber of trials per quarter
    # quit if n_trials not what expected
    if n_trials != DEFAULT_N_TRIALS:
        sys.exit('%d trials different from %d expected' %(n_trials, DEFAULT_N_TRIALS))
        
    ### Add columns to df for behavioral calculations ###
    #-----
    # Trial: this is the trial type (sustained vs impulsive)
    #-----
    # get trial numbers corresponding to sustained/impulsive
    ## QUESTION: what is the best way to do this?
    #sus_trials = np.array((range(1, n_subtrials + 1))) # first half is sustained, 1-indexed
    #imp_trials = np.array((range(1, n_subtrials + 1) + n_subtrials)) # second half is impulsive, 1-indexed
    sus_trials = np.array((range(1, (df['Trial'].max()/2) + 1))) # first half is sustained, 1-indexed
    imp_trials = np.array((range(1, (df['Trial'].max()/2) + 1) + (df['Trial'].max()/2))) # second half is impulsive, 1-indexed
    
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
    code1_idx = df[code1_rows].index # indices of code1 rows
    code1_rtidx, code1_zeroidx = get_code_rtidx(df, code1_idx)
    #get rts, apply corrections if needed
    code1_rts = df.loc[code1_rtidx, 'TTime']
    code1_rts[code1_zeroidx] = 0
    # picture 2 (non-target)
    code2_idx = df[code2_rows].index # indices of code2 rows
    code2_rtidx, code2_zeroidx = get_code_rtidx(df, code2_idx)
    code2_rts = df.loc[code2_rtidx, 'TTime']
    code2_rts[code2_zeroidx] = 0
    
    #-----
    # CorrectRT: these are response times for Pictures with Code = 1 stimulus
    # (analogous to Ttime in Erwin's calculations)
    #-----
    # add rts to df, in row where picture was presented
    df.loc[code1_idx, 'CorrectRT'] = np.array(code1_rts)
    
    #-----
    # CommissionsRT: these are RTs of responses for Pictures with Code = 2 stimulus
    #-----
    # add rts to df, in row where picture was presented
    # first make sure not a corner case where last trial was Code = 2 and no resp
    last_row_num = df.shape[0]-1
    last_trial_type = df.loc[last_row_num]['Code']
    #if trial code = 2, append a 0 onto code2_rts so the indexing is the correct shape
    if last_trial_type == '2':
        code2_rts_new = np.append(code2_rts, 0)
        df.loc[code2_idx, 'CommissionsRT'] = np.array(code2_rts_new)
    else:
        df.loc[code2_idx, 'CommissionsRT'] = np.array(code2_rts)
    
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
    # first make sure not a corner case where last trial was Code = 2 and no resp
    #if trial code = 2, append a 0 onto code2_rts so the indexing is the correct shape
    if last_trial_type == '2':
        code2_rts_new = np.append(code2_rts, 0)
        code2_ant = (code2_rts_new > 0) & (code2_rts_new < 1500)
    else:
        code2_ant = (code2_rts > 0) & (code2_rts < 1500)
    df.loc[code2_idx, 'Anticipatory_NonTarget'] = np.array(code2_ant)
    
    ### Clean up CorrectRT ###
    # set RTs < 150ms to nan across all trials
    # cleans up: anticipatory RTs (RTs between 0 and 1500)
    # cleans up: omission error RTs (RTs == 0)
    all_ant = df.loc[:, 'CorrectRT'] < 1500
    # set these values to nans
    df.loc[all_ant, 'CorrectRT'] = np.nan
    ## do the same to CommissionRT ##
    all_com_ant = df.loc[:, 'CommissionsRT'] < 1500
    df.loc[all_com_ant, 'CommissionsRT'] = np.nan
    
    #-----
    # 'MultipleResponse': booleans of whether an RT exists after code1_rtidx
    # QUESTION: if removing outliers, do this after?
    #-----
    df.loc[:, 'MultipleResponse'] = np.nan
    # loop through code1_rtidx
    for rt_row in code1_rtidx:
        # get correct RT info
        corr_rt_row = rt_row - 1 # where the 'CorrectRT' lives
        corr_rt = df.loc[corr_rt_row, 'CorrectRT']

        # get post rt row info, if a response - add it
        post_rt_row = rt_row + 1 # row after the RT
        
        # only run if this is a correct RT trial and not the end of df
        if not np.isnan(corr_rt) and post_rt_row < df.shape[0]:
            # if a response - add it
            post_code = df.loc[post_rt_row, 'Code']
            post_trial = df.loc[post_rt_row, 'Trial']
            if post_code == '100' and post_trial not in FIXATION_TRIALS:
                df.loc[corr_rt_row, 'MultipleResponse'] = True

    #-----
    # PostCommissionsRT: these are correct target RTs after a commission
    #-----
    df.loc[:, 'PostCommissionsRT'] = np.nan
    # get rows and indices of commission errors
    comm_rows = df.loc[:, 'CommissionsRT'] > 0
    comm_idx = df[comm_rows].index
    # get rows and indices of correct RT
    corrt_rows = df.loc[:, 'CorrectRT'] > 0
    corrt_idx = df[corrt_rows].index
    
    # loop through comm_idx, to find next correct rt
    for cidx, comm_row in enumerate(comm_idx):
        comm_trial = df.loc[comm_row, 'Trial']
        
        # don't run if it's the last trial 
        if comm_trial != imp_trials[-1]:
            # index, trial of next correct rt
            next_corrt_row = corrt_idx[np.where(corrt_idx > comm_row)[0][0]]
            next_corrt_trial = df.loc[next_corrt_row, 'Trial']

            # make sure the commission and next correct RT are in the same task block
            for fidx, fix_trial in enumerate(FIXATION_TRIALS):
                if not comm_trial < next_corrt_trial < fix_trial:
                    # if corrt is before the next commission, add it to the commission row
                    if cidx+1 < len(comm_idx):
                        if next_corrt_row < comm_idx[cidx+1]:
                            post_comm_rt = df.loc[next_corrt_row, 'CorrectRT']
                            df.loc[comm_row, 'PostCommissionsRT'] = post_comm_rt
                            
                    # if it's the last item in comm_row list, just take the next target
                    elif cidx+1 == len(comm_idx):
                        post_comm_rt = df.loc[next_corrt_row, 'CorrectRT']
                        df.loc[comm_row, 'PostCommissionsRT'] = post_comm_rt
                        
        
    ### Clean up Outliers ###
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
            if var != 'exgauss':
                df_out[hdr] = val
            elif var == 'exgauss':
                # loop through hdr and add relevant val
                for idx, exg in enumerate(hdr):
                    df_out[exg] = val[idx]

              
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
        rt_z = acs_zscore(float(df_out['rt_mean_sustained']), float(rtnorm_all['half1_mean'].loc[rtnorm_idx]),
                          float(rtnorm_all['half1_std'].loc[rtnorm_idx]))
        # 2: dprime - half 2 (impulsive)
        dp_z = acs_zscore(float(df_out['dprime_impulsive']), float(dpnorm_all['half2_mean'].loc[dpnorm_idx]),
                          float(dpnorm_all['half2_std'].loc[dpnorm_idx]))
        # 3: rt std (total)
        var_z = acs_zscore(float(df_out['rt_stdev_total']), float(varnorm_all['total_mean'].loc[varnorm_idx]),
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
    1/0
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
