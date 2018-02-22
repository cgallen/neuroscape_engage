#!/usr/bin/env python

"""Generic utilities that may be needed by the other modules.
"""
# Written by Courtney Gallen, February 2018

#-----------------------------------------------------------------------------
# IMPORTS
#-----------------------------------------------------------------------------
import os,string,sys
from os.path import join as pjoin
import re
import numpy as np
from scipy.stats import norm
Z = norm.ppf

#-----------------------------------------------------------------------------
# FUNCTIONS
#-----------------------------------------------------------------------------
def get_sub_info(df, task):
    '''
    Get subject and session info from datafile name
    '''
    # get the subject string input by experimenter
    if task == 'tova':
        sub_str = df['Subject'][0]
    elif task == 'filter':
        sub_str = df
    else:
        raise NotImplementedError('%s sub naming not implemented' %task)
    
    # get sub num (integers in string)
    try:
        sub_num = map(int, re.findall(r'\d+', sub_str))
        # make sure something wonky hasn't happened and you only have one sess
        if len(sub_num) != 1:
            print 'WARNING: multiple sub matches, %s' %(sub_num)
        sub_num = sub_num[0]
            
        # get sess info
        pre_list = ['pre', 'Pre', 'PRE']
        post_list = ['post', 'Post', 'POST']
        fu_list = ['fu', 'FU', 'followup', 'FOLLOWUP']
        sess = []
        if any(x in sub_str for x in pre_list):
            sess.append('pre')
        if any(x in sub_str for x in post_list):
            sess.append('post')
        if any(x in sub_str for x in fu_list):
            sess.append('fu')
        # if no session info, set to open quotes
        if len(sess) == 0:
            sess = ''
        # make sure something wonky hasn't happened and you only have one sess
        elif len(sess) > 1:
            sys.exit('too many matches for sess %s' %(sess))
        else:
            sess = sess[0]

    except TypeError:
        # if the subject string is just an integer, make that sub_num
        try:
            sub_str = int(sub_str)
            sub_num = sub_str
            sess = ''
        except ValueError:
            sys.exit('ppt-id formatting not implemented, %s' %(sub_str))
            
    return sub_num, sess


def return_indices_of_a(a, b):
    '''
    find indices of a that match b
    '''
    b_set = set(b)
    return [i for i, v in enumerate(a) if v in b_set]


def calc_dprime(hit_rate, fa_rate):
    '''calculate dprime, adapted from:
    http://lindeloev.net/calculating-d-in-python-and-php/
    with some modifications specified by the TOVA manual
    '''
    # Calculate hitrate and avoid d' infinity
    hitRate = hit_rate
    if hitRate == 1: hitRate = 0.99999
    if hitRate == 0: hitRate = 0.00001
 
    # Calculate false alarm rate and avoid d' infinity
    faRate = fa_rate
    if faRate == 1: faRate = 0.99999
    if faRate == 0: faRate = 0.00001
 
    # Return d'
    dprime = Z(hitRate) - Z(faRate)
    
    return dprime

