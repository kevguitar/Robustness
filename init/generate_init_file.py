"""
generate_init_file.py (05/09/2018)

Generates an init file with unit covariance matrix and the simplest possible
inital quantites and stores it under "./new_[datdim].init.txt".
"""
import argparse
import numpy as np

# Reads the data dimension as a parser argument passed in the shell.
parser = argparse.ArgumentParser(description='Generates an initialization file \
                                 of specified data dimension.')
parser.add_argument('datadim', type=str, help='The dimension of data.')
args = parser.parse_args()

easy_covariance = True
easy_g_matrix = True
lower_limit = 0
upper_limit = 10

datdim = int(args.datadim)
pardim = 2
data_mean = 0.0
prior_params = [0, 0]
prior_variance = 0.01

output_file = open('new_%i.init' % datdim, 'w')

string0 = '// This file lists all the initial parameters needed to execute ' \
           + 'all four\noperation modes in robustness_main.cc. The program ' \
           + 'reads from here before\nexecution. \n'

string0 += '\n// This semicolon is required for correct parsing.\n'
string0 += ';\n'

string1 = 'datdim = ' + str(datdim) + ';\n'

# Generates the string of the diagonal data_covar_array.
string1 += '\ndata_covar_array = {\n'
if easy_covariance:
    for i in range(datdim/2):
        string1 += '{'
        for j in range(datdim/2 - 1):
            if j == i:
                string1 += '1,'   
            else:
                string1 += '0,'
        if i == (datdim/2 - 1):        
            string1 += '1}\n};\n'
        else:
            string1 += '0},\n' 
else: # tridiagonal matrix
    for i in range(datdim/2):
        string1 += '{'
        for j in range(datdim/2 - 1):
            if j == i:
                string1 += '1,'
           
            elif (j == i+1) or (j == i-1): 
                string1 += '0.2,'   
            else:
                string1 += '0,'
        if i == (datdim/2 - 1):        
            string1 += '1}\n};\n'  
        elif i == (datdim/2 - 2):
            string1 += '0.2},\n'
        else:
            string1 += '0},\n' 


# Generates the string of the uniform data_mean_array.
string1 += '\ndata_mean_array =\n{'
for i in range(datdim/2 - 1):
    string1 += str(data_mean) + ','
string1 += str(data_mean) + '};\n'

string1 += '\npardim = ' + str(pardim) + ';\n'

# Generates the string of the uniform prior_params_array.
string2 = '\nprior_params_array =\n{'
for i in range(len(prior_params) - 1):
    string2 += str(prior_params[i]) + ','
string2 += str(prior_params[-1]) + '};\n'   

# Generates the string of the diagonal prior_matrix.
string2 += '\nprior_matrix = \n{'
for i in range(pardim):
    string2 += '{'
    for j in range(pardim - 1):
        if j == i:
            string2 += str(prior_variance) + ','
        else:
            string2 += '0,'
    if i == (pardim - 1):        
        string2 += str(prior_variance) + '} };\n'
    else:
        string2 += '0},\n'

# Generates the string of the precision matrix. It contains ones everywhere 
# except on the diagonals, where it puts zeros.
string2 += '\ng_matrix_array = \n{'
for i in range(pardim):
    string2 += '{'
    if easy_g_matrix:
        for j in range(datdim/2 - 1):
            if j == i:
                string2 += '0,'
            else:
                string2 += '1,'
        if i == (pardim - 1):   
            string2 += '1 } };\n'
        else:
            string2 += '1},\n'

    else: # non-trivial g_matrix
        if i == (pardim - 1):  
            for j in range(datdim/2 - 1):
                string2 += str(int(np.log(j + 1))) + ','  
            string2 += str(int(np.log(datdim/2))) + ' } };\n'
        else:
            for j in range(datdim/2 - 1):
                string2 += str(int(np.log(datdim/2 - j))) + ',' 
            string2 += '0},\n'

string3 = '\nlower_limit = ' + str(lower_limit) + ';\n'
string3 += 'upper_limit = ' + str(upper_limit) + ';\n'
string3 += 'num_steps_analyt_distribution = 1e4;\n'
string3 += 'num_sampled_data = 1e5\n'
string3 += '\n// No semicolon here'

output_file.write(string0 + string1 + string2 + string3)        
