
import argparse
import numpy as np


def generate_covariance_string(datdim, cov=None, easy_diagonal=False,
                               easy_tridiagonal=False):
    """
    Generates the string of the diagonal data_covar_array.
    """
    string = '\ndata_covar_array = {\n'
    if cov is not None:
        for i in range(int(datdim/2)):
            string += '{'
            for j in range(int(datdim/2 - 1)):
                string += str(cov[i, j]) + ", "
            if i == (int(datdim/2 - 1)):  # last element        
                string += str(cov[i, i]) + "}\n};\n"
            else:  # last element in row
                string += str(cov[i, int(datdim/2)] - 1) + '},\n' 

    elif easy_diagonal:
        for i in range(int(datdim/2)):
            string += '{'
            for j in range(int(datdim/2 - 1)):
                if j == i:
                    string += '1,'   
                else:
                    string += '0,'
            if i == int(datdim/2 - 1):        
                string += '1}\n};\n'
            else:
                string += '0},\n'

    elif easy_tridiagonal:
        for i in range(int(datdim/2)):
            string += '{'
            for j in range(int(datdim/2 - 1)):
                if j == i:
                    string += '1,'
            
                elif (j == i+1) or (j == i-1): 
                    string += '0.2,'   
                else:
                    string += '0,'
            if i == (datdim/2 - 1):        
                string += '1}\n};\n'  
            elif i == (datdim/2 - 2):
                string += '0.2},\n'
            else:
                string += '0},\n' 
    else:
        raise ValueError("You have to either provide a covariance matrix, "
                         "or select 'easy_diagonal' or easy_tridiagonal'.")        

    return string


def generate_data_mean_string(datdim, data_mean):
    """
    Generates the string of the data_mean_array.
    """
    string = '\ndata_mean_array =\n{'
    for i in range(int(datdim/2 - 1)):
        string += str(data_mean[i]) + ','
    string += str(data_mean[-1]) + '};\n'

    return string


def generate_prior_mean_string(prior_params):
    """
    Generates the string of the prior mean values (prior_params_array).
    """
    string = '\nprior_params_array =\n{'
    for i in range(len(prior_params) - 1):
        string += str(prior_params[i]) + ','
    string += str(prior_params[-1]) + '};\n'

    return string


def generate_prior_matrix_string(pardim, prior_mat):
    """
    Generates the string of the prior_matrix.
    """
    string = '\nprior_matrix = \n{'
    for i in range(pardim):
        string += '{'
        for j in range(pardim - 1):
            string += str(prior_mat[i, j]) + ','
        if i == (pardim - 1):  # last element      
            string += str(prior_mat[i, i]) + "}\n};\n"
        else:  # last element in line
            string += str(prior_mat[i, pardim - 1]) + '},\n'
    return string


def generate_design_matrix_string(pardim, datdim, g_matrix=None,
                                  easy_g_matrix=False):
    """
    Generates the string of the design matrix. It contains ones everywhere 
    except on the diagonals, where it puts zeros.
    """
    string = '\ng_matrix_array = \n{'
    for i in range(pardim):
        string += '{'
        if g_matrix is not None:
            for j in range(int(datdim/2) - 1):
                string += str(g_matrix[i, j]) + ','
            if i == (pardim - 1):  # last element
                string += str(g_matrix[i, int(datdim/2) - 1]) + '}\n};\n'
            else:  # end of line
                print("g_matrix", g_matrix.shape,
                      "datdim", datdim, 
                      "idx", i, int(datdim/2) - 1)
                string += str(g_matrix[i, int(datdim/2) - 1]) + '},\n'

        elif easy_g_matrix:
            for j in range(int(datdim/2) - 1):
                if j == i:
                    string += '0,'
                else:
                    string += '1,'
            if i == (pardim - 1):  # last element
                string += '1 }\n};\n'
            else:  # end of line
                string += '1},\n'

        else:
            raise ValueError("You have to either provide a g_matrix, or "
                             "set easy_g_matrix to true.")
    return string


def main(config):
    """
    Generates an init file with unit covariance matrix and the simplest possible
    inital quantites and stores it under "./new_[datdim].init.txt".
    """
    import yaml
    with open(config.globals, "r") as f:
        args = yaml.load(f, Loader=yaml.SafeLoader)

    lower_limit = args["lower_limit"]
    upper_limit = args["upper_limit"]

    datdim = int(args["datadim"])
    pardim = args["pardim"]
    data_mean = np.load(args["data_mean_dir"])["data_mean"]
    prior_mean = list(args["prior_mean_dict"].values())
    covmat = np.load(args["covmat_dir"])["covariance"]
    g_matrix = np.load(args["g_matrix_dir"])["g_matrix"]
    prior_matrix = np.load(args["prior_matrix_dir"])["prior_matrix"]

    output_file = open('new_%i.init' % datdim, 'w')

    string0 = '// This file lists all the initial parameters needed to execute ' \
            + 'all four\noperation modes in robustness_main.cc. The program ' \
            + 'reads from here before\nexecution. \n'

    string0 += '\n// This semicolon is required for correct parsing.\n'
    string0 += ';\n'

    string1 = 'datdim = ' + str(datdim) + ';\n'
    string1 += generate_covariance_string(datdim, cov=covmat)
    string1 += generate_data_mean_string(datdim, data_mean)
    string1 += '\npardim = ' + str(pardim) + ';\n'

    string2 = generate_prior_mean_string(prior_mean)
    string2 += generate_prior_matrix_string(pardim, prior_matrix)
    string2 += generate_design_matrix_string(pardim, datdim, g_matrix=g_matrix)

    string3 = '\nlower_limit = ' + str(lower_limit) + ';\n'
    string3 += 'upper_limit = ' + str(upper_limit) + ';\n'
    string3 += 'num_steps_analyt_distribution = 1e4;\n'
    string3 += 'num_sampled_data = 1e5\n'
    string3 += '\n// No semicolon here'

    output_file.write(string0 + string1 + string2 + string3)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Init file generator")
    parser.add_argument("--globals", type=str, help="Path to the yaml file")

    args = parser.parse_args()
    print("args.globals", args.globals)

    main(args)
