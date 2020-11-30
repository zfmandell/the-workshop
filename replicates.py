# runs on python3
# requires numpy, pandas, matplotlib

# will take the result of a differential expression analysis, and a file containing TPM expression values of replicates
# and output a graph of the mean of the log10 TPM values as compared to their log2FC values
# use CSV files as input files
# first input is TPM file
# second input is DE file
# third input is name of graph to be saved

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import sys

#defining a function that we can use later on in the main method, fyle is the input for the function
def read_TPM(fyle):
    return_dict = {}
    """
    data structures used : list, dictionary

    list is exactly what it sounds like, a list of things, the order is always retained
    list = [a,b,c,d]
            0,1,2,3

    list at 0 = a (list[0] = a)
    list at 1 = b (list[1] = b)
    list at 1 through 2 = b,c (list[1:3] = b,c)

    dictionary has the organization Key:Value

    dict = {a:1,b:2,c:3}

    dict at a = 1 (dict[a] = 1)
    dict at c = 2 (dict[b] = 2)

    asking for a key will give you its value
    there can only be one of each key
    very quick for finding things
    """
    with open(fyle,'r') as inp:
        #as you iterate through a file, the location will be remembered and only cycled through once
        firstline = inp.readline()
        #will iterate through all remaining lines in the fyle
        for line in inp:
            """
            start) 'a,2,3,4\n'
            1) 'a,2,3,4'
            2) ['a','2','3','4']
            transcript = 'a'
            TPM_values = ['2','3','4']

            #strip(input) will remove input, no input will remove newline special character
            #split(input) will divide line by input, returning a list of values between input
            """
            contents = line.strip().split(',')
            transcript = contents[0]
            TPM_values = contents[1:]
            """
            1) take all TPM values, convert each from string to float (numeric data type)
            from :  ['2','3','4']
            to : [2,3,4]
            2) log10 transform each value
            3) calculate the mean of the 3 values
            4) add the mean to a dict, with key : value (transcript:mean)

            list comprehension - [function(x) for x in list]
            make list
            for x in list:
                do thing
                add to list
            return new list

            mapping - map(function,list)
            make list
            for x in list
                do thing
                add to list
            return new list
            """
            expression = np.mean([np.log10(x) for x in map(float,TPM_values)])
            #make new key:value pair in return dict
            return_dict[transcript] = expression

    #once the function has ran, return dict will be returned
    return return_dict

def read_DE(fyle):
    return_dict = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            """
            try except switches
            will try and execute the code within the try block
            if successful, great, move on to next iteration of for loop
            if not succesful, check the exception that was raised
            if the exception matches the exception specified in except
            in this case 'ValueError'
            will execute the code in the except block
            """
            try:
                transcript = line.strip().split(',')[0]
                DE = abs(float(line.strip().split(',')[2]))
            except ValueError:
                pass
            return_dict[transcript] = DE
    return return_dict

def main():
    #sys.argv[] will allow you to specify arguments at the CLI
    TPM = sys.argv[1]
    DE = sys.argv[2]
    output = sys.argv[3]

    TPM_dict = read_TPM(TPM)
    DE_dict = read_DE(DE)

    """
    need final dictionary to look like {['TPM values'] = [list of TPM values]
                                  ['DE values'] = [list of DE values]
                                 }
    each list need to reference the same transcript at each index
    ie)
    [list of TPM values] at index 0 references transcript A
    [list of DE values] at index 0 references transcript A

    this will allow easy conversion from a dictionary to a pandas dataframe that can be loaded into matplotlib functions
    so we first need to make the two lists
    and then add them into a final dict as key:value pairs with the above keys
    """

    #initialize final dict and two lists
    final_dict,sub_TPM,sub_DE = {},[],[]
    #iterate through key:value pairs of the DE dict
    for key,value in DE_dict.items():
        #append the value to the new DE list
        sub_DE.append(value)
        #access the value at the transcript (key) of the TPM dict, add it to the new TPM list
        sub_TPM.append(TPM_dict[key])
    #make new key:value pairs where key = column name, value = column contents
    final_dict['log10 mean TPM in WT'] = sub_TPM
    final_dict['log2FC'] = sub_DE

    # will convert the dict into a pandas dataframe
    df = pd.DataFrame(final_dict)

    #will plot the dataframe as a scatter plot, x axis and y axis label will be column names, dots will be blue
    df.plot.scatter(x='log10 mean TPM in WT',
                   y='log2FC',
                   c='blue')

    #will save the plot to the working directory as name output
    plt.savefig(output)

# calls main method, will only run this file if running program directory, not as an additional module
if __name__ == '__main__':
    main()
