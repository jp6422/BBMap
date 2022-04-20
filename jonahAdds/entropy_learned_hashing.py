import pandas as pd
import numpy as np
import math
import sys


def read_input(filepath):
    """
    reads in csv for now of FastQC overrepresented sequences format

    input: filepath to csv format

    output: dataframe of (STRING sequence, FLOAT relative_abundance)

    can make this more robust to dif filetypes or actually incorporate FastQC file like we might need in bash script
    """
    #band aid fix...easier to still have in np array format for now tho
    outdata = pd.read_csv(filepath)
    outdata2 = outdata.to_numpy()

    return (outdata2)

def estimate_weighted_entropy(data_array, base_mapper):
    """
    gives an estimate for entropy weighted by relative abundance. The goal will be to use this to determine optimal
    partial key based on overrepresented seqs

    data_array must have been transformed to relative abundance and have characters replaced with integers in seq

    populate matrix with relative probabilities corresponding to each sequence and then calucate entropy by byte
    """

    #initialize array to store most entropic bytes
    #shape is length of longest sequence, number of possible bases G,C,A,T or nothing
    entropy_matrix = np.zeros([len(data_array[0,0]),5])

    
    #pass over each row and add weighted entropy to bucket for each characrter, O(N*max_len)
    for sequence in range(data_array.shape[0]):
        
        #get relative probability, and populate proper array locations
        relative_probability = data_array[sequence,2]
        i=0
        for char in data_array[sequence,0]:

            entropy_matrix[i,base_mapper(char)] += float(relative_probability)
            i+=1

    #calculate renyi entropy by byte in entropy_array, we now have a non-normalized discrete distribution for each byte based on training data
    entropy_array = np.zeros(len(data_array[0,0]))

    #calcualte and append byte entropy to new array, O(max_len)
    for byte in range(len(data_array[0,0])):
        entropy=0
        
        #create variable for row sum to establish probability of each distribution
        byte_sum = sum(entropy_matrix[byte,:])

        for bucket in range(5):

            #calculate sum of squares for entropy sum
            if entropy_matrix[byte,bucket] !=0:
                entropy += (entropy_matrix[byte,bucket]/byte_sum)**2

        #insert negative log of sum of squares into entropy_array
        entropy_array[byte]=-math.log(entropy,2)

    #return sorted entropy array (desc) and array of bytes by entropy (desc) 
    return entropy_array[[np.argsort(-entropy_array)]], np.array(range(len(data_array[0,0])))[[np.argsort(-entropy_array)]]

   

def base_mapper(base):

    if base=='G':
        return 1
    elif base=='C':
        return 2
    elif base=='A':
        return 3
    elif base=='T':
        return 4
    else:
        return 0


def weighted_next_byte (data_matrix):
    """
    next byte alg from hentschel paper implemented with abundance values

    sequences should be sorted by relative abundance in descending order to avoid having to do lookup again
    """

    #initialize comparison variables
    min_coll,min_i = math.inf(),-1

    #create partial key (list of indices)
    partial_key=[]

    #for now assume all same length
    for i in range(len(data_matrix[0,0])):

        #create new partial key to examine
        iter_key = partial_key+[i]
        #initialize variables
        count_table, num_coll = dict(), 0

        for j in range(data_matrix.shape[0]):

            #create table key of the proper indices retrieved from sequence we are examining
            table_key = tuple(data_matrix[0][x] for x in iter_key)

            if table_key in count_table:
                count_table[table_key]+=float(data_matrix[i,2])


def main(csv_filepath, array_curoff):

    data_array = read_input(csv_filepath)

    a,b = estimate_weighted_entropy(data_array)

    return b[:array_cutoff]


if __name__ == "__main__":

    csv_filepath = sys.argv[1]
    array_cutoff = sys.argv[2]

    data_array = read_input(csv_filepath)

    entropies, bytes = estimate_weighted_entropy(data_array, base_mapper)

    main (csv_filepath, array_cutoff)



    

    


