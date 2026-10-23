import nd2
import numpy as np
from scipy.io import savemat

# import nd2 file list for conversion
file_list = open('nd2_path_list.txt','r')
nd2_list = file_list.readlines()
nd2_list = [x.strip() for x in nd2_list]
# print(nd2_list)
# print(nd2_list[0])

for nd2_file_path in nd2_list:
    my_array = nd2.imread(nd2_file_path)
    mat_file_path = nd2_file_path
    mat_file_path = mat_file_path[0:-3]
    mat_file_path = mat_file_path + 'mat'
    savemat(mat_file_path, {'mov': my_array})