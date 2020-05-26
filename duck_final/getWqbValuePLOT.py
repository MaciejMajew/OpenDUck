import numpy as np
import os
import matplotlib.pyplot as plt

def get_Wqb_value(file_duck_dat):
    f = open(file_duck_dat,'r')
    data = []
    for line in f:
        a = line.split()
        data.append([float(a[1]), float(a[3]), float(a[5]), float(a[8])])
    f.close()
    data = np.array(data[1:])
    Work = data[:,3]
    #split it into segments of 200 points 
    num_segments = int(len(data)/200) 
    num_segments = int(len(data)/200) 
    #alayze each segment to see if minimum in the segment is the local minimum
    #local minimum is the point with the lowest value of 200 neighbouring points
    #first local minumum is miminum used later to duck analysis
    for segment in range(num_segments):
        #detecting minium inthe segment
        sub_data = data[segment * 200 : (segment + 1) * 200]
        sub_Work = sub_data[:,3] 
        index_local = np.argmin(sub_Work)
        #segment of 200 points arround detected minimum
        index_global = index_local + segment * 200
        if index_global > 100:
            sub2_data = data[index_global - 100 : index_global + 101]
        else:
            sub2_data = data[0 : index_global + 101]
        sub2_Work = sub2_data[:,3]
        index_local2 = np.argmin(sub2_Work)
        if index_global < 100:
            if index_local2 == index_global:
                
                Wqb_min_index = index_global
            break
        else:
            if index_local2 == 100:
                Wqb_min_index = index_global
                break
    
    Wqb_min = Work[Wqb_min_index]
    sub_max_data = data[Wqb_min_index:]
    sub_max_Work = sub_max_data[:,3]
    
    
    Wqb_max = max(sub_max_Work)
    
    Wqb_value = Wqb_max - Wqb_min
    return(Wqb_value, data, Wqb_min)

def get_Wqb_value_all():
    file_list = []
    for fil in os.listdir(os.getcwd()):
        if fil[-3:] == 'dat':
            file_list.append(fil)

    Wqb_values = []
    plt.figure(figsize = (7,7))
    for fil in file_list:
        Wqb_data = get_Wqb_value(fil)
        Wqb_values.append(Wqb_data[0])
        plt.plot(10*Wqb_data[1][:,0], Wqb_data[1][:,3]-Wqb_data[2])

    plt.xlabel('HB Distance (A)')
    plt.ylabel('Work (kcal/mol)')
    plt.savefig('wqb_plot.png')
    Wqb = min(Wqb_values)
    return(Wqb)

if __name__ == '__main__':
    print(get_Wqb_value_all())