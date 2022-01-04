import csv
import numpy as np
import math
import matplotlib.pyplot as plt
import statistics 
import scipy.stats
import scipy

# 
# EMBRY-RIDDLE AERONAUTICAL UNIVERSITY
# MA412 - PROBABILITY AND STATISTICS
# FINAL PROJECT
# Jose Nicolas Gachancipa
# Spring 2018
# 
print('\n#####################\n')
print('Embry-Riddle Aeronautical University')
print('Department of Mathematics')
print('MA412 - Probability and Statistics')      
print('Code developed by: Jose Nicolas Gachancipa')
print('Purpose: Establish the elevation threshold for a GPS receiver using CMC.')
print('\n#####################\n')

# Inputs.  
min_value = 0
max_value = 90
groups_range = range(90,1000,100)        # [Start, End, Increment]
significance_level = 0.001
reverse = 0
levene = 0
brown_forsynthe = 1
    
# Open the Excel file and extract the info.
x_axis_column = 0
y_axis_column = 1
header_cutoff_row = 0
with open('CMC_DATA.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    count = 1
    x_axis = []
    y_axis = []
    for row in csv_reader:
        if count>header_cutoff_row:
            x_axis.append(float(row[x_axis_column]))
            y_axis.append(float(row[y_axis_column]))
        count = count + 1

# Sort the values using the x_axis as a baseline.
original_y_axis = [t for _,t in sorted(zip(x_axis,y_axis))]
original_x_axis = sorted(x_axis)
plt.plot(x_axis,y_axis,'o')

#  Subdivide the array. 
def divide_1(x,y,number_of_groups):
    divisions = np.linspace(min_value,max_value,number_of_groups+1)
    output_x = []
    output_y = []
    for i in range(len(divisions)-1):
        lower_limit = divisions[i]
        upper_limit = divisions[i+1]
        x_vector = []
        y_vector = []
        count = 0
        for value in x:
            if value>=lower_limit and value<upper_limit:
                x_vector.append(value)
                y_vector.append(y[count])
            count = count + 1
        if len(x_vector) > 0:
            output_x.append(x_vector)
            output_y.append(y_vector)
    return [output_x,output_y]

# Define functions.
# Function:
def divide_2(vector, number_of_divisions):
    a = np.linspace(0,len(vector),number_of_divisions+1)
    out = []
    for i in range(len(a)-1):    
        out.append(vector[math.floor(a[i]):math.floor(a[i+1])])
    return out

# Run iteratevely for multiple divisions.
all_cut_off = []
for groups in groups_range:

    # Divide the data into the number of groups specified by the user.
    [input_x_axis, input_y_axis] = divide_1(original_x_axis,original_y_axis,groups)
    #input_x_axis = divide_2(original_x_axis, groups)
    #input_y_axis = divide_2(original_y_axis, groups)
    if reverse == 1:
        input_x_axis = input_x_axis[::-1]
        input_y_axis = input_y_axis[::-1]
    
    # Run the test progressively, until the variances are not equal.
    for i in range(2,len(input_x_axis)+1):
        
        # Set the arrays correspondingly first.
        divided_x_axis = input_x_axis[:i]
        divided_y_axis = input_y_axis[:i]
        x_axis = [j for i in divided_x_axis for j in i]
        y_axis = [j for i in divided_y_axis for j in i]
    
        # Compute the means and lengths of the lists as needed.
        mean = sum(y_axis)/len(y_axis)
        Nis = [len(x) for x in divided_x_axis]
        medians = [statistics.median(j) for j in divided_y_axis]
        means = [sum(j)/len(j) for j in divided_y_axis]
        N = len(x_axis)   # Total number of cases in all groups.
        k = len(divided_x_axis)        # Number of groups.
        
        # Compute the degrees of freedom.
        v1 = N-k
        v2 = k-1
        
        # Find (Zi.) and (Zi..).
        Zi_dots = []
        zi_doubledot = 0
        count = 0
        for group in divided_y_axis:   
            if levene == 1:
                local_parameter = means[count] # Use the mean for levene's test.
            elif brown_forsynthe == 1:
                local_parameter = medians[count] # Use the median for brown-forsynthe.
            element_differences = []
            for element in group:
                element_differences.append(abs(element - local_parameter))
            zi_doubledot = zi_doubledot + sum(element_differences)
            average_differences = sum(element_differences)/len(element_differences)
            Zi_dots.append(average_differences)    
            count = count + 1
        zi_doubledot = zi_doubledot/len(x_axis)
        
        # Compute the numerator of the function using (Zi.) and (Zi..).
        count = 0
        numerator = 0
        for zi_dot in Zi_dots:
            local_length = Nis[count]
            numerator = numerator + (local_length*((zi_dot-zi_doubledot)**2))
            count = count + 1
        numerator = numerator*(v1)
        
        # Compute the denominator.
        count = 0
        denominator = 0
        for group in divided_y_axis: 
            if levene == 1:
                local_parameter = means[count] # Use the mean for levene's test.
            elif brown_forsynthe == 1:
                local_parameter = medians[count] # Use the median for brown-forsynthe.
            element_differences = 0
            local_zdot = Zi_dots[count]
            for element in group:
                element_differences = element_differences + ((abs(element - \
                                            local_parameter)-local_zdot)**2)
            denominator = denominator + element_differences
            count = count + 1
        denominator = denominator*(v2)
        
        # Determine the f critical value.
        F = scipy.stats.f.ppf(q=1-significance_level, dfn=v1, dfd=v2)
        
        # Compute W.
        W = numerator/denominator  
        
        # Run until we reject the null hypothesis.
        # The null hypothesis establishes that the variances of all selected groups are equal.
        # The null hypothesis is rejected when W is bigger than F.
        # W = Test value
        # F = Critical falue of the F distrbution
        if W>F:
            cutoff_threshold =  divided_x_axis[-1][-1]
            print('k:',groups,'| Cut-off elevation threshold: ',cutoff_threshold,'degrees.')
            all_cut_off.append(cutoff_threshold)
            break
    
# Plot.
r = [j for i in divided_x_axis[:-1] for j in i] 
q = [j for i in divided_y_axis[:-1] for j in i]
plt.plot(r,q,'ro') 
plt.xlabel('Elevation (degrees)')
plt.ylabel('Code-Minus-Carrier (CMC)')
plt.show()    
plt.plot(groups_range,all_cut_off)    
plt.ylim((0,90))
plt.xlabel('# of Windows')
plt.ylabel('Elevation Threshold (degrees)')
plt.show()   
        
    
    
