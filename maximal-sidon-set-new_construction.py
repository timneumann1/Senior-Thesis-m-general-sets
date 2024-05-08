"""
In this script, we determine large 4-general sets / Sidon sets over AG(n,3), using the new construction
developed in the senior thesis. 
@author: Tim Neumann
"""

import galois
import numpy as np
import subprocess
import os
import argparse

def main():
    points = sidon_recursion(n)
    for p in points:
        sidon_set.append(p)
    
# Defining the recursion    
def sidon_recursion(n):
    # Base Case 1: n is even
    if (n % 2 == 0):
        return sidon(n)
    # Base Case 2: n = 1
    if (n == 1):
        return [[0],[1]]
    # Entering the recursion
    else:
        n = n-1
        list = []
        new_points_rec = sidon(n)
        for p in new_points_rec:
            list.append( np.concatenate(([0],p)) )
        
        new_points_rec = sidon_recursion(int(n/2))
        for p in new_points_rec:
            list.append( np.concatenate(([1],p,np.zeros(n//2, dtype = int))) )
        return list
    
def sidon(n):
    
    # In even dimension n, we construct a maximal Sidon set of size 3^(n/2)
    
    GF = galois.GF(3**int(n/2), repr="int")
    sidon_points = []

    for i in range(GF.order):
        
        element = GF(i)
        element_vec = element.vector()
        square_element = element ** 2
        square_element_vec = square_element.vector()
        
        sidon_element = np.concatenate((np.array(square_element_vec)[::-1],np.array(element_vec)[::-1]))
        sidon_points.append(sidon_element)

    return sidon_points

if __name__ == '__main__':
    
    # Parsing terminal input
    parser = argparse.ArgumentParser(description="Finding m-general sets in AG(n,q) and PG(n,q)")
    parser.add_argument("-n","--n", type=int, required = True, help="dimension over which we want to find maximal Sidon set in AG(n,3) the space over F_q")
   
    args = parser.parse_args()

    global n
    n = args.n
    N = n
    
    global sidon_set
    sidon_set = []
    
    main()
    
    # We use the points of the new construction to initialize the search for 4-general sets over AG(n,3)
    sidon_set = ''.join(map(str, sidon_set))
    N = str(N)
    command = ["python", "m-general-sets-AG-max_search.py", "-m 4", "-q 3", "-n", N , "-p", sidon_set]

    os.environ["PYTHONUNBUFFERED"] = "TRUE"

    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True) as process:
        for line in process.stdout:
            print(line, end='')  # Print the output from stdout in real-time
        for line in process.stderr:
            print(line, end='')
