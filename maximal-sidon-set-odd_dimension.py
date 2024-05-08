"""
In this script, we determine large 4-general sets / Sidon sets over AG(n,3), using the (x,x^2) construction in (even) dimension n-1. 
@author: Tim Neumann
"""

import galois
import numpy as np
import subprocess
import os
import argparse

def main():

    if (n % 2 != 0):
        print("Please enter odd integer n.")
        exit()
    
    GF = galois.GF(3**int(n/2), repr="int")
    print(GF.properties)

    sidon_set = []

    # adding points (0,x,x^2) for x in field of order n/2
    for i in range(GF.order):
        element = GF(i)
        element_vec = element.vector()
        square_element = element ** 2
        square_element_vec = square_element.vector()
        
        # Concatenate the vectors and add a zero to embed in respective space
        sidon_element = np.concatenate(([0],np.array(square_element_vec)[::-1],np.array(element_vec)[::-1]))
        sidon_set.append(sidon_element)
        
    # We use the points of the new construction to initialize the search for 4-general sets over AG(n,3)   
    
    print(f"Size of Sidon set for initialization: {len(sidon_set)}")
 
    sidon_set = ''.join(map(str, sidon_set))
    print(sidon_set)
    N = str(n+1)

    command = ["python", "m-general-sets-AG-max_search.py", "-m 4", "-q 3", "-n", N, "-p", sidon_set ]

    os.environ["PYTHONUNBUFFERED"] = "TRUE"

    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True) as process:
        for line in process.stdout:
            print(line, end='') 
        for line in process.stderr:
            print(line, end='')

if __name__ == '__main__':
    
    # Parsing terminal input
    parser = argparse.ArgumentParser(description="Finding m-general sets in AG(n,q) and PG(n,q)")
    parser.add_argument("-n","--n", type=int, required = True, help="dimension over which we want to find maximal Sidon set in AG(n,3) the space over F_q")
   
    args = parser.parse_args()

    global n
    n = args.n - 1

    main()
    
