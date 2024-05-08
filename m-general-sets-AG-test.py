"""
In this script, we test whether or not a given set is m-general over affine space AG(n,q).
@author: Tim Neumann
"""

from helperFunctions import *
import argparse

def f_is_0(subsetOfPointSet):
    '''
    Function that determines if the function f defined in Theorem 2.5 of "Improved bounds on sizes of generalized caps in AG(n,q)"
    by Tait and Won evaluates to zero when given the input of some m+1 points, where combinations of those points are governed
    by C_0^(m+1) as defined in the aforementioned paper.
    
    Input: set of points in affine space of dimension n+1
    Output: truth value 
    '''

    # iterate through all possible linear combinations
    for c in affine_comb_Coeff:
        if affine_comb(c,subsetOfPointSet,q) == [0]*(n):
            return True
    return False

def weakly_avoids_0(set, m):
    '''
    Function that determines whether or not a given set weakly avoids zero. 
    The nested for-loops iterate over every m-tuple in the given set and checks whether or not the points in the set
    affinely combine to zero.
    
    Input: set of points in affine space
    Output: truth value for weak avoidance
    '''
    
    global affine_comb_Coeff
    
    if (m == 3):
        
        for i in range(0,len(set)-2):                
            for j in range(i+1,len(set)-1):  
                for k in range(j+1, len(set)):
                    subset = [set[i],set[j],set[k]]
                    if f_is_0(subset) == True:
                        return False
        return True        
    
    elif (m == 4):
        
        for i in range(0,len(set)-3):
            for j in range(i+1,len(set)-2):
                for k in range(j+1, len(set)-1):
                    for l in range(k+1, len(set)):
                        subset = [set[i],set[j],set[k],set[l]]
                        if f_is_0(subset) == True:
                            return False
        return True
        
    elif (m == 5):
        
        for i in range(0,len(set)-4):
            for j in range(i+1,len(set)-3):
                for k in range(j+1, len(set)-2):
                    for l in range(k+1, len(set)-1):
                        for m in range(l+1, len(set)):
                            subset = [set[i],set[j],set[k],set[l],set[m]]
                            if f_is_0(subset) == True:
                                return False
        return True
        
    else:
        return "Please set m to an integer in {3,4,5}"
      
    
def checkPoints(set, m):
    '''
    Function that determines if set of points in affine space forms an m-general set.
    
    Input: set of m-general points, parameter m for m-general set
    Output: possible choices up to the knowledge of all points in the set
    '''
    
    weakly_avoids = weakly_avoids_0(set,m)
     
    return weakly_avoids
        
def main():
    '''
    Main function that initializes the parameters and starts the search for m-general sets in AG(n,q) 
    '''
    
    # The original general set (consisting of affinely independent vectors in F_q^n) must have
    # cardinality at least m-1, since the inequality q^n >= m-1 must be satisfied in order for the 
    # program to work (otherwise the for loops will not be triggered). Note however that this is
    # not a very restrictive constraint for small m, since if q^n < m-1 for m = 3,4,5 and some prime q,
    # then n is 1 or 2, and these cases are trivial
    
    if (q**n < m-1):
        print("Please choose n such that q^n >= m-1 is satisfied.")
    
    global affine_comb_Coeff
    affine_comb_Coeff = coeff_vector_C0(m,q)
    
    result = checkPoints(pointsAG, m)
    print(f"The set of points {pointsAG} is a {m}-general set in AG({n},{q})? {result}. ")
    
if __name__ == '__main__':
    
    # Parse terminal input
    parser = argparse.ArgumentParser(description="Finding m-general sets in AG(n,q) and PG(n,q)")

    parser.add_argument("-m", "--m", type=int, help="m for m-general set")
    parser.add_argument("-q","--q", type=int, help="q for order of the underlying field")
    parser.add_argument("-n","--n", type=int, help="dimension of the space over F_q")
    # Points need to be added as decimals
    # SAMPLE INPUT: -p 0 28 55 90 124 148 171 
    parser.add_argument("-p","--points", type=int, nargs="+", action="append", help="points in affine space AG(n,q) to check for m-general property")
    
    # Parsing terminal input
    args = parser.parse_args()

    global m # positive integer that defines the m in m-general sets
    m = args.m
    global q # positive prime integer that denotes the size of the underlying field
    q = args.q
    global n # positive integer greater than 1 (we are not interested in projective lines) that encodes the dimension of the space over F_q that we are finding m-general sets on 
    n = args.n
    global points
    points = args.points
    global pointsAG
    pointsAG = []
    
    for point in points[0]:
        pointsAG.append(decimal_to_qary_vector(point,n,q))
    
    print(f"Points to be checked: {pointsAG}")
    
    main()

