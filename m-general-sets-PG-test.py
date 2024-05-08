"""
In this script, we test if a set is m-general over projective space PG(n,q).
@author: Tim Neumann
"""

from helperFunctions import *
import argparse

def create_PG_from_AG(n,q):
    '''
    Function that creates points in projective space by adding subspaces at infinity for the affine space of equal dimension.
    (Note that we could also create these points by applying the eqiuvalence relation on the vector space of higher dimension.)
    We write the q^n+q^(n-1)+...+q+1 elements in PG that are created with this function using homogeneous coordinates.
        
    Input: n (dimension of the space) and q (order of the underlying field for affine space)
    Output: projective points in homogeneous coordinates
    '''
    
    pointsPG = []
    for i in range(n+1):
        for j in range(q**i):
            qary = [None] * (n+1)
            quotient = j
            for k in range(1,i+1):
                qary[i-k] = int(quotient % q)
                quotient = int(quotient // q)
            qary[i] = 1
            for l in range(i+1,n+1):
                qary[l] = 0
            pointsPG.append(qary[::-1])       
    return pointsPG

def point_PG_to_AG(pointPG,q): 
    '''
    This function transforms a point in projective space PG(n,q), expressed in its homogeneous coordinates,
    to a point in affine space AG(n+1,q). This is done by rescaling the homogeneous coordinates based on a non-zero anchor.
    Note that by this way of expressing points in projective space, the transformation viewed as a function
    is injective but not surjective, i.e., no two points in projective space will give the same vector in affine space, but not every
    point in affine space will be covered by this construction. Most importantly, it is well-defined.
    
    Input: point in projective space, order of the underlying field
    Output: point in affine space
    '''
        
    pointAG = [None]*len(pointPG)
    homogeneous_anchor = 0 # this anchor will be changed since we exclude the zero vector in the input domain
    for i in range(len(pointPG)-1,-1,-1):
        pointAG[i] = pointPG[i]
        if pointAG[i] != 0:
            homogeneous_anchor = pointAG[i]
            loc_anchor = i
                
    for i in range(loc_anchor,len(pointPG)):
        pointAG[i] = (pointPG[i] * inverse_finite_field(homogeneous_anchor,q) ) % q

    return pointAG

def f_is_0(subsetOfPointSet):
    '''
    Function that determines if the there exists a (nonzero) linear combination of the points in the set that sums to zero, which implies
    that the set is linearly dependent.
    
    Input: set of points in affine space of dimension n+1
    Output: truth value 
    '''

    for c in coefficientVec:    # iterate through all possible linear combinations
        if linear_comb(c,subsetOfPointSet,q) == [0]*(n+1):
            return True

    return False
    
def linearly_independent(setOfPoints, m):
    '''
    Function that determines if set of points in vector space F_q^{n+1} is linearly independent.
    
    The nested for-loops traverse through all subsets of size m of the given set of points. Now testing if those m points
    are linearly independent verifies the m-general property in projective space.
    
    Input: set of points in vector space F_q^{n+1}
    Output: truth value for linear independence
    '''
    
    if (m == 3):

        for i in range(0,len(setOfPoints)):
            for j in range(i+1,len(setOfPoints)):
                for k in range(j+1, len(setOfPoints)):
                    subset = [setOfPoints[i],setOfPoints[j],setOfPoints[k]]      
                    if f_is_0(subset) == True:
                        return False
        return True

    elif (m == 4):

        for i in range(0,len(setOfPoints)):
            for j in range(i+1,len(setOfPoints)):
                for k in range(j+1, len(setOfPoints)):
                    for l in range(k+1, len(setOfPoints)):
                        subset = [setOfPoints[i],setOfPoints[j],setOfPoints[k],setOfPoints[l]] 
                        if f_is_0(subset) == True:
                            return False
        return True

    elif (m == 5):
            
        for i in range(0,len(setOfPoints)):
            for j in range(i+1,len(setOfPoints)):
                for k in range(j+1, len(setOfPoints)):
                    for l in range(k+1, len(setOfPoints)):
                        for m in range(l+1, len(setOfPoints)):
                            subset = [setOfPoints[i],setOfPoints[j],setOfPoints[k],setOfPoints[l],setOfPoints[m]] 
                            if f_is_0(subset) == True:
                                return False
        return True
    
    else:
        return "Please set m to an integer in {3,4,5}"
      
def checkPoints(set, m):
    '''
    Function that determines if set of points in projective space forms an m-general set.
    
    Input: set of points in projective space, parameter m for m-general set
    Output: possible choices up to the knowledge of all points in the set
    '''
        
    new_set = []
    for point in set:
        new_set.append(point_PG_to_AG(point,q))    
    return linearly_independent(new_set,m)
       
def main():
    '''
    Main function that initializes the parameters and starts the check for m-general property in PG(n,q) 
    '''
    
    # Determining coefficient vector to test for linear independence
    global coefficientVec
    coefficientVec = createCoefficientVec(m,q)
    
    result = checkPoints(points, m)
    print(f"The set of points {points} is a {m}-general set in PG({n},{q})? {result}. ")
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Finding m-general sets in AG(n,q) and PG(n,q)")

    parser.add_argument("-m", "--m", type=int, help="m for m-general set")
    parser.add_argument("-q","--q", type=int, help="q for order of the underlying field")
    parser.add_argument("-n","--n", type=int, help="dimension of the space over F_q")
    # Points need to be added in homogeneous coordinates
    # SAMPLE INPUT: -m 4 -q 3 -n 3 -p 0 0 0 1 -p 0 0 1 0 -p 0 1 0 0 -p 1 0 0 0 -p 1 1 1 1 
    parser.add_argument("-p","--points", type=int, nargs="+", action="append", help="points in projective space PG(n,q) to check for m-general property")

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
    
    print(f"Points to be checked: {points}")
    
    main()
    
    

