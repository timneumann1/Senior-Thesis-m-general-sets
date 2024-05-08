"""
In this script, we determine complete m-general sets over projective space PG(n,q).
@author: Tim Neumann
"""

from helperFunctions import *
import argparse

def create_PG(n,q):
    '''
    Function that creates the q^n+q^(n-1)+...+q+1 points in projective space by adding all possible
    homogeneous coordinates of projective space.
        
    Input: n (dimension of the space) and q (order of the underlying field for affine space)
    Output: projective points in homogeneous coordinates
    '''
    
    pointsPG = []
    for i in range(n+1):
        for j in range(q**i):
            qary = [None] * (n+1)
            quotient = j
            for k in range (0,i):
                qary[k] = int(quotient % q)
                quotient = int(quotient // q)
            qary[i] = 1
            for l in range(i+1,n+1):
                qary[l] = 0
            pointsPG.append(qary[::-1])   
    return pointsPG

# DISCLAIMER: This function is not used but could be used for other purposes.
def point_VS_to_PG(pointVS,q): 
    '''
    This function transforms a point in the vector space F_q^(n+1) to a point in projective space PG(n,q),
    expressed in its homogeneous coordinates, by projecting according to the usual equivalence relation.
    Note that by this way of expressing points in projective space, the transformation viewed as a function
    is surjective but not injective, since every point in projective space will be mapped onto by some point in affine space, but
    two linearly dependent vectors can be mapped to the same point. Most importantly, the map is well-defined only
    if we exclude the zero vector from the domain (that is, the parameter pointVS cannot assume the zero vector as value).
        
    Input: point in vector space
    Output: point in projective space
    '''
    
    if (pointVS == [0]*len(pointVS)):
        return "The zero vector cannot be used as input to this function"
    
    pointPG = [None]*len(pointVS)
    homogeneous_anchor = 0 # this anchor will be changed since we exclude the zero vector in the input domain
    loc_anchor = len(pointVS)-1
    for i in range(len(pointVS)-1,-1,-1):
        pointPG[i] = pointVS[i]
        if pointVS[i] != 0:
            homogeneous_anchor = pointVS[i]
            loc_anchor = i
    for i in range(loc_anchor,len(pointPG)):
        pointPG[i] = (pointPG[i] * inverse_finite_field(homogeneous_anchor,q) ) % q

    return pointPG

def point_PG_to_VS(pointPG,q): 
    '''
    This function transforms a point in projective space PG(n,q), expressed in its homogeneous coordinates,
    to a point in the vector space F_q^(n+1). This is achieved by rescaling the homogeneous coordinates based on a non-zero anchor.
    Note that by this way of expressing points in projective space, the transformation viewed as a function
    is injective but not surjective, i.e., no two points in projective space will yield the same vector in the higher-dimensional
    vector space, but not every point in the vector space will be covered by this construction. 
    
    Input: point in projective space, order of the underlying field
    Output: point in vector space
    '''
        
    pointVS = [None]*len(pointPG)
    homogeneous_anchor = 0 # this anchor will be changed during execution since we exclude the zero vector in the input domain
    for i in range(len(pointPG)-1,-1,-1):
        pointVS[i] = pointPG[i]
        if pointVS[i] != 0:
            homogeneous_anchor = pointVS[i]
            loc_anchor = i
                
    for i in range(loc_anchor,len(pointPG)):
        pointVS[i] = (pointPG[i] * inverse_finite_field(homogeneous_anchor,q) ) % q

    return pointVS

def f_is_0(subsetOfPointSet, reduc):
    '''
    Function that determines if the there exists a (nonzero) linear combination of the points in the set that sums to zero, which implies
    that the set is linearly dependent.
    
    Input: set of points in vector space of dimension n+1, 
           variable that indicates if we are testing for linear independence of smaller subsets 
    Output: truth value 
    '''

    if reduc == 0:
        for c in coefficientVec:     # Iterate through all possible linear combinations
            if linear_comb(c,subsetOfPointSet,q) == [0]*(n+1):
                return True
            
    elif reduc == 1:
        for c in coefficientVec_m_minusOne:     # Iterate through all possible linear combinations
            if linear_comb(c,subsetOfPointSet,q) == [0]*(n+1):
                return True
            
    elif reduc == 2:
        for c in coefficientVec_m_minusTwo:     # Iterate through all possible linear combinations
            if linear_comb(c,subsetOfPointSet,q) == [0]*(n+1):
                return True

    return False
    
def linearly_independent(setOfPoints, m, reduc = 0):
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
                    if f_is_0(subset, reduc) == True:
                        return False
        return True

    elif (m == 4):

        for i in range(0,len(setOfPoints)):
            for j in range(i+1,len(setOfPoints)):
                for k in range(j+1, len(setOfPoints)):
                    for l in range(k+1, len(setOfPoints)):
                        subset = [setOfPoints[i],setOfPoints[j],setOfPoints[k],setOfPoints[l]] 
                        if f_is_0(subset, reduc) == True:
                            return False
        return True

    elif (m == 5):
            
        for i in range(0,len(setOfPoints)):
            for j in range(i+1,len(setOfPoints)):
                for k in range(j+1, len(setOfPoints)):
                    for l in range(k+1, len(setOfPoints)):
                        for m in range(l+1, len(setOfPoints)):
                            subset = [setOfPoints[i],setOfPoints[j],setOfPoints[k],setOfPoints[l],setOfPoints[m]] 
                            if f_is_0(subset, reduc) == True:
                                return False
        return True
    
    else:
        return "Please set m to an integer in {3,4,5}"
      
def possibleChoices_fun(genSet, possibleChoices, m, first_iter):
    '''
    Function that determines all projective points that can be added to a given m-general set while maintaining its 
    m-general property.
    
    Input: set of m-general points,
           list of all possible choices that could potentially be added to the set 
           parameter m for m-general set,
           boolean variable indicating whether or not discarding lexicographically earlier points
    Output: possible choices up to the knowledge of all points in the set
    '''
        
    record = [] # vector that will contain all points that are marked for removal
    
    for possible_point in possibleChoices:
        new_genSet = [] 
        for i in range(0,len(genSet)):
            new_genSet.append(point_PG_to_VS(genSet[i],q))
        new_genSet.append(point_PG_to_VS(possible_point,q))
        
        # we need to verify that the newly constructed set does not violate m-general, or if possible, (m-1)- or (m-2)-general property
        if (linearly_independent(new_genSet,m) == False) or (m-1 > 2 and linearly_independent(new_genSet,m-1,reduc = 1) == False) or (m-2 > 2 and linearly_independent(new_genSet,m-2,reduc = 2) == False) :
            record.append(possible_point)
    
    # remove all points that are lexicographically earlier than the point added last
    
    if not(opt == "min") and first_iter == False:
        last = qary_vector_to_decimal(genSet[-1],q)
        for index in range(1,last): # the zero vector cannot be transformed to projective coordinates
            record.append(point_VS_to_PG(decimal_to_qary_vector(index,n+1,q),q))
    else:
        pass
  
    # removing all points that are flagged for removal
    
    for entry in possibleChoices[:]:
        if entry in record:
            possibleChoices.remove(entry)
    
    return possibleChoices
    
def buildCap(genSet, possibleChoices, m, first_iter = False):
    '''
    Function that recursively tests all elements that could potentially be added to a set of points while 
    maintaining its m-general property, appending the point if possible. 
    
    Input: set of points that form a m-general set
           possible choices that can potentially be added to this set while maintaining m-general property,
           parameter m for m-general sets,
           boolean variable indicating whether or not discarding lexicographically earlier points

    Output: void
    '''
    
    global maxSize_mgeneralSet
    global max_mgeneralSet
    global minSize_mgeneralSet
    global min_mgeneralSet

    possible = possibleChoices_fun(genSet,possibleChoices,m,first_iter)

    # Base Case: There are no possible points to add without violating the m-general property
    if len(possible) == 0:
        
        if len(genSet) > maxSize_mgeneralSet and opt == "max":
            maxSize_mgeneralSet = len(genSet)
            max_mgeneralSet = genSet.copy()
            print(f"A complete {m}-general set of (max) size {maxSize_mgeneralSet} is given by {max_mgeneralSet}.")
        
        if len(genSet) == maxSize_mgeneralSet:
            print(f"A complete {m}-general set of size {len(genSet)} is given by {genSet}.")
            
        if len(genSet) < minSize_mgeneralSet and opt == "min":
            minSize_mgeneralSet = len(genSet)
            min_mgeneralSet = genSet.copy()
            print(f"A complete {m}-general set of (min) size {minSize_mgeneralSet} is given by {min_mgeneralSet}.")
        
        if len(genSet) == minSize_mgeneralSet:
            print(f"A complete {m}-general set of size {len(genSet)} is given by {genSet}.")
        
        return 
            
    else:
        for i in range (0,len(possible)): 
            
            genSet.append(possible[i])
            newPossible = possible.copy()
            newPossible.remove(possible[i])
            buildCap(genSet, newPossible, m)
            genSet.pop()
        
def main():
    '''
    Main function that initializes the parameters and starts the search for m-general sets in PG(n,q) 
    '''
    
    if (n < m-2):
        print("Please choose n such that n >= m-2 is satisfied.")
        
    global maxSize_mgeneralSet # variable that captures the size of the largest complete m-general set found
    maxSize_mgeneralSet = 0
    global max_mgeneralSet # variable that captures the elements of the largest complete m-general set found
    max_mgeneralSet = []

    global minSize_mgeneralSet # variable that captures the size of the largest complete m-general set found
    minSize_mgeneralSet = float("inf")
    global min_mgeneralSet # variable that captures the elements of the largest complete m-general set found
    min_mgeneralSet = []
    
    # Determining coefficient vector to test for linear independence
    global coefficientVec
    coefficientVec = createCoefficientVec(m,q)
    
    global coefficientVec_m_minusOne
    if m == 4 or m == 5:
        coefficientVec_m_minusOne = createCoefficientVec(m-1,q)
    
    global coefficientVec_m_minusTwo
    if m == 5:
        coefficientVec_m_minusTwo = createCoefficientVec(m-2,q)
    
    # Initializing points in projective space by embedding the affine space of equal dimension 
    points = create_PG(n,q)    

    if mode == "zero":
    # Initializing the search with an empty set 
        start = []
            
    elif mode == "init":
    # Adding as many points as possible to the start list for which we are guaranteed to have a m-general set
        start = [points[0]]
        i = 0
        while ( int((q**(i+1)-1)/(q-1)) < len(points)):
            start.append(points[int((q**(i+1)-1)/(q-1))])
            i += 1        

    possible = points.copy()
    for entry in start:
        possible.remove(entry)
    
    print(f"\n Starting the search for {m}-general sets over PG({n},{q})...\n")

    buildCap(start,possible,m, first_iter = True)
    
    if opt == "max":
        print(f"\n A complete {m}-general set of maximal size {maxSize_mgeneralSet} in PG({n},{q}) is {max_mgeneralSet}.")
    elif opt == "min":
        print(f"\n A complete {m}-general set of minimal size {minSize_mgeneralSet} in PG({n},{q}) is {min_mgeneralSet}.")

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Finding m-general sets in PG(n,q)")

    parser.add_argument("-m", "--m", type=int, required = True, help="m for m-general set")
    parser.add_argument("-q","--q", type=int, required = True, help="q for order of the underlying field")
    parser.add_argument("-n","--n", type=int, required = True, help="dimension of the space over F_q")
    parser.add_argument("-mode","--mode", type=str, required = True, help="zero or init: mode of initialization")
    parser.add_argument("-opt","--opt", type=str, required = True, help="min or max: mode of optimization")

    args = parser.parse_args()

    global m # positive integer that defines the m in m-general sets
    m = args.m
    global q # positive prime integer that denotes the size of the underlying field
    q = args.q
    global n # positive integer greater than 1 (we are not interested in projective lines) that encodes the dimension of the space over F_q that we are finding m-general sets on 
    n = args.n
    global mode 
    mode = args.mode
    global opt
    opt = args.opt
    
    main()
    
    

