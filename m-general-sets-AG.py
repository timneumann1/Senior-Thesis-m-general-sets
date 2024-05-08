"""
In this script, we determine complete m-general sets over affine space AG(n,q).
@author: Tim Neumann
"""

from helperFunctions import *
import argparse

def possibleChoices_fun(genSet, possibleChoices, newChoice, m, lex = True):
    '''
    Function that determines all points that can be added to a given m-general set while maintaining its 
    m-general property.
    The nested for-loops iterate over every m-2 tuple in the given set and append it with the last vector that has been 
    added to the set. For the resulting set of m-1 vectors, we then determine every vector among the possible vectors that 
    can be obtained as an affine combination of the vectors in the set, which are marked for removal.
       
    
    Input: set of m-general points, 
           list of all possible choices that could potentially be added to the set (up to the knowledge of adding the last point),
           last entry that was added to the set,
           parameter m for m-general sets, 
           boolean variable indicating whether or not discarding lexicographically earlier points
    Output: possible choices up to the knowledge of all points in the set
    '''
    
    global affine_comb_Coeff
   
    record = [] # vector that will contain all points that are marked for removal
    
    new_genSet = []
    for i in range(0,len(genSet)):
        new_genSet.append(genSet[i])
    
    if (m == 3):
        
        for i in range(0,len(new_genSet)):                
                        
            points = [new_genSet[i],newChoice]
            
            for l in range(len(affine_comb_Coeff)):
                aff_comb = affine_comb(affine_comb_Coeff[l],points,q)
                record.append(qary_vector_to_decimal(aff_comb,q))
    
    elif (m == 4):
        
        for i in range(0,len(new_genSet)-1):
            for j in range(i+1,len(new_genSet)):
                
                points = [new_genSet[i],new_genSet[j],newChoice]
                
                for l in range(len(affine_comb_Coeff)):
                    aff_comb = affine_comb(affine_comb_Coeff[l],points,q)
                    record.append(qary_vector_to_decimal(aff_comb,q))
        
    elif (m == 5):
        
        for i in range(0,len(new_genSet)-2):
            for j in range(i+1,len(new_genSet)-1):
                for k in range(j+1, len(new_genSet)):
                
                    points = [new_genSet[i],new_genSet[j],new_genSet[k],newChoice]

                    for l in range(len(affine_comb_Coeff)):
                        aff_comb = affine_comb(affine_comb_Coeff[l],points,q)
                        record.append(qary_vector_to_decimal(aff_comb,q))
        
    else:
        return "Please set m to an integer in {3,4,5}"
      
    if lex == True and not (opt == "min"): # if flag is set to True, we may remove all points that are lexicographically earlier than the point added last
        last = qary_vector_to_decimal(newChoice,q)
        for index in range(last):
            record.append(index)
              
    # Removing all points that are marked for removal from the list of possible choices, 
    # storing the remaining choices in a new vector
    
    possibleChoices = sorted(set(possibleChoices).difference(record))

    return possibleChoices
    
def buildCap(mgenSet, possibleChoices, m):
    '''
    Function that recursively traverses the search space of all relevant t-element subsets of points in F_q^n and determines
    all m-general sets in affine space by determining which points may be added to a given set while maintaining 
    the m-general property of the set, starting from a set of n affinely independent points. 
    
    Input: set of points that form a m-general set,
           possible choices that can potentially be added to this set while maintaining m-general property,
           parameter m for m-general sets
    Output: void
    '''
    
    global maxSize_mgeneralSet
    global max_mgeneralSet
    global minSize_mgeneralSet
    global min_mgeneralSet

    # Base Case: There are no possible points to add without violating the m-general property
    if len(possibleChoices) == 0:
        
        if len(mgenSet) > maxSize_mgeneralSet and opt == "max":
            maxSize_mgeneralSet = len(mgenSet)
            max_mgeneralSet = mgenSet.copy()
            print(f"A complete {m}-general set of (max) size {maxSize_mgeneralSet} is given by {set_to_decimal(max_mgeneralSet,q)}.")
        
        if len(mgenSet) == maxSize_mgeneralSet:
            print(f"A complete {m}-general set of size {len(mgenSet)} is given by {set_to_decimal(mgenSet,q)}.")
            
        if len(mgenSet) < minSize_mgeneralSet and opt == "min":
            minSize_mgeneralSet = len(mgenSet)
            min_mgeneralSet = mgenSet.copy()
            print(f"A complete {m}-general set of (min) size {minSize_mgeneralSet} is given by {set_to_decimal(min_mgeneralSet,q)}.")
        
        if len(mgenSet) == minSize_mgeneralSet:
            print(f"A complete {m}-general set of size {len(mgenSet)} is given by {set_to_decimal(mgenSet,q)}.")

        return 
    
    # Recursion Step: Traverse the list of possible points to check
    else:
        for i in range (0,len(possibleChoices)): 
                    
            mgenSet.append(decimal_to_qary_vector(possibleChoices[i],n,q))
            newPossible = possibleChoices.copy()
            newPossible.remove(possibleChoices[i])
            
            possible = possibleChoices_fun(mgenSet,newPossible,decimal_to_qary_vector(possibleChoices[i],n,q),m)
            
            buildCap(mgenSet, possible,m)
            
            mgenSet.pop()
        
def main():
    '''
    Main function that initializes the parameters and starts the search for m-general sets in AG(n,q) 
    '''
    
    # The original general set (consisting of affinely independent vectors in F_q^n) must have
    # cardinality at least m-1, since the inequality q^n >= m-1 must be satisfied in order for the 
    # program to work (otherwise the for loops will not be triggered). Note however that this is
    # not a very restrictive constraint for small m, since if q^n < m-1 for m = 3,4,5 and some prime q,
    # then n is 1 or 2, and these cases are trivial
    
    if (n < m-2):
        print("Please choose n such that n >= m-2 is satisfied.")
    
    global maxSize_mgeneralSet # variable that captures the size of the largest complete m-general set found
    maxSize_mgeneralSet = 0
    global max_mgeneralSet # variable that captures the elements of the largest complete m-general set found
    max_mgeneralSet = []
    global minSize_mgeneralSet # variable that captures the size of the smallest complete m-general set found
    minSize_mgeneralSet = float("inf")
    global min_mgeneralSet # variable that captures the elements of the smallest complete m-general set found
    min_mgeneralSet = []

    global affine_comb_Coeff
    affine_comb_Coeff = coeff_vectors(m,q)
    
    if mode == "init":
    # Initializing the search with a set of n+1 points that is guaranteed to be m-general 
        start = [decimal_to_qary_vector(0,n,q)] 
        for i in range(n):
            start.append(decimal_to_qary_vector(q**i,n,q))
            
    elif mode == "zero":
    # Initializing the search with an empty set 
        start = []
    
    else:
        pass
    
    # Determining initial possible additions
    possible = list(range(q**n))
    possible = sorted ( set(possible).difference(set_to_decimal(start,q)) )
       
    possible = possibleChoices_fun(start, possible, decimal_to_qary_vector(0,n,q), m, lex = False) 
    for i in range(n): # checking all affine combinations of initial set of points
        possible = possibleChoices_fun(start, possible, decimal_to_qary_vector(q**(i),n,q), m, lex = False) 
           
    print(f"\n Starting the search for {m}-general sets over AG({n},{q})...\n")

    buildCap(start,possible,m)
    
    if opt == "max":
        print(f"\n A complete {m}-general set of maximal size {maxSize_mgeneralSet} in AG({n},{q}) is {set_to_decimal(max_mgeneralSet,q)}.")
    elif opt == "min":
        print(f"\n A complete {m}-general set of minimal size {minSize_mgeneralSet} in AG({n},{q}) is {set_to_decimal(min_mgeneralSet,q)}.")

if __name__ == '__main__':
    
    # Parsing terminal input
    parser = argparse.ArgumentParser(description="Finding m-general sets in AG(n,q)")

    parser.add_argument("-m", "--m", type=int, required = True, help="m for m-general set")
    parser.add_argument("-q","--q", type=int, required = True, help="q for order of the underlying field")
    parser.add_argument("-n","--n", type=int, required = True, help="dimension of the space over F_q")
    parser.add_argument("-mode","--mode", type=str, required = True, help="init or zero: mode of initialization")
    parser.add_argument("-opt","--opt", type=str, required = True, help="min or max: mode of optimization")

    args = parser.parse_args()

    global m # positive integer that defines the m in m-general sets
    m = args.m
    global q # positive prime integer that denotes the size of the underlying field
    q = args.q
    global n # positive integer that encodes the dimension of the space over F_q that we are finding m-general sets on 
    n = args.n
    global mode 
    mode = args.mode
    global opt
    opt = args.opt
    
    main()

