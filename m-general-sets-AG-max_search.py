"""
In this script, we determine complete m-general sets over affine space AG(n,q).
@author: Tim Neumann
"""

from helperFunctions import *
import argparse
import sys
import re
import ast

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
      
    if lex == True: # if flag is set to True, we remove all points that are lexicographically earlier than the point added last
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

    # Base Case: There are no possible points to add without violating the m-general property
    if len(possibleChoices) == 0:
        
        if len(mgenSet) > maxSize_mgeneralSet:
            maxSize_mgeneralSet = len(mgenSet)
            max_mgeneralSet = mgenSet.copy()
            print(f"A complete {m}-general set of (max) size {maxSize_mgeneralSet} is given by {set_to_decimal(max_mgeneralSet,q)}.")
            
        if len(mgenSet) == maxSize_mgeneralSet:
            print(f"A complete {m}-general set of size {len(mgenSet)} is given by {mgenSet}.")

        return 
    
    # Recursion Step: Traverse the list of possible points to check
    else:
        for i in range (0,len(possibleChoices)): 
                    
            mgenSet.append(decimal_to_qary_vector(possibleChoices[i],n,q))
            possible = possibleChoices_fun(mgenSet,possibleChoices,decimal_to_qary_vector(possibleChoices[i],n,q),m)
    
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
    
    if (q**n < m-1):
        print("Please choose n such that q^n >= m-1 is satisfied.")
    
    global maxSize_mgeneralSet # variable that captures the size of the largest complete m-general set found
    maxSize_mgeneralSet = 0
    global max_mgeneralSet # variable that captures the elements of the largest complete m-general set found
    max_mgeneralSet = []

    global affine_comb_Coeff
    affine_comb_Coeff = coeff_vectors(m,q)
    
    start = []
    possible = list(range(q**n))

    if points == None:
        # THE BELOW SET MUST BE m-general FOR THE RESPECTIVE INPUT VARIABLES
        # e.g., for m = 4, q = 3, n = 7: [0, 28, 55, 90, 124, 148, 171, 202, 232, 258, 277, 322, 346, 362, 386, 425, 444, 465, 501, 538, 547, 587, 600, 633, 670]
        init_set = [] # Insert initialization set in decimal 
        
        for el in init_set:
            start.append(decimal_to_qary_vector(el,n,q))
       
    else: 
        for el in points:
            print(f"Element of initialization:{el}")
            print(el)
            start.append(el)
        
    print(f"Initialization set: {set_to_decimal(start,q)}")
    possible = sorted ( set(possible).difference(set_to_decimal(start,q)) )

    for el in start: # checking all affine combinations of initial set of points
        possible = possibleChoices_fun(start, possible, el, m, lex = False) 
         
    print(f"\n Starting the search for {m}-general sets over AG({n},{q})...\n")

    buildCap(start,possible,m)
    
    print(f"\n A complete {m}-general set of maximal size {maxSize_mgeneralSet} in AG({n},{q}) is {set_to_decimal(max_mgeneralSet,q)}.")

if __name__ == '__main__':
    
    # Terminal input
    parser = argparse.ArgumentParser(description="Finding m-general sets in AG(n,q) and PG(n,q)")

    parser.add_argument("-m", "--m", type=int, required = True, help="m for m-general set")
    parser.add_argument("-q","--q", type=int, required = True, help="q for order of the underlying field")
    parser.add_argument("-n","--n", type=int, required = True, help="dimension of the space over F_q")
    # Input points need to be given as list of n-tuples over F_q (in order to enable smooth compatibility with the output of maximal-sidon-set-F3.py)
    # e.g., -p "[0 0 0 0 0][0 1 0 1 0][0 1 0 2 0][0 1 1 0 1][0 2 0 1 1][0 2 2 2 1][0 1 1 0 2][0 2 2 1 2][0 2 0 2 2][1 0 0 0 0][1 1 1 0 0][1 1 2 0 0]"
    # Alternatively, initial points can be specified as list of integers in the main() function above.
    parser.add_argument("-p","--points", type=str, required = False, help="points in affine space AG(n,q) to check for m-general property")

    args = parser.parse_args()

    global m # positive integer that defines the m in m-general sets
    m = args.m
    global q # positive prime integer that denotes the size of the underlying field
    q = args.q
    global n # positive integer greater than 1 (we are not interested in projective lines) that encodes the dimension of the space over F_q that we are finding m-general sets on 
    n = args.n
    global points
    
    # Input processing
    if args.points:
        points = args.points
        pattern = r'\[(.*?)\]'
        matches = re.findall(pattern, points)
        points = [[int(num) for num in match.split()] for match in matches]

    else:
        points = None

    main()
    
