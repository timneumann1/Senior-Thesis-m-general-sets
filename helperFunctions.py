"""
In this script, we provide auxiliary functions that can be used for the search of complete m-general sets
over affine space AG(n,q) and projective space PG(n,q).
@author: Tim Neumann
"""

def coeff_vectors(m,q):
    '''
    Function that determines all coefficient vectors such that the sum of some m-1 vectors
    satisfies the constraints imposed by affinity
    
    Input: parameters m (for m-general set) and q (order of the underlying field)
    Output: coefficient vector
    '''
    
    coeffVecs = []
    for decVec in range(q**(m-1)):
        qaryVec = decimal_to_m_minus_one_vector(decVec,m,q)
        if ( sum_vector_components_q(decimal_to_m_minus_one_vector(decVec,m,q),q) == 1 ):
            coeffVecs.append(qaryVec)
    return coeffVecs

def coeff_vector_C0(m,q):
    '''
    Function that determines all coefficient vectors such that the sum of some m vectors
    satisfies the constraints needed to test for affine independence
    
    Input: parameters m (for m-general set) and q (order of the underlying field)
    Output: coefficient vector
    '''
    
    coeffVecs = []
    for decVec in range(1, q**m): # we have to exclude the zero vector
        qaryVec = decimal_to_m_vector(decVec,m,q)
        if ( sum_vector_components_q(decimal_to_m_vector(decVec,m,q),q) == 0 ):
            coeffVecs.append(qaryVec)
    return coeffVecs

def createCoefficientVec(m,q):
    '''
    Function that determines the coefficient vectors required to test for linear independence of a set of size m.
    
    Input: parameter m for size of the and q for size of underlying field
    '''

    C_0 = []
    for decVec in range(1,q**m): # we start the loop at 1, since we are not interested if linear combinations sum to 0 trivially
        qaryVec = decimal_to_m_vector(decVec,m,q)
        C_0.append(qaryVec)
    return C_0

def sum_vector_components_q(vector,q):
    '''
    Function that return the sum of the components of a given vector modulo the order of the underlying field
    
    Input: vector and order of the field
    Output: sum of the components
    '''
    
    sum = 0
    for i in range(0,len(vector)):
        sum += vector[i]
    return sum % q

def multiply_scalar_into_vector(scalar, vector, q):
    '''
    Function that multiplies a scalar into each component of a vector over some underlying field
    
    Input: scalar and vector for the operation described above, order of the underlying field
    Output: resulting vector 
    '''
    
    result = vector.copy()
    for i in range(0,len(vector)):
        result[i] = (scalar*vector[i]) % q
    return result

# DISCLAIMER: affine_comb and linear_comb are identical, but they receive different names in order to 
#             clearly state the context in which they are used (the coefficient vector for affine combinations
#             must satisfy an additional constraint)

def affine_comb(coeffVect, vectors, q):
    ''' 
    Function that determines an affine combinations of a set of vectors, using a coefficient vector
    that controls that the required constraints imposed by affinity are obeyed.
    Recall: an affine combination is a linear combination of vectors in which the coefficients sum to 1
    
    Input: coefficient vector with entries that sum to 1 over the underlying field, set of vectors,
           order of the underlying field
    Output: affine combination of the vectors under multiplication of the provided coefficients
    '''
    
    assert len(coeffVect) == len(vectors)
    newVecs = []
    for i in range(len(vectors)):
        newVecs.append(multiply_scalar_into_vector(coeffVect[i],vectors[i],q))
    affineCombination = add_vectors_mod_q(newVecs,q)
    return affineCombination

def linear_comb(coeffVect, vectors, q):
    ''' 
    Function that determines the linear combinations of a set of vectors, using a coefficient vector
    
    Input: coefficient vector with entries not equally zero, set of vectors,
           order of the underlying field
    Output: linear combination of the vectors under multiplication of the provided coefficients
    '''
    
    assert len(coeffVect) == len(vectors)
    newVecs = []
    for i in range(len(vectors)):
        newVecs.append(multiply_scalar_into_vector(coeffVect[i],vectors[i],q))
    linearCombination = add_vectors_mod_q(newVecs,q)
    return linearCombination

def add_vectors_mod_q(vectors, q):
    '''
    Function that adds a set of vectors over some underlying field of order q
    
    Input: set of vectors and size of the underlying field 
    Output: sum of the vectors over the field
    '''
    
    result = vectors[0].copy()
    for i in range(0,len(vectors[0])):
        result[i] = 0
        for j in range(0,len(vectors)):
            result[i] += vectors[j][i]
        result[i] = (result[i]) % q
    return result

def decimal_to_m_minus_one_vector(decim, m, q):
    '''
    Function that transforms a given decimal integer to a vector of size m-1 over an underlying field of order q
    
    Input: decimal integer, integer parameter m and order of underlying field 
    Output: vector of size m-1 over underlying field
    '''
    
    qary = [None] * (m-1)
    quotient: int
    quotient = decim
    for i in range(1,m):
        qary[m-1-i] = int(quotient % q)
        quotient = int(quotient // q)
    return qary  

def decimal_to_m_vector(decim, m, q):
    '''
    Function that transforms a given decimal integer to a vector of size m over an underlying field of order q
    
    Input: decimal integer, integer parameter m and order of underlying field 
    Output: vector of size m+1 over underlying field
    '''
    qary = [None] * (m)
    quotient: int
    quotient = decim
    for i in range(1,m+1):
        qary[m-i] = int(quotient % q)
        quotient = int(quotient // q)
    return qary  

def decimal_to_m_plus_one_vector(decim, m, q):
    '''
    Function that transforms a given decimal integer to a vector of size m+1 over an underlying field of order q
    
    Input: decimal integer, integer parameter m and order of underlying field 
    Output: vector of size m+1 over underlying field
    '''
    qary = [None] * (m+1)
    quotient: int
    quotient = decim
    for i in range(1,m+2):
        qary[m+1-i] = int(quotient % q)
        quotient = int(quotient // q)
    return qary  

def decimal_to_qary_vector(decim,n,q):
    '''
    Function that transforms a given decimal integer to a vector of size n over an underlying field of order q
    
    Input: decimal integer, integer parameter n for dimension of the vector and order of underlying field 
    Output: vector of size n over underlying field
    '''
    
    qary = [None] * n
    quotient: int
    quotient = decim
    for i in range(1,n+1):
        qary[n-i] = int(quotient % q)
        quotient = int(quotient // q)
    return qary  

def set_to_decimal(set,q):
    '''
    Function that transforms a set of arbitrary size, containing vectors of size n over an underlying field of order q
    into a set of decimal integers
    
    Input: set of vectors and order of the underlying field
    Output: decimal integer
    '''
    setResult = []
    for i in range(len(set)):
        setResult.append(qary_vector_to_decimal(set[i],q))
    return setResult

def qary_vector_to_decimal(vector,q):
    '''
    Function that transforms a given vector of size n over an underlying field of order q into a decimal integer
    
    Input: vector and order of underlying field,  
    Output: decimal integer
    '''
    decimal = 0
    for i in range(len(vector)):
        decimal += vector[i]*q**(len(vector)-(i+1))
    return decimal 

def inverse_finite_field(a, q): 
    '''
    Function that determines the inverse of a field element a over a finite Galois field of prime-power order q.
    Note that in order for the inverse to be uniquely defined, we need gcd(a,q)=1, which is guaranteed for q being prime.
    
    Input: element a and order of the field q
    Output: inverse of the element a over the field
    '''
     
    for x in range(q):
        if (a * x) % q == 1:
            return x

