def PatternCount(text, pattern):
    count = 0
    for i in range(len(text) - len(pattern) + 1):
        if text[i:i + len(pattern)] == pattern:
            count = count + 1
    return count

def FrequencyMap(text, k):
    freq = {}
    n = len(text)
    for i in range(n - k + 1):
        pattern = text[i:i + k]
        freq[pattern] = 0
        for j in range(n - k + 1):
            if text[j:j + k] == pattern:
                freq[pattern] = freq[pattern] + 1
    return freq

def FrequentWords(text, k):
    words = []
    freq = FrequencyMap(text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words

def Reverse(pattern):
    rev = ""
    for char in pattern:
        rev = char + rev
    return rev

def Complement(pattern):
    complement = ""
    dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    for char in pattern:
        complement = complement + dict[char]
    return complement

def ReverseComplement(pattern):
    rev_pattern = Reverse(pattern)
    comp_pattern = Complement(rev_pattern)
    return comp_pattern

def PatternMatching(pattern, genome):
    positions = []
    k = len(pattern)
    n = len(genome)
    for i in range(n - k + 1):
        genome_pattern = genome[i:i + k]
        if genome_pattern == pattern:
            positions.append(i)
    return positions

def SymbolArray(genome, symbol):
    array = {}
    n = len(genome)
    ExtendedGenome = genome + genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(ExtendedGenome[i:i+(n//2)], symbol)
    return array

def FasterSymbolArray(genome, symbol):
    array = {}
    n = len(genome)
    ExtendedGenome = genome + genome[0:n//2]
    # look at the first half of genome to compute first array value
    array[0] = PatternCount(genome[0:n//2], symbol)
    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]
        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

def SkewArray(genome):
    skew = [0]
    n = len(genome)
    for i in range(1, n + 1):
        skew.append(skew[i - 1])
        if genome[i - 1] == "C":
            skew[i] = skew[i] - 1
        if genome[i - 1] == "G":
            skew[i] = skew[i] + 1
        else:
            skew[i] = skew[i]
    return skew

def MinimumSkew(genome):
    # generate an empty list positions
    pos = []
    # set a variable equal to SkewArray(Genome)
    skew = SkewArray(genome)
    # find the minimum value of all values in the skew array
    skew_min = min(skew)
    # range over the length of the skew array and add all positions achieving the min to positions
    for i in range(len(skew) - 1):
        if skew[i] == skew_min:
            pos.append(i)
    return pos

def HammingDistance(p, q):
    hamming = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            hamming += 1
        else:
            hamming = hamming
    return hamming

def ApproximatePatternMatching(text, pattern, d):
    pos = [] # initializing list of positions
    for i in range(0, (len(text) - len(pattern)) + 1):
        hamming = HammingDistance(pattern, text[i:i+len(pattern)])
        if hamming <= d:
            pos.append(i)
    return pos

def ApproximatePatternCount(pattern, text, d):
    count = 0 # initializing count of positions
    for i in range(0, (len(text) - len(pattern)) + 1):
        hamming = HammingDistance(pattern, text[i:i+len(pattern)])
        if hamming <= d:
            count += 1
    return count
