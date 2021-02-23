import random

def Count(motifs):
    count = {}
    k = len(motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    t = len(motifs)
    for i in range(t):
        for j in range(k):
            symbol = motifs[i][j]
            count[symbol][j] += 1
    return count

def Profile(motifs):
    profile = {}
    k = len(motifs[0])
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
             profile[symbol].append(0)
    t = len(motifs)
    for i in range(t):
        for j in range(k):
            symbol = motifs[i][j]
            profile[symbol][j] += 1/t
    return profile

def Consensus(motifs):
    k = len(motifs[0])
    count = Count(motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

def Score(motifs):
    score = 0
    consensus = Consensus(motifs)
    k = len(motifs[0])
    t = len(motifs)
    for i in range(t):
        for j in range(k):
            symbol = motifs[i][j]
            if symbol != consensus[j]:
                score += 1
    return score

def Pr(text, profile):
    product = 1
    for index, nucleotide in enumerate(text):
        product *= profile[nucleotide][index]
    return product

def ProfileMostProbableKmer(text, k, profile):
    kmers = {}
    prob_kmer = []
    n = len(text)
    for i in range(n - k + 1):
        kmer = text[i:i + k]
        Pr_kmer = Pr(kmer, profile)
        kmers[kmer] = Pr_kmer
    max_kmer = max(kmers.values())
    for kmer in kmers:
        if kmers[kmer] == max_kmer:
            prob_kmer.append(kmer)
    return prob_kmer[0]

def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for m in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][m:m+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

def Entropy(profile):
    entropy = 0
    for i in range(len(profile)):
        for j in profile[i]:
            if j == 0:
                entropy += 0
            else:
                entropy += -(j*(math.log(j,2)))
    return entropy

def CountWithPseudocounts(motifs):
    count = {}
    k = len(motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1)
    t = len(motifs)
    for i in range(t):
        for j in range(k):
            symbol = motifs[i][j]
            count[symbol][j] += 1
    return count

def ProfileWithPseudocounts(motifs):
    k = len(motifs[0])
    profile = CountWithPseudocounts(motifs)
    for i in range(k):
        n_symbol = 0
        for symbol in "ACGT":
            n_symbol += profile[symbol][i]
        for symbol in "ACGT":
            profile[symbol][i] = profile[symbol][i]/n_symbol
    return profile

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for m in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][m:m+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

def Motifs(profile, Dna):
    motifs = []
    t = len(Dna)
    for i in range(t):
        motif = ProfileMostProbableKmer(Dna[i], k, profile)
        motifs.append(motif)
    return motifs

def RandomMotifs(Dna, k, t):
    motifs = []
    for i in range(t):
        start = random.randint(0, len(Dna[i]) - k)
        motif = Dna[i][start:start + k]
        motifs.append(motif)
    return motifs

def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs

def Normalize(probabilities):
    normalized_prob = {}
    k = sum(probabilities.values())
    for kmer in probabilities:
        normalized_prob[kmer] = probabilities[kmer] / k
    return normalized_prob

def WeightedDie(probabilities):
    p = random.uniform(0, 1)
    for kmer in probabilities:
        p -= probabilities[kmer]
        if p <= 0:
            return kmer

def ProfileGeneratedString(text, profile, k):
    n = len(text)
    probabilities = {}
    for i in range(0, n-k+1):
        probabilities[text[i:i+k]] = Pr(text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

def GibbsSampler(Dna, k, t, N):
    motifs = RandomMotifs(Dna, k, t)
    best_motifs = motifs
    for j in range(N):
        i = random.randint(0, t-1)
        motifs_i = motifs[:i] + motifs[i+1:]
        profile = ProfileWithPseudocounts(motifs_i)
        motifs[i] = ProfileGeneratedString(Dna[i], profile, k)
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs
        return best_motifs
