def cluster_sequence(alg,collection):
    
    algs = ["antikt","eekt"] 
    algo = algs[alg]

    global Clustering
    global FlavourHelper

    if algo==0:
        Clustering = ExclusiveJetClusteringHelper(collection, njets, tag)