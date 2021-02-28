# Ref: (For sgt)
#https://towardsdatascience.com/sequence-embedding-for-clustering-and-classification-f816a66373fb
# Not finished

import pandas as pd
from sgt import SGT
import numpy as np

def main():

    # Read sample data (www.uniprot.org)
    # Data taken from: https://github.com/cran2367/sgt/tree/master/python

    protein_data = pd.read_csv('data/protein_classification.csv')
    print(protein_data)

    # Data has the following entries:
    # Entry,Entry name,Status,Protein names,Gene names,Organism,Length,Sequence,Function [CC],Features,Taxonomic lineage (all),Protein families
    # Preprocessing
    sequences_data = protein_data['Sequence']
    print(sequences_data)

    sequences_splitted = [list(x) for x in sequences_data]
    print(sequences_splitted[0])
    #print(sequences_splitted)

    corpus = sequences_splitted.copy()

    # Convert 'letters' into integer
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    mapping = {w: i for i, w in enumerate(letters)}
    print(mapping)
    print(mapping['A'])

    sequences_length = len(sequences_splitted)

    for sequence_index in range (sequences_length):
        sequence_length = len(corpus[sequence_index])
        for index in range (sequence_length):
            ix2pr = mapping[corpus[sequence_index][index]]
            corpus[sequence_index][index] = ix2pr

    # Convert it into id - sequence
    corpus = pd.DataFrame({'sequence' : sequences_splitted, 'id': enumerate(sequences_splitted)})
    print(corpus)

    # Create 'k-mers' corpus: (non-overlapping 3-grams)
    # TODO - index-3mer mapping
    # TODO - Corpus 3mer conversion


    # SGT embedding of the sequence
    generator = SGT(kappa=1,
           lengthsensitive=False,
           flatten=False,
           mode='default')

    embedding = generator.fit_transform(corpus)

if __name__ == '__main__':
    main()