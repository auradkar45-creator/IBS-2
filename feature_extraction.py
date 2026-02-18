import pandas as pd
from Bio import SeqIO
from Bio import pairwise2
from collections import Counter
from itertools import product
from sklearn.model_selection import train_test_split


#Files

HEALTHY_FASTA = "rice_1000_full_genes.fasta"
DISEASED_FASTA = "rice_diseased_clean.fasta"

FULL_DATASET_CSV = "rice_full_dataset.csv"
TRAIN_FEATURES = "train_features.csv"
TEST_FEATURES = "test_features.csv"
ALIGNMENT_OUTPUT = "alignment_example.txt"

#K Value for k-mer analysis
K = 3


#Create Dataset

rows = []

healthy_records = list(SeqIO.parse(HEALTHY_FASTA, "fasta"))
diseased_records = list(SeqIO.parse(DISEASED_FASTA, "fasta"))

for rec in healthy_records:
    rows.append([str(rec.seq), 0])

for rec in diseased_records:
    rows.append([str(rec.seq), 1])

df = pd.DataFrame(rows, columns=["DNA_SEQUENCE", "LABEL"])
df.to_csv(FULL_DATASET_CSV, index=False)

print("Dataset Created")
print("Total Samples:", len(df))
print("Healthy Samples:", len(healthy_records))
print("Diseased Samples:", len(diseased_records))


#Redundancy Check
unique_count = len(set(df["DNA_SEQUENCE"]))
duplicates = len(df) - unique_count

print("Unique Sequences:", unique_count)
print("Duplicate Sequences:", duplicates)

if duplicates == 0:
    print("Result: No redundancy detected.")
else:
    print("Warning: Redundant sequences found.")


#Feature Extraction

alphabet = ['A','T','G','C']
all_kmers = [''.join(p) for p in product(alphabet, repeat=K)]

def kmer_features(seq):
    counts = Counter(seq[i:i+K] for i in range(len(seq)-K+1))
    total = sum(counts.values())
    return [counts[k]/total if total > 0 else 0 for k in all_kmers]

def gc_content(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)

features = []
for seq in df["DNA_SEQUENCE"]:
    features.append(kmer_features(seq) + [gc_content(seq)])

feature_df = pd.DataFrame(features, columns=all_kmers + ["GC_Content"])
feature_df["LABEL"] = df["LABEL"]

print("Total Features per Sequence:", len(all_kmers) + 1)


#Splitting the train, test data
train_df, test_df = train_test_split(
    feature_df,
    test_size=0.2,
    stratify=feature_df["LABEL"],
    random_state=42
)

train_df.to_csv(TRAIN_FEATURES, index=False)
test_df.to_csv(TEST_FEATURES, index=False)

print("\nTrain Samples:", len(train_df))
print("Test Samples:", len(test_df))


#Alignment Checking
if len(healthy_records) > 0 and len(diseased_records) > 0:

    seq1 = str(healthy_records[0].seq)
    seq2 = str(diseased_records[0].seq)

    alignment = pairwise2.align.globalxx(seq1, seq2)[0]

    # Calculating mutation count
    matches = alignment[2]
    length = len(seq1)
    differences = length - matches

    with open(ALIGNMENT_OUTPUT, "w") as f:
        f.write(str(alignment))

    print("Alignment saved to:", ALIGNMENT_OUTPUT)
    print("Sequence Length:", length)
    print("Matching Positions:", matches)
    print("Differing Positions:", differences)

else:
    print("Alignment skipped, FASTA files missing.")

print("Completed!")
