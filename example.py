from UfukMSA import mafft

sequences = {
    "seq1": "ACGTACGT",
    "seq2": "ACGTTCGT",
    "seq3": "ACGT-CGT",
}

result = mafft(sequences)
for name, aligned in result.items():
    print(f"{name}: {aligned}")