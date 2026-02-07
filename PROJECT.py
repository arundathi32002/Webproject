#Biopython project - Arundathi. A 

#STEP 1. RETRIEVAL OF SEQUENCES
# UNIPROT- FASTA files of hypothetical proteins from five model organisms (Human, Rat, Roundworm, Rice, and Thale cress) were downloaded.

Seq1 = "Human - B3KN02"
Seq2 = "Rat - M0RCA7"
Seq3 = "Roundworm - P34459"
Seq4 = "Rice - Q60EB8 "
Seq5 = "Thale cress - A0A7G2EZ38"
print("Hypothetical protein FASTA files from five different model organisms were downloaded from UniProt and are listed below:")
print(Seq1)
print(Seq2)
print(Seq3)
print(Seq4)
print(Seq5)

#STEP 2. Sequence Quality Analysis

#2.1 First import the FASTA files

from Bio import SeqIO

fasta_files = [
    "HUMAN.fasta",
    "RAT.fasta",
    "ROUNDWORM.fasta",
    "RICE.fasta",
    "THALE_CRESS.fasta"
]
print("-" *100)

#2.2 Quality parameters 

from collections import Counter

for file in fasta_files:
    print("File:", file)
    print("-" *50)
    for record in SeqIO.parse(file, "fasta"):
        print("QUALITY PARAMETES")
        print("ID:", record.id)
        print("Description:", record.description)
        print("Sequence:", record.seq)
        print("Length:", len(record.seq))
        print("Amino acid composition:",Counter(record.seq)) 
       

#STEP 3. Filtering and validation 
#3.1 Basic criteria 

        valid_aa = set("ARNDCQEGHILKMFPSTWYV")
        min_length = 100
        min_valid_percent = 95
        seq = str(record.seq).upper()
        length = len(seq)
        valid_count = sum(1 for aa in seq if aa in valid_aa)
        valid_percent = (valid_count / length) * 100
        print("Valid amino acids (%):", f"{valid_percent:.2f}")
        
#3.2 Decision 
        if length < min_length:
            print("Decision: Rejected - sequence too short")
        elif valid_percent < min_valid_percent:
            print("Decision: Rejected - low sequence quality")
        else:
            print("Decision: Accepted - Suitable for downstream analysis ")
        print("-" *100)


#STEP 5. Functional annotation
 
        from Bio.Blast import NCBIXML
        # List of BLAST result files
        blast_files = [
            "HUMAN_result.xml",
            "RAT_result.xml",
            "ROUNDWORM_result.xml",
            "RICE_result.xml",
            "THALE_CRESS_result.xml"
        ]

        for file in blast_files:
            print("\nBLAST Result File:", file)
            print("-" * 50)

            with open(file) as b:
                blast_record = NCBIXML.read(b)

            if blast_record.alignments:
                best_alignment = blast_record.alignments[0]   # Best hit
                best_hsp = best_alignment.hsps[0]              # Best HSP

                print("Best Hit Title:")
                print(best_alignment.title)
                print("Alignment Length:", best_alignment.length)
                print("Bit Score:", best_hsp.score)
                print("E-value:", best_hsp.expect)

                print("\nAlignment Details:")
                print("Query  :", best_hsp.query)
                print("Match  :", best_hsp.match)
                print("Subject:", best_hsp.sbjct)
            else:
                print("No significant BLAST hits found")
                print("-" * 100)