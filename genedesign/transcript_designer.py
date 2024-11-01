from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
import re
import random
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.seq_utils.Translate import Translate

class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and chooses an RBS using the highest CAI codon for each amino acid.
    """

    def __init__(self):
        self.aminoAcidToCodon = {}
        self.rbsChooser = None
        self.forbiddenSeqChecker = None
        self.codonChecker = None
        self.promoterChecker = None
        random.seed(10)

    def initiate(self) -> None:
        """
        Initializes the codon table, the RBS chooser and the checkers.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        self.forbiddenSeqChecker = ForbiddenSequenceChecker()
        self.forbiddenSeqChecker.initiate()

        self.codonChecker = CodonChecker()
        self.codonChecker.initiate()

        self.promoterChecker = PromoterChecker()
        self.promoterChecker.initiate()

        # Codon table with highest freq codon for each amino acid (for E. coli)
        
        with open('genedesign/data/codon_usage.txt', 'r') as file:

            for line in file:
                parts = line.split()
                codon = parts[0]
                amino_acid = parts[1]
                freq = float(parts[2])  # Extract frequency percentage as needed

                # Store codons with their frequencies in a temporary dictionary
                if amino_acid not in self.aminoAcidToCodon:
                    self.aminoAcidToCodon[amino_acid] = []
                self.aminoAcidToCodon[amino_acid].append((codon, freq))

            # # Calculate relative adaptive parameter and store in self.aminoAcidToCodon
            # for amino_acid, codons in codon_data.items():
            #     max_freq = max(freq for _, freq in codons)  # Find max frequency for the amino acid
            #     self.aminoAcidToCodon[amino_acid] = [(codon, freq / max_freq) for codon, freq in codons]

    def guided_random_codon_selection(self, aa):
        codons, weights = zip(*self.aminoAcidToCodon[aa])
        choice = random.choices(codons, weights=weights, k=1)[0]
        if weights[codons.index(choice)] >= 0.01:
            return choice
        else:
            return self.guided_random_codon_selection(aa)
    
    
    def codon_generator(self, amino_acids, samples=100):
        sampled_sequences = []
        for _ in range(samples):
            sampled_seq = ''.join(self.guided_random_codon_selection(aa) for aa in amino_acids)
            sampled_sequences.append(sampled_seq)
        return sampled_sequences
    
    def window_scorer(self, seq):
        score = 0
        if not self.forbiddenSeqChecker.run(seq)[0]:
            score += 10
        if not self.codonChecker.run(seq)[0]:
            score += 1
        if not hairpin_checker(seq)[0]:
            score += 1
        if not self.promoterChecker.run(seq):
            score += 1
        return score

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA and selects an RBS.
        
        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.
        
        Returns:
            Transcript: The transcript object with the selected RBS and translated codons.
        """
        peptide += '*'
        cds = '' #Initialize cds string

        #Iterate through the peptide string to 
    
        final_sequence = ''

        window_size = 3

        # Iterate over amino acid sequence in sliding windows
        for i in range(0, len(peptide), window_size):
            # Define upstream, current window, and downstream context
            preamble = final_sequence[:36]
            current_window = peptide[i:i + window_size]
            downstream_context = peptide[i + window_size:i + window_size + 6]

            # Use codon_generator to sample codons for the current window
            sampled_codons = self.codon_generator(current_window + downstream_context)

            # Score each sampled codon sequence
            scored_codons = [(seq[:9], self.window_scorer(preamble + seq)) for seq in sampled_codons]

            # Choose the best codon sequence based on score
            best_codon_seq = min(scored_codons, key=lambda x: x[1])[0]

            # Retain only the middle 3 codons and add to final sequence
            final_sequence += best_codon_seq

        cds += final_sequence

            # Choose an RBS
        selectedRBS = self.rbsChooser.run(cds, ignores)

        # Return the Transcript object
        return Transcript(selectedRBS, peptide, codons=[cds[i:i+3] for i in range(0, len(cds), 3)])
    

if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MYPFIRTARMTV"
    
    designer = TranscriptDesigner()
    designer.initiate()
    translate = Translate()
    translate.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)
    
    # Print out the transcript information
    print(transcript)

    cds = ''.join(transcript.codons)
    print(translate.run(cds))

