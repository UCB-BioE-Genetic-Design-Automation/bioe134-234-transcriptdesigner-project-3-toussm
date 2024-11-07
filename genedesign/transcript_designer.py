from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
import random
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.seq_utils.Translate import Translate
from genedesign.checkers.gc_checker import GCContentChecker

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
        self.selectedRBS = None
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

        self.gcChecker = GCContentChecker().is_gc_content_acceptable

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
                if codon in self.codonChecker.rare_codons:
                    freq = freq / 5 # Lower freq value of rarecodons so they are much less likely to be chosen
                    self.aminoAcidToCodon[amino_acid].append((codon, freq))
                else:
                    self.aminoAcidToCodon[amino_acid].append((codon, freq))


    def guided_random_codon_selection(self, aa):
        codons, weights = zip(*self.aminoAcidToCodon[aa])
        return random.choices(codons, weights=weights, k=1)[0]
    
    def window_scorer(self, seq):
        score = 0
        if not self.forbiddenSeqChecker.run(seq)[0]:
            score += 15
        if not self.codonChecker.run(seq)[0]:
            score += 5
        if not self.gcChecker(seq)[0]:
            score += 5
        if not hairpin_checker(seq)[0]:
            score += 5
        if not self.promoterChecker.run(seq):
            score += 5
        return score
    
    def codon_generator(self, amino_acids, preamble, samples=50):
        sampled_sequence = [('', 50)]
        for _ in range(samples):
            sampled_seq = ''.join(self.guided_random_codon_selection(aa) for aa in amino_acids)
            score = self.window_scorer(preamble + sampled_seq)
            if score < 5:
                return sampled_seq
            elif score < sampled_sequence[0][1]:
                sampled_sequence.insert(0, (sampled_seq, score))
        return sampled_sequence[0][0]
    
    

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA and selects an RBS.
        
        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.
        
        Returns:
            Transcript: The transcript object with the selected RBS and translated codons.
        """
        # Add stop codon
        peptide += '*'
    
        final_sequence = ''

        window_size = 3
        
        # Iterate over amino acid sequence in sliding windows
        for i in range(0, len(peptide), window_size):
            # Define upstream, current window, and downstream context
            preamble = final_sequence[-9:]
            current_window = peptide[i:i + window_size]
            downstream_context = peptide[i + window_size:i + window_size + 6]

            # Sample codons for the current window plus downstream context
            sampled_codon = self.codon_generator(current_window + downstream_context, preamble=preamble)

            # Add best window is added to final sequence
            final_sequence += sampled_codon[:9]
        
        # Choose an RBS    
        self.selectedRBS = self.rbsChooser.run(final_sequence, ignores)


        # Return the Transcript object
        return Transcript(self.selectedRBS, peptide, codons=[final_sequence[i:i+3] for i in range(0, len(final_sequence), 3)])
    

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

