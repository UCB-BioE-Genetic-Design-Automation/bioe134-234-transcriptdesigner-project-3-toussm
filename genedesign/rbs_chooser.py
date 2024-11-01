from genedesign.models.rbs_option import RBSOption
from genedesign.seq_utils.Translate import Translate
from genedesign.seq_utils.hairpin_counter import hairpin_counter
from genedesign.seq_utils.calc_edit_distance import calculate_edit_distance
from typing import Set
from Bio import SeqIO
from collections import defaultdict
import pandas as pd

class RBSChooser:
    """
    A class to choose the best RBS for a given CDS sequence.
    """
    
    rbs_options: Set[RBSOption] = set()

    def initiate(self) -> None:
        """
        Initialization method for RBSChooser.
        """
        
        # Function to extract UTR, gene, and CDS information from the GenBank file
        def extract_genes_info(genbank_file):
            gene_dict = defaultdict(dict)  # Dictionary to store gene info
            for record in SeqIO.parse(genbank_file, "genbank"):
                for feature in record.features:
                    if feature.type == "gene":
                        locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                        gene_name = feature.qualifiers.get("gene", [None])[0]

                        # CDS information
                        cds_feature = None
                        for cds in record.features:
                            if cds.type == "CDS" and cds.qualifiers.get("locus_tag") == [locus_tag]:
                                cds_feature = cds
                                break

                        if cds_feature:
                            start, end = cds_feature.location.start, cds_feature.location.end
                            strand = cds_feature.location.strand
                            if strand == 1:  # Forward strand
                                utr_start = max(0, start - 50)
                                utr_seq = record.seq[utr_start:start]
                            else:  # Reverse strand, we need to reverse complement
                                utr_start = end
                                utr_seq = record.seq[utr_start:utr_start + 50].reverse_complement()

                            cds_seq = cds_feature.extract(record.seq)
                            # Save the gene information in the dictionary
                            gene_dict[locus_tag] = {
                                "gene": gene_name,
                                "UTR": utr_seq,
                                "CDS": cds_seq
                            }
            return gene_dict

        # Example usage
        genbank_file = "genedesign/data/sequence.gb"  # Now using the correct file name
        genes_info = extract_genes_info(genbank_file)

        # Function to prune the dataset and return top 5% abundance values as a dictionary
        def prune_proteomics_data(file_path):
            # Read the dataset as a pandas DataFrame, skipping metadata lines starting with '#'
            data = pd.read_csv(file_path, comment='#', sep='\t', names=['locus_tag', 'abundance'])

            # Remove the '511145.' prefix from the locus tags
            data['locus_tag'] = data['locus_tag'].str.replace('511145.', '', regex=False)

            # Sort the data by abundance in descending order
            data_sorted = data.sort_values(by='abundance', ascending=False)

            # Calculate the top 5% of the dataset
            top_5_percent_count = int(len(data_sorted) * 0.05)

            # Prune the dataset to the top 5%
            top_5_percent_data = data_sorted.head(top_5_percent_count)

            # Convert the pruned data into a dictionary with locus_tag as key and abundance as value
            pruned_data_dict = dict(zip(top_5_percent_data['locus_tag'], top_5_percent_data['abundance']))

            return pruned_data_dict

        # Path to your dataset file
        file_path = 'genedesign/data/511145-WHOLE_ORGANISM-integrated.txt'

        # Call the function and get the pruned data dictionary
        pruned_data_dict = prune_proteomics_data(file_path)

        # Function to merge the top 5% proteomics data with gene sequence data
        def merge_data(top_5_dict, gene_dict):
            merged_dict = {}
            for locus_tag in top_5_dict:
                if locus_tag in gene_dict:
                    merged_dict[locus_tag] = {
                        "abundance": top_5_dict[locus_tag],
                        "gene_info": gene_dict[locus_tag]
                    }
            return merged_dict


        # Merge the two datasets
        merged_data = merge_data(pruned_data_dict, genes_info)

        translate = Translate()
        translate.initiate()

        for locus_tag, info in merged_data.items():
            rbs_option = RBSOption(
                utr=info['gene_info']['UTR'],
                cds=info['gene_info']['CDS'],
                gene_name=info['gene_info']['gene'],
                first_six_aas=translate.run(info['gene_info']['CDS'][:18])
            )

            self.rbs_options.add(rbs_option)


    def run(self, cds: str, ignores: Set[RBSOption]) -> RBSOption:
        """
        Chooses the best RBS option for the provided CDS sequence, ignoring specified options and refining based on secondary structure and peptide similarity.

        Parameters:
            cds (str): The coding sequence for which to choose the best RBS.
            ignores (Set[RBSOption]): A set of RBSOptions to ignore during selection.

        Returns:
            RBSOption: The best RBS option for the given CDS.
        """
        translate = Translate()
        translate.initiate()

        peptide_cds = translate.run(cds[:18])  # Get the first six amino acids from the input CDS

        # Step 2: Filter out RBS options that are in the ignores set
        valid_rbs_options = [rbs for rbs in self.rbs_options if rbs not in ignores]

        # Step 3 & 4: Score each valid RBS option based on both secondary structure and peptide similarity
        scored_rbs_options = []
        for rbs in valid_rbs_options:
            # Recalculate potential hairpins between the RBS and the input CDS
            combined_sequence = rbs.utr + cds
            hairpin_count = hairpin_counter(combined_sequence)[0]

            # Calculate peptide similarity (edit distance)
            edit_distance = calculate_edit_distance(peptide_cds, rbs.first_six_aas)

            # Score calculation: minimize both hairpin count and edit distance
            # Assign weights to hairpin count and edit distance

            score = (hairpin_count * 2) + (edit_distance * 1.5)

            # Store the RBS and its score
            scored_rbs_options.append((rbs, score))

        # Sort by score (ascending, lower score is better)
        scored_rbs_options.sort(key=lambda x: x[1])

        # Return the best RBS based on the combined score
        if scored_rbs_options:
            return scored_rbs_options[0][0]  # Return the RBSOption with the lowest score
        else:
            return None  # Return None if no valid RBSOption found

if __name__ == "__main__":
    # Example usage of RBSChooser
    cds = "ATGGTAAGAAAACAGTTGCAGAGAGTTGAATT"

    # Initialize the chooser
    chooser = RBSChooser()
    chooser.initiate()

    # Choose RBS with no ignores
    ignores = set()
    selected1 = chooser.run(cds, ignores)
    
    # Add the first selection to the ignore list
    ignores.add(selected1)
    
    # Choose another RBS option after ignoring the first
    selected2 = chooser.run(cds, ignores)

    # Print the selected RBS options
    print("Selected1:", selected1)
    print("Selected2:", selected2)
