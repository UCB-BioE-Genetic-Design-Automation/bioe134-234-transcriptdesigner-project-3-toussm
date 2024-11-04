class GCContentChecker:
    """
    Efficiently checks if the GC content of a DNA sequence is within an acceptable range.
    """

    def __init__(self, min_gc=40, max_gc=60):
        """
        Initializes the GC content checker with minimum and maximum GC content thresholds.

        Parameters:
            min_gc (int): Minimum acceptable GC content percentage.
            max_gc (int): Maximum acceptable GC content percentage.
        """
        self.min_gc = min_gc
        self.max_gc = max_gc

    def is_gc_content_acceptable(self, sequence):
        """
        Efficiently checks if the GC content of a given DNA sequence falls within the acceptable range.

        Parameters:
            sequence (str): DNA sequence to check.

        Returns:
            bool: True if GC content is within range, False otherwise.
            float: The GC content percentage.
        """
        length = len(sequence)
        if length == 0:
            return False, 0.0  # Empty sequence has undefined GC content

        # Efficient count of 'G' and 'C' bases
        gc_count = sequence.count('G') + sequence.count('C')
        gc_content = (gc_count / length) * 100

        # Check if GC content is within the acceptable range
        is_acceptable = self.min_gc <= gc_content <= self.max_gc
        return is_acceptable, gc_content