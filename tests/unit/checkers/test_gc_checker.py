import pytest
from genedesign.checkers.gc_checker import GCContentChecker

@pytest.fixture
def gc_checker():
    checker = GCContentChecker()
    return checker

def test_valid_gc_content(gc_checker):
    # Sequences with GC content within 40% - 60%
    valid_gc_seqs = [
        "ATGCGCATATGCGCGTATAT",  # ~50% GC
        "GCGTATGCGCGCATATATAA",  # ~55% GC
        "ATGCATGCATGCATGCATGC",  # 50% GC
        "ATGCGCATGCGTATGCATAT",  # ~50% GC
        "GCGCGTATATCGCGCGTATG",  # ~55% GC
    ]

    print(">> Testing sequences with GC content between 40% and 60% (expected True)")
    for seq in valid_gc_seqs:
        result, content = gc_checker.is_gc_content_acceptable(seq)
        print(f"result: {result} on {seq} (GC content: {content:.2f}%)")
        assert result == True

def test_invalid_gc_content(gc_checker):
    # Sequences with GC content outside 40% - 60%
    invalid_gc_seqs = [
        "ATATATATATATATATATAT",  # ~0% GC
        "CCCCCCCCCCCCCCCCCCCC",  # 100% GC
        "GCGCGCGCGCGCGCGCGCGC",  # 100% GC
        "ATATATATGCATATATATAA",  # ~20% GC
        "CGCGCGCGCGCGCGCGATGC",  # ~90% GC
    ]

    print("\n>> Testing sequences with GC content outside 40% - 60% (expected False)")
    for seq in invalid_gc_seqs:
        result, content = gc_checker.is_gc_content_acceptable(seq)
        print(f"result: {result} on {seq} (GC content: {content:.2f}%)")
        assert result == False