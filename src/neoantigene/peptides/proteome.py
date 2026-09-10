"""Reference proteome lookup.

Expects an Ensembl `Homo_sapiens.GRCh38.pep.all.fa[.gz]`, whose headers carry
both the protein and transcript identifiers.
"""

from __future__ import annotations

from pathlib import Path
from typing import Self

from ..io.fasta import read_fasta
from ..io.ids import strip_version

_SEPARATOR = "*"


class ProteomeIndex:
    def __init__(self) -> None:
        self._by_transcript: dict[str, str] = {}
        self._by_protein: dict[str, str] = {}
        self._by_symbol: dict[str, str] = {}
        self._concatenated: str | None = None

    @classmethod
    def from_fasta(cls, path: Path) -> Self:
        index = cls()
        for header, raw_sequence in read_fasta(Path(path)):
            sequence = raw_sequence.rstrip("*")
            if not sequence:
                continue
            tokens = header.split()
            protein_id = strip_version(tokens[0]) if tokens else ""
            attributes: dict[str, str] = {}
            for token in tokens[1:]:
                if ":" in token:
                    key, value = token.split(":", 1)
                    attributes.setdefault(key, value)
            if protein_id:
                index._by_protein[protein_id] = sequence
            transcript_id = strip_version(attributes.get("transcript"))
            _keep_longest(index._by_transcript, transcript_id, sequence)
            _keep_longest(index._by_symbol, attributes.get("gene_symbol"), sequence)
        return index

    def get(self, transcript: str | None = None, gene: str | None = None) -> str | None:
        key = strip_version(transcript)
        if key:
            sequence = self._by_transcript.get(key) or self._by_protein.get(key)
            if sequence:
                return sequence
        if gene:
            return self._by_symbol.get(gene)
        return None

    def contains_peptide(self, peptide: str) -> bool:
        """Exact-match check against the whole reference proteome.

        Backed by a lazily built concatenated sequence: substring search is
        cheap enough for the few hundred peptides that survive gating, and
        avoids holding a multi-gigabyte k-mer set in memory.
        """
        if not peptide:
            return False
        if self._concatenated is None:
            sequences = self._by_protein.values() or self._by_transcript.values()
            self._concatenated = _SEPARATOR.join(sequences)
        return peptide in self._concatenated

    def __len__(self) -> int:
        return len(self._by_transcript)


def _keep_longest(index: dict[str, str], key: str | None, sequence: str) -> None:
    if not key:
        return
    previous = index.get(key)
    if previous is None or len(sequence) > len(previous):
        index[key] = sequence
