/**
 * A mutant peptide with its mutated residue marked, and the wild-type below.
 *
 * Showing the pair stacked is the point: the whole agretopicity argument is
 * "this residue changed and the binding changed with it", and that is only
 * legible when the two sequences line up character for character.
 */
export function MutantPeptide({
  mutant,
  wildtype,
  position,
}: {
  mutant: string;
  wildtype: string | null;
  position: number;
}) {
  // `mutation_position` is 1-based over the peptide.
  const index = position - 1;

  return (
    <div className="font-mono text-[13px] leading-tight">
      <div className="tracking-[0.08em]">
        {mutant.split("").map((residue, i) => (
          <span
            key={i}
            className={
              i === index ? "font-semibold text-primary" : undefined
            }
          >
            {residue}
          </span>
        ))}
      </div>
      {wildtype && (
        <div className="tracking-[0.08em] text-muted-foreground">
          {wildtype.split("").map((residue, i) => (
            <span key={i} className={i === index ? "underline" : undefined}>
              {residue}
            </span>
          ))}
        </div>
      )}
    </div>
  );
}
