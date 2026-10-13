# Research notes

This folder is the record of what we have actually checked. It is written in
plain language so it can be read without the chat that produced it.

Markdown is the format because it stays in the repository, opens in any
editor, and can be searched. Each file is one subject. Start here, then open
the file for the question you have.

The notes distinguish three kinds of statement:

- **Counted.** We opened the file or ran the comparison.
- **Reported.** A paper or a company said it. We name who.
- **Not shown.** A number people quote that this record does not support.

Nothing here changes the ranking rule in `config/default.yaml`. That rule
stays frozen until a later group of patients is used to train, and a group
after that is used to test.

## Files

| File | Question it answers |
| --- | --- |
| [the-study.md](the-study.md) | What the first study is, what it can tell, and how the company would get paid |
| [moderna.md](moderna.md) | What Moderna showed, and which of their lists are public |
| [public-exams.md](public-exams.md) | Which public patient sets can grade a picker, and which cannot |
| [field-audit.md](field-audit.md) | What the multi-agent audit checked in other people's data and code |

Longer lab notes, with file paths and line references, stay in `docs/`:

- `docs/AUDIT.md` is the source record for the multi-agent audit.
- `docs/BENCHMARK.md` is the Ott and CU04 runs.
- `docs/HARNESS.md` is how a comparison is scored.
- `docs/CITATIONS.md` is why each weight in the frozen rule has the sign it has.

`share/index.html` is the one-page explanation for a co-founder or an
investor. It is a summary of the study, not this record.
