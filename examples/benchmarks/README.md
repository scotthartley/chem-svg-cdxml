# Benchmark prompts (M5)

`prompts.yaml` holds one-sentence requests a use session should turn into an approved drawing,
each with a reference spec that renders well. `tests/test_benchmarks.py` checks that every
reference exists and validates.

## Running a trial

1. Freeze the skill: `git archive HEAD | tar -x -C build/m5_trials/<round>/repo` (with `fonts/`),
   or the working tree (`git ls-files -co --exclude-standard | rsync -a --files-from=- . <repo>`).
   Remove `examples/` and `tests/` from the copy and start an empty `issues/inbox.md`.
2. For each prompt, start a fresh agent in `build/m5_trials/<round>/<id>/`. Its brief holds only
   the prompt, the path to the frozen `SKILL.md`, and the use-session rules (SKILL.md "Rules").
   It must not read `src/`, `examples/`, `tests/` or the gallery.
3. The agent appends friction to the frozen `issues/inbox.md` and writes `trial.md`:
   - the number of `validate.py` runs and render iterations,
   - each validator message it hit (verbatim),
   - what it guessed,
   - the final files (`.yaml`, `.svg`, `.png`, `.cdxml`).
4. Compare the round's counts with the previous round; triage the inbox items (CLAUDE.md).

A prompt whose answer is nearly a SKILL.md or schema example only tests copying (r2: synthesis,
polymers, polymerization, suzuki_cycle). Write new prompts unlike every doc example.

`build/` is not committed; summarize each round in `references/lessons.md`.
