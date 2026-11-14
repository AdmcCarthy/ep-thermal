# Status of ep-thermal

**State.** Cut out of the E&P data engine on the date of the first commit: 1 modules, 1 test files. The engine switched over on 2026-11-11: it depends on this library (its old module paths are aliases of the modules here), and this repository is the source of truth.

## Decisions (Adam, 2026-11-11)
* Domain libraries, one per topic; this one is `ep-thermal`. Figures belong to their domain.
* Built public-ready; made public on 2026-11-11 after the gate passed (0 blockers). Licence MIT with acknowledgement.
* One-way design: the library never imports a knowledge store; assumptions travel as packaged defaults with a level, a one-line basis and an author-and-year citation where there is one.
* The hydrocarbon colours of the Sodir Geostandard are used as they are.

## Changes made when the code was cut out
* none.

## Tests that moved only in part (the rest use the engine and stay there)
* `test_temperature.py: kept 10, 1 stay in the engine`

## Couplings still to cut
* none: every `src` import was mapped.

## Next steps
* Add the knowledge identifiers used here to `KNOWLEDGE.md` and export the defaults from the knowledge store.
* Keep `scripts/publication_gate.py` at zero blockers (it runs in CI); read the warnings (Sodir, NPD, Equinor Design System names) when adding documents.

## Conductivity, quartz cement and maturity (2026-11-14)
* Three modules and three default files added on the owner's decision that the values are valid for use for now: `conductivity`, `quartz_cement`, `maturity` (38 tests in all). Details and limits in `docs/conductivity_cement_maturity.md`.
* To verify when data are available: conductivity against measured values, the cement band against core cement with a known temperature history, EASY%Ro against its printed table and measured Ro. The EASY%Ro constants were transcribed from an open-source code; the paper has not been read.
* Not released: no version bump or tag yet (changelog entry is under Unreleased).
