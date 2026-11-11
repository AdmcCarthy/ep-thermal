# Status of ep-thermal

**State.** Cut out of the E&P data engine on the date of the first commit: 1 modules, 1 test files. The engine switched over on 2026-11-11: it depends on this library (its old module paths are aliases of the modules here), and this repository is the source of truth.

## Decisions (Adam, 2026-11-11)
* Domain libraries, one per topic; this one is `ep-thermal`. Figures belong to their domain.
* Private until the publication gate passes; built public-ready. Licence MIT with acknowledgement.
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
* Run `scripts/publication_gate.py` and clear the warnings before publishing.
