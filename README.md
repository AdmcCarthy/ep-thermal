# ep-thermal: How hot the rock is, from a bottom-hole reading to a formation temperature

Formation temperature: Horner extrapolation, gradients and the seabed temperature.

**Status: private, built public-ready.** The library follows the rules of a public one (no host-repository imports, no data, no private paths, an automated publication gate in `scripts/publication_gate.py`) and is published only when that gate passes. Licence: MIT with acknowledgement (`LICENSE`).

## Install

```
uv add git+https://github.com/AdmcCarthy/ep-thermal
```

Depends on: `ep-physics-base`. Optional extras:

none

## What is in it

| Module | What it does |
| :--- | :--- |
| `temperature` | Formation temperature from bottom-hole temperatures: seabed temperature, depth below seabed, Harrison and Horner corrections, geothermal gradients. |

Functions are pure: they take frames and parameters and return `(result, events)` where they change or flag data; nothing reads or writes files except the configuration loader.

## Settings and assumptions

Every setting has a packaged default in `src/ep_thermal/defaults/`. A project can override a file by providing `config/<same name>` in its working directory; the project file wins. A default is an assumption unless it says otherwise: its level (`placeholder`, `inferred`, `published`, `official`, `exact`) and a one-line basis travel with it, and a citation (author and year) is given where one exists. Libraries never import a knowledge store: a host that has one passes its values in as an overlay.

## Develop

```
uv sync --all-extras
uv run pytest            # 1 test files
uv run ruff check .
uv run python scripts/publication_gate.py
```

Tests use small synthetic data generated on the fly. No network, no real data.
