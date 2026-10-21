# EPSR_python_translation

A line-by-line Python translation of EPSR26 (Empirical Potential Structure Refinement,
A. K. Soper), written to be read and modified rather than to be fast. It keeps the
`epsrshell` interface and is validated byte for byte against the Fortran on measured
neutron data for water, amorphous silica and Cu-BTC.

Status: setting up. See [PLAN.md](PLAN.md).

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e '.[dev]'
```
