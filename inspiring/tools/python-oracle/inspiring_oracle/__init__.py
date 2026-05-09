"""InspiRING.Pack reference oracle.

Implements Algorithm 1 of the InsPIRe paper (eprint 2025/1352) verbatim from
the project SPEC.md, at tiny parameters (d in {8, 16}). Used as the byte-equal
correctness ground truth for the Rust crate.

Modules are added stage by stage; see ../../../.cursor/plans/phase_2_python_oracle_*.plan.md.
"""

from inspiring_oracle.params import RlweParams, ORACLE_TINY, ORACLE_SMALL

__all__ = ["RlweParams", "ORACLE_TINY", "ORACLE_SMALL"]
