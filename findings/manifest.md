---
schema_version: 1
generated: 2026-09-22
next_id: 3
engine_version: 0.5.1
---

# Findings manifest

Derived index of the findings graph — regenerable from the finding files. One row per finding.

| ID | Slug | Title | Status | Phase | Kind | Entities | Relationships | Updated | Data version |
|----|------|-------|--------|-------|------|----------|---------------|---------|--------------|
| 1 | run-order-aliased-with-condition | Run order is aliased with condition: every control was run before every raloxifene-d0 sample within each batch | candidate | exploratory | caveat | — | relates_to:2 | 2026-09-22 | sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74 |
| 2 | two-acquisition-batches-6-vs-2 | Unequal two-batch structure of unknown nature: 6 runs on 2021-05-06 vs 2 runs on 2022-03-18 | candidate | exploratory | caveat | — | relates_to:1 | 2026-09-22 | sha256:bc6b73d30e40d6ec190f8cd4494ec58d984a475f670d5aab5d8b238974d1ba74 |
