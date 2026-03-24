# Phase 10 Full Space-Time Framework

Phase 10 extends the Phase 9 minimal differentiable loop into a fuller
framework with:

- a full objective layer,
- dense and sparse backward differentiation,
- `xi=(q_bar,tau)` reparameterized optimization,
- benchmark and ablation runners,
- a unified framework suite.

Main entry points:

- `phase_10/blom_full_backward_diff.py`
- `phase_10/blom_space_time_optimizer.py`
- `phase_10/blom_benchmark_suite.py`
- `phase_10/blom_phase10_framework_suite.py`

Recommended commands:

```bash
python3 -m phase_10.blom_full_backward_diff
python3 -m phase_10.blom_space_time_optimizer
python3 -m phase_10.blom_benchmark_suite
python3 -m phase_10.blom_phase10_framework_suite
python3 -m unittest discover -s phase_10 -p 'test*.py'
```

Scope note:

- This framework is intended to be benchmarkable and optimizer-ready.
- General even `k` is exposed through the public interfaces, but explicit
  placeholders are still allowed wherever the current implementation cannot yet
  support a requested `k` faithfully.
