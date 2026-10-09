# Retired E07 test evidence

These text files preserve the original strict-xfail test and nine Phase 2F-6
passing diagnostic cases. They are not collected by pytest and are not executable
supported transport code. The implementation was removed in Phase 2F-7.

At master `62a51f9005f0cf1541f5c7996123b6b077edf847`, the original E07 test
fails with builtin-object subscripting. The nine diagnostic cases pass by
observing known broken behavior after test-only bypasses. Their outputs and
mathematical interpretation are recorded in `../PHASE2F6.md`. Git history
preserves the original implementation and complete runnable tests.

Retirement deliberately removes one expected-failure case and nine diagnostic
cases from collection; it does not convert them to passing mathematical
regressions or skips. `tests/test_core_packaging.py` checks that the retired
module cannot be imported or shipped, while compatibility and core math tests
continue to cover the supported functionality.
