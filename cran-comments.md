## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new release.

## Notes

* On macOS with a very recent Apple clang (21), a compiler warning originating
  from R's own header (`R_ext/Boolean.h`: unknown warning group
  `-Wfixed-enum-extension`) may be reported. It is unrelated to package code and
  does not occur on the standard CRAN toolchains.
