# Validation record

Validated locally on 2026-10-16 using C++17 / Apple Clang and MATLAB R2024a.

- `make all test test-multipath`: all four C++ study executables compiled; 4,717 OFDM checks and 892 multipath checks passed.
- `make sanitize sanitize-multipath`: the same existing C++ suites passed AddressSanitizer and UndefinedBehaviorSanitizer checks.
- `./bin/ofdm smoke 262144 64`: 18 reduced-bit C++ BER points completed, with correct timing acquisition in every point; counts and BER arithmetic checked.
- `unpack_evidence.py` / `summarize_results.py`: 132,000 coarse and 396,000 fine timing estimates summarized into 66 and 198 conditions; 24 paired multipath timing conditions checked; all 18 paired M40/M64 full BER counts matched.
- Retained full AWGN evidence: all 18 points per accumulation setting use 33,554,432 bits; exact errors/bits accounting independently checked. These retained full experiments were not rerun during publication preparation.
- MATLAB regenerated 10 main and 4 multipath figures from C++ measurements. `verify_ofdm_figures` reopened all 14 saved figures and checked 140 measured series against CSVs; the PNG exports were visually reviewed.
- Core C++ implementation and existing C++ tests are byte-identical to the source implementation used to prepare this public snapshot.

GitHub Actions independently builds and tests C++ on Linux and macOS, runs the reduced-bit simulation, and checks retained data. It does not rerun the full BER/multipath studies or MATLAB verification.
