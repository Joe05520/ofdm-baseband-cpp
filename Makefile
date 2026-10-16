CXX = c++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -Wpedantic
SOURCES = src/fft.cpp src/ofdm.cpp
TEST_DIR := $(if $(wildcard tests/test_ofdm.cpp),tests,src/tests)

.PHONY: all test run sanitize clean test-multipath run-multipath sanitize-multipath
all: bin/ofdm bin/sync_study bin/multipath_study bin/timing_compare
bin:
	mkdir -p bin
bin/ofdm: $(SOURCES) src/main.cpp src/fft.hpp src/ofdm.hpp src/sync_config.hpp | bin
	$(CXX) $(CXXFLAGS) $(SOURCES) src/main.cpp -o $@
bin/test_ofdm: $(SOURCES) $(TEST_DIR)/test_ofdm.cpp src/fft.hpp src/ofdm.hpp | bin
	$(CXX) $(CXXFLAGS) -Isrc $(SOURCES) $(TEST_DIR)/test_ofdm.cpp -o $@
bin/sync_study: $(SOURCES) src/sync_study.cpp src/fft.hpp src/ofdm.hpp | bin
	$(CXX) $(CXXFLAGS) -Isrc $(SOURCES) src/sync_study.cpp -o $@
test: bin/test_ofdm
	./bin/test_ofdm
bin/multipath_study: $(SOURCES) src/multipath.cpp src/multipath_study.cpp src/multipath.hpp | bin
	$(CXX) $(CXXFLAGS) -Isrc $(SOURCES) src/multipath.cpp src/multipath_study.cpp -o $@
bin/test_multipath: $(SOURCES) src/multipath.cpp $(TEST_DIR)/test_multipath.cpp src/multipath.hpp | bin
	$(CXX) $(CXXFLAGS) -Isrc $(SOURCES) src/multipath.cpp $(TEST_DIR)/test_multipath.cpp -o $@
test-multipath: bin/test_multipath
	./bin/test_multipath
run-multipath: bin/multipath_study
	./bin/multipath_study results/multipath
sanitize-multipath: | bin
	$(CXX) -O1 -g -std=c++17 -fsanitize=address,undefined -fno-omit-frame-pointer -Isrc $(SOURCES) src/multipath.cpp $(TEST_DIR)/test_multipath.cpp -o bin/test_multipath_sanitize
	./bin/test_multipath_sanitize
run: bin/ofdm
	./bin/ofdm results
sanitize: | bin
	$(CXX) -O1 -g -std=c++17 -fsanitize=address,undefined -fno-omit-frame-pointer -Isrc $(SOURCES) $(TEST_DIR)/test_ofdm.cpp -o bin/test_ofdm_sanitize
	./bin/test_ofdm_sanitize
clean:
	rm -f bin/ofdm bin/sync_study bin/multipath_study bin/timing_compare bin/test_ofdm bin/test_ofdm_sanitize bin/test_multipath bin/test_multipath_sanitize

bin/timing_compare: $(SOURCES) src/multipath.cpp src/timing_compare.cpp src/multipath.hpp | bin
	$(CXX) $(CXXFLAGS) -Isrc $(SOURCES) src/multipath.cpp src/timing_compare.cpp -o $@
