all: build/CMakeCache.txt run_make
.PHONY: all

# Run cmake if not yet done or if CMakeLists.txt has changed.
build/CMakeCache.txt: CMakeLists.txt
	@echo "Running cmake"
	@mkdir -p build
	@cd build && cmake ..

run_make: build/CMakeCache.txt
	@echo "Running make"
	$(MAKE) -C build
.PHONY: run_make

update:
	@touch src/CMakeLists.txt
	@touch test/src/CMakeLists.txt
	$(MAKE) -C build
.PHONY: update

unittest: update
	@./test/bin/epa_test
.PHONY: test

clean:
	@echo "Cleaning"
	@rm -rf build
	@rm -rf bin
	@rm -rf test/bin
.PHONY: clean

#======================================
#===		Test commands follow				===
#======================================
EPABIN=./bin/epa-ng
TEST=test/data/lucas
TREE=$(TEST)/tree.newick
REF=$(TEST)/reference.fasta
QRY=$(TEST)/query.fasta.bfast
INFO=$(TEST)/infofile
BINFILE=$(TEST)/epa_binary_file
OUTDIR=/tmp/epa

BINARY_WRITE= -t $(TREE) -s $(REF) -B -w $(OUTDIR) --verbose $(F)
BINARY_READ=-b $(BINFILE) -q $(QRY) -w $(OUTDIR) -g 0.99 --verbose $(F)
NORM_TEST=-t $(TREE) -s $(REF) -q $(QRY) --model $(INFO)  -w $(OUTDIR)  --verbose $(F)

test: #update
	mkdir -p $(OUTDIR)
	rm -f $(OUTDIR)/*
	$(EPABIN) $(NORM_TEST) --threads 4
.PHONY: test

bintest: update
	mkdir -p $(OUTDIR)
	rm -f $(OUTDIR)/*
	$(EPABIN) $(BINARY_WRITE)
	$(EPABIN) $(BINARY_READ)
.PHONY: bintest

mpi_test: #update
	mkdir -p $(OUTDIR)
	rm -f $(OUTDIR)/*
	mpirun -n 2 $(EPABIN) $(NORM_TEST) --threads 2
.PHONY: mpi_test
