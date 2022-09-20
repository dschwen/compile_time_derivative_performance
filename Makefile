all: test-clang test-gcc test-nvc

COMMON = -I$(MOOSE_DIR)/framework/include/utils -O2 --std=c++17 $(shell $(LIBMESH_DIR)/bin/libmesh-config --include)
clean:
	rm test-*

test-clang: test.C
	clang++ -o test-clang test.C $(COMMON)

test-gcc: test.C
	g++ -o test-gcc test.C $(COMMON)

test-nvc: test.C
	nvc++ -o test-nvc test.C $(COMMON)
