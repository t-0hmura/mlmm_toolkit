.PHONY: all clean

all:
	python -c "from hessian_ff.native import build_native_extensions; \
	    build_native_extensions(verbose=True, force_rebuild=True)"

clean:
	rm -rf .build_nonbonded .build_bonded .build_analytical_hessian
