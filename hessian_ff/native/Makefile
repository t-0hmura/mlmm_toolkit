.PHONY: all clean

all:
	python -c "from hessian_ff.native import build_native_extensions; \
	    result = build_native_extensions(verbose=True, force_rebuild=True); \
	    assert all(value == 'ok' for value in result.values()), result"

clean:
	rm -rf .build_nonbonded .build_bonded .build_analytical_hessian
