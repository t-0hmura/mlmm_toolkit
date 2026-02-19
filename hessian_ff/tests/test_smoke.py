def test_public_api_symbols() -> None:
    import hessian_ff

    required = [
        "clear_runtime_cache",
        "load_system",
        "load_coords",
        "ForceFieldTorch",
        "system_summary",
        "torch_energy",
        "torch_energy_batch",
        "torch_force",
        "torch_force_batch",
        "torch_hessian",
        "verify_openmm",
    ]
    for name in required:
        assert hasattr(hessian_ff, name)
