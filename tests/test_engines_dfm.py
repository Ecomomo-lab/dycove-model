
def test_dfm_engine_imports():
    from dycove import DFM_hydro
    assert hasattr(DFM_hydro, "DFM")
    assert hasattr(DFM_hydro, "DFMEngine")