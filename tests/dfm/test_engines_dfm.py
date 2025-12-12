from pytest import mark, approx

@mark.dfm
def test_dfm_engine_imports():
    from dycove import DFM_hydro
    assert hasattr(DFM_hydro, "DFM")
    assert hasattr(DFM_hydro, "DFMEngine")

@mark.dfm
def test_engine_abstractmethods(required_engine_methods):
    from dycove import DFM_hydro
    for method in required_engine_methods:
        assert hasattr(DFM_hydro.DFMEngine, method)