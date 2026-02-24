from pytest import mark, approx

@mark.dfm
@mark.unit
def test_dfm_engine_imports():
    from dycove import DFM_hydro
    assert callable(DFM_hydro.DFM)
    assert callable(DFM_hydro.DFMEngine)

@mark.dfm
@mark.unit
def test_engine_abstractmethods(required_engine_methods):
    from dycove import DFM_hydro
    for method in required_engine_methods:
        assert hasattr(DFM_hydro.DFMEngine, method)