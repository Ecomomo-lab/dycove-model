
def test_anuga_engine_imports():
    from dycove import ANUGA_hydro
    assert hasattr(ANUGA_hydro, "ANUGA")
    assert hasattr(ANUGA_hydro, "AnugaEngine")