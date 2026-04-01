"""Tests for my_pack modules."""


def test_counter_increments():
    from my_pack import Counter

    counter = Counter(name="test", min_dt=0.1)
    counter.reset()

    counter.advance_to(0.1)
    counter.advance_to(0.2)

    outputs = counter.get_outputs()
    assert outputs["count"].value["count"] == 2


def test_counter_reset():
    from my_pack import Counter

    counter = Counter(min_dt=0.1)
    counter.advance_to(0.1)
    counter.reset()
    counter.advance_to(0.2)

    outputs = counter.get_outputs()
    assert outputs["count"].value["count"] == 1


def test_counter_visualize():
    from my_pack import Counter

    counter = Counter(name="viz_test", min_dt=0.1)
    assert counter.visualize() is None

    counter.advance_to(0.1)
    vis = counter.visualize()
    assert vis is not None
    assert vis["render"] == "timeseries"
    assert "viz_test" in vis["data"]["series"][0]["name"]


def test_accumulator_accumulates():
    from my_pack import Accumulator
    from biosim import BioSignal

    acc = Accumulator(initial=10.0, min_dt=0.1)
    acc.set_inputs({"value": BioSignal(source="src", name="value", value=5.0, time=0.1)})
    acc.advance_to(0.1)
    acc.set_inputs({"value": BioSignal(source="src", name="value", value=3.0, time=0.2)})
    acc.advance_to(0.2)

    outputs = acc.get_outputs()
    assert outputs["total"].value["total"] == 18.0
