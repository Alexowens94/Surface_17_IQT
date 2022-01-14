from error import *
from logical_qubit import *
from projectq import MainEngine
from projectq.backends import Simulator
from projectq.cengines import DummyEngine
import pytest

@pytest.fixture
def main_engine():
    return MainEngine(backend=Simulator(), engine_list=[DummyEngine()])

def test_noisy_gate():
    # Test that when you add a noisy gate it is applying
    # the gate and your command list length goes up by 2
    dummy_engine = DummyEngine(save_commands=True)
    eng = MainEngine(backend=Simulator(), engine_list=[dummy_engine])
    data = eng.allocate_qureg(1)
    ancilla = eng.allocate_qureg(1)
    rxx_gate = NoisyGate(Rxx(pi / 2))
    rxx_gate.apply((data[0], ancilla[0]))
    rx_gate = NoisyGate(Rx(pi / 2))
    rx_gate.apply((data[0]))
    All(Measure) | data
    All(Measure) | ancilla
    assert len(dummy_engine.received_commands) == 6

test_noisy_gate()