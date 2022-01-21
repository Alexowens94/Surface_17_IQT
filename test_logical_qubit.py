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
    rxx_gate = NoisyGate(Rxx(pi / 2), (data[0], ancilla[0]))
    rxx_gate.apply()
    rx_gate = NoisyGate(Rx(pi / 2), (data[0]))
    rx_gate.apply()
    All(Measure) | data
    All(Measure) | ancilla
    assert len(dummy_engine.received_commands) == 6

def test_stabilizer_cycle():
    backend = Simulator()
    dummy_engine = DummyEngine(save_commands=True)
    eng = MainEngine(backend, [dummy_engine])
    data = eng.allocate_qureg(9)
    ancilla = eng.allocate_qureg(8)
    syndrome0 = StabiliserCycle(data=data, ancilla=ancilla, location=[0], error_model=pdd_error_model)
    assert len(syndrome0.all_gate_locations) == 60
    assert len(syndrome0.gate_locations[Rxx]) == 24
    assert len(syndrome0.gate_locations[Rx]) == 12
    assert len(syndrome0.gate_locations[Ry]) == 24
    assert len(syndrome0.error_locations[YCtrlError]) == 24
    assert len(syndrome0.error_locations[XCtrlError]) == 36
    assert len(syndrome0.error_locations[DephasingError]) == 60

test_noisy_gate()

test_stabilizer_cycle()