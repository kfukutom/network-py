import pytest
from your_urban_mapping_module import UrbanNetwork # REQUIRED <--- Change this to the name of your module

def test_urban_network_creation():
    network = UrbanNetwork()
    assert network is not None

def test_add_node():
    network = UrbanNetwork()
    network.add_node(1, 2, 1)
    assert network.node_count == 1

def test_add_edge():
    network = UrbanNetwork()
    network.add_node(1, 2, 1)
    network.add_node(2, 3, 2)
    network.add_edge(1, 2, 1)
    assert network.edge_count == 1