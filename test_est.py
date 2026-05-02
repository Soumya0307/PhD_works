import cirq
from pyLIQTR.utils.resource_analysis import estimate_resources

q = cirq.LineQubit(0)
c1 = cirq.Circuit(cirq.rz(0.1)(q))
c2 = cirq.Circuit(cirq.rz(0.01)(q))

print("c1 T-count:", estimate_resources(c1, profile=False).get('T', 0))
print("c2 T-count:", estimate_resources(c2, profile=False).get('T', 0))
