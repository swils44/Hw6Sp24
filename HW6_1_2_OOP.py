from HW6_1_OOP import ResistorNetwork, Resistor, VoltageSource, Loop
from scipy.optimize import fsolve

# region class definitions
class ResistorNetwork2(ResistorNetwork):
    def __init__(self):
        super().__init__()

    def GetElementDeltaV(self, elementName):
        element = self.GetElementByName(elementName)
        if isinstance(element, Resistor):
            return element.GetVoltageDrop()
        elif isinstance(element, VoltageSource):
            return element.Value
        else:
            return 0  # Return 0 for unknown elements

    # Rest of the ResistorNetwork2 class code...

    def AnalyzeCircuit(self):
        """
        Use fsolve to find currents in the resistor network.
        1. KCL:  The total current flowing into any node in the network is zero.
        2. KVL:  When traversing a closed loop in the circuit, the net voltage drop must be zero.
        :return: a list of the currents in the resistor network
        """
        # need to set the currents to that Kirchoff's laws are satisfied
        i0 = [0.1, 0.1, 0.1]  # define an initial guess for the currents in the circuit
        i = fsolve(self.GetKirchoffVals, i0)
        # print output to the screen
        print("I1 = {:0.1f}".format(i[0]))
        print("I2 = {:0.1f}".format(i[1]))
        print("I3 = {:0.1f}".format(i[2]))
        return i

    def GetKirchoffVals(self, i):
        """
        This function uses Kirchoff Voltage and Current laws to analyze this specific circuit
        KVL:  The net voltage drop for a closed loop in a circuit should be zero
        KCL:  The net current flow into a node in a circuit should be zero
        :param i: a list of currents relevant to the circuit
        :return: a list of loop voltage drops and node currents
        """
        # set current in resistors in the top loop.
        self.GetResistorByName('ad').Current = i[0]  # I_1 in diagram
        self.GetResistorByName('bc').Current = i[0]  # I_1 in diagram
        self.GetResistorByName('cd').Current = i[2]  # I_3 in diagram
        # set current in resistor in bottom loop.
        self.GetResistorByName('ce').Current = i[1]  # I_2 in diagram
        # calculate net current into node c
        Node_c_Current = sum([i[0], i[1], -i[2]])

        KVL = self.GetLoopVoltageDrops()  # two equations here
        KVL.append(Node_c_Current)  # one equation here
        return KVL

def main():
    net = ResistorNetwork2()
    net.BuildNetworkFromFile('ResistorNetwork_2.txt')
    iVals = net.AnalyzeCircuit()
    print(f"Calculated currents: {iVals}")

if __name__ == "__main__":
    main()