from HW6_1_OOP import ResistorNetwork, Resistor, VoltageSource, Loop
from scipy.optimize import fsolve

# region class definitions
class ResistorNetwork2(ResistorNetwork):
    def __init__(self):
        super().__init__()

        self.Loops = []  # initialize an empty list of loop objects in the network
        self.Resistors = []  # initialize an empty a list of resistor objects in the network
        self.VSources = []  # initialize an empty a list of source objects in the network


    def AnalyzeCircuit(self):
        """
        Use fsolve to find currents in the resistor network.
        1. KCL:  The total current flowing into any node in the network is zero.
        2. KVL:  When traversing a closed loop in the circuit, the net voltage drop must be zero.
        :return: a list of the currents in the resistor network
        """
        # need to set the currents to that Kirchoff's laws are satisfied
        i0 = [1, 1, 1, 1, 1]  # define an initial guess for the currents in the circuit
        i = fsolve(self.GetKirchoffVals, i0)
        # print output to the screen
        print("I1 = {:0.1f}A".format(i[0]))
        print("I2 = {:0.1f}A".format(i[1]))
        print("I3 = {:0.1f}A".format(i[2]))
        print("I4 = {:0.1f}A".format(i[3]))
        print("I5 = {:0.1f}A".format(i[4]))
        return i

    def GetKirchoffVals(self, i):
        """
        This function uses Kirchoff Voltage and Current laws to analyze this specific circuit
        KVL:  The net voltage drop for a closed loop in a circuit should be zero
        KCL:  The net current flow into a node in a circuit should be zero
        :param i: a list of currents relevant to the circuit
        :return: a list of loop voltage drops and node currents
        """
        # set current in resistors in the top loop
        self.GetResistorByName('ad').Current = i[0]  # I_1 in diagram
        self.GetResistorByName('bc').Current = i[0]  # I_2 in diagram

        # set current in resistors in the bottom loop
        self.GetResistorByName('ce').Current = i[4]  # Current 5
        self.GetResistorByName('de').Current = i[3]  # Current 4

        # calculate current through resistor 'cd'
        self.GetResistorByName('cd').Current = i[2] # Current 3

        # calculate net current into node c
        Node_c_Current = sum([i[4], -i[2], i[0]])

        # calculate net current into node d
        Node_d_Current = sum([-i[1], i[2], -i[0], i[3]])


        KVL = self.GetLoopVoltageDrops()  # three equations here
        KVL.append(Node_c_Current)  # one equation here
        KVL.append(Node_d_Current)  # one equation here

        return KVL

def main():
    net = ResistorNetwork2()
    net.BuildNetworkFromFile('ResistorNetwork_2.txt')
    iVals = net.AnalyzeCircuit()
    print(f"Calculated currents: {[f'{i:.1f}A' for i in iVals]}")

if __name__ == "__main__":#
    main()