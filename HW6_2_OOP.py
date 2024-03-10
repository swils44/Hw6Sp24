# region imports
import numpy as np
import math
from scipy.optimize import fsolve
import random as rnd
# endregion

# region class definitions
class Fluid():
    #region constructor
    def __init__(self, mu=0.00089, rho=1000):
        '''
        default properties are for water
        :param mu: dynamic viscosity in Pa*s -> (kg*m/s^2)*(s/m^2) -> kg/(m*s)
        :param rho: density in kg/m^3
        '''
        self.mu = mu  # dynamic viscosity
        self.rho = rho  # density
        self.nu = self.mu / self.rho  # kinematic viscosity in units of m^2/s
    #endregion

class Node():
    #region constructor
    def __init__(self, Name='a', Pipes=[], ExtFlow=0):
        '''
        A node in a pipe network.
        :param Name: name of the node
        :param Pipes: a list of pipes connected to this node
        :param ExtFlow: any external flow into (+) or out (-) of this node in L/s
        '''
        self.name = Name
        self.pipes = Pipes
        self.extFlow = ExtFlow
    #endregion

    #region methods/functions
    def getNetFlowRate(self):
        '''
        Calculates the net flow rate into this node in L/s
        :return:
        '''
        Qtot = self.extFlow  # Start with the external flow
        for p in self.pipes:
            Qtot += p.getFlowIntoNode(self.name)  # Add the flow from each connected pipe
        return Qtot
    #endregion

class Loop():
    #region constructor
    def __init__(self, Name='A', Pipes=[]):
        '''
        Defines a loop in a pipe network.
        :param Name: name of the loop
        :param Pipes: a list/array of pipes in this loop
        '''
        self.name = Name
        self.pipes = Pipes
    #endregion

    #region methods/functions
    def getLoopHeadLoss(self):
        '''
        Calculates the net head loss as I traverse around the loop, in m of fluid.
        :return:
        '''
        deltaP = 0  # Initialize to zero
        startNode = self.pipes[0].startNode  # Begin at the start node of the first pipe
        for p in self.pipes:
            deltaP += p.getFlowHeadLoss(startNode)  # Calculates the head loss in the pipe
            # Move to the next node
            startNode = p.endNode if startNode != p.endNode else p.startNode
        return deltaP
    #endregion

class Pipe():
    #region constructor
    def __init__(self, Start='A', End='B', L=100, D=200, r=0.00025, fluid=Fluid()):
        '''
        Defines a generic pipe with orientation from lowest letter to highest, alphabetically.
        :param Start: the start node (string)
        :param End: the end node (string)
        :param L: the pipe length in m (float)
        :param D: the pipe diameter in mm (float)
        :param r: the pipe roughness in m  (float)
        :param fluid:  a Fluid object (typically water)
        '''
        self.startNode = min(Start, End)  # makes sure to use the lowest letter for startNode
        self.endNode = max(Start, End)  # makes sure to use the highest letter for endNode
        self.length = L
        self.r = r
        self.fluid = fluid  # the fluid in the pipe

        # other calculated properties
        self.d = D / 1000.0  # diameter in m
        self.relrough = r / self.d  # relative roughness for use later
        self.A = math.pi * (self.d / 2) ** 2  # pipe cross sectional area
        self.Q = 10  # working in units of L/s, just an initial guess
        self.vel = self.V()  # calculate the initial velocity of the fluid
        self.reynolds = self.Re()  # calculate the initial Reynolds number
    #endregion

    #region methods/functions
    def V(self):
        '''
        Calculate average velocity in the pipe for volumetric flow self.Q
        :return: the average velocity in m/s
        '''
        self.vel = (self.Q / 1000) / self.A  # Convert Q from L/s to m^3/s and divide by area
        return self.vel

    def Re(self):
        '''
        Calculate the Reynolds number under current conditions.
        :return: the Reynolds number
        '''
        self.reynolds = (self.fluid.rho * self.V() * self.d) / self.fluid.mu
        return self.reynolds

    def FrictionFactor(self):
        '''
        This function calculates the friction factor for a pipe based on the
        notion of laminar, turbulent and transitional flow.
        :return: the (Darcy) friction factor
        '''
        Re = self.Re()
        if Re < 2000:
            return 64 / Re
        elif Re > 4000:
            # Using Haaland's formula as an approximation for the Colebrook equation
            return (1. / (-1.8 * np.log10((6.9 / Re) + ((self.relrough / 3.7) ** 1.11)))) ** 2
        else:
            # Transitional flow: using a weighted average of laminar and turbulent flow friction factors
            lam = 64 / Re
            turb = (1. / (-1.8 * np.log10((6.9 / Re) + ((self.relrough / 3.7) ** 1.11)))) ** 2
            return lam + ((Re - 2000) / (4000 - 2000)) * (turb - lam)

    def frictionHeadLoss(self):
        '''
        Use the Darcy-Weisbach equation to find the head loss through a section of pipe.
        :return: the head loss in m of water
        '''
        g = 9.81  # m/s^2, acceleration due to gravity
        f = self.FrictionFactor()  # Darcy friction factor
        v = self.V()  # average velocity
        hl = (f * self.length * v**2) / (2 * g * self.d)  # Darcy-Weisbach head loss equation
        return hl

    def getFlowHeadLoss(self, s):
        '''
        Calculate the head loss for the pipe.
        :param s: the node I'm starting with in a traversal of the pipe
        :return: the signed headloss through the pipe in m of fluid
        '''
        nTraverse = 1 if s == self.startNode else -1
        nFlow = 1 if self.Q >= 0 else -1
        return nTraverse * nFlow * self.frictionHeadLoss()

    def Name(self):
        '''
        Gets the pipe name.
        :return: the pipe name
        '''
        return self.startNode + '-' + self.endNode

    def oContainsNode(self, node):
        '''
        Does the pipe connect to the node?
        :param node: the node to check
        :return: True if the node is one of the endpoints of the pipe
        '''
        return self.startNode == node or self.endNode == node

    def printPipeFlowRate(self):
        '''
        Prints the flow rate in the pipe.
        '''
        print('The flow in segment {} is {:0.2f} L/s'.format(self.Name(), self.Q))

    def getFlowIntoNode(self, n):
        '''
        Determines the flow rate into node n.
        :param n: the node name
        :return: the flow rate into the node (signed)
        '''
        if n == self.startNode:
            return -self.Q
        elif n == self.endNode:
            return self.Q
        else:
            raise ValueError('Node is not connected to this pipe.')
    #endregion

# ... (other class definitions) ...

class PipeNetwork():
    #region constructor
    def __init__(self, Pipes=[], Loops=[], Nodes=[], fluid=Fluid()):
        '''
        The pipe network is built from pipe, node, loop, and fluid objects.
        :param Pipes: a list of pipe objects
        :param Loops: a list of loop objects
        :param Nodes: a list of node objects
        :param fluid: a fluid object
        '''
        self.loops = Loops
        self.nodes = Nodes
        self.Fluid = fluid
        self.pipes = Pipes
    #endregion

    #region methods/functions
    def findFlowRates(self):
        '''
        a method to analyze the pipe network and find the flow rates in each pipe
        given the constraints of: i) no net flow into a node and ii) no net pressure drops in the loops.
        :return: a list of flow rates in the pipes
        '''
        N = len(self.nodes) + len(self.loops)  # the number of equations
        Q0 = np.full(N, 10.0)  # initial guess for the flow rates

        def fn(Q):
            '''
            Callback function for fsolve. The mass continuity equations at the nodes and the loop equations
            are functions of the flow rates in the pipes.
            :param Q: an array of flow rates in the pipes
            :return: an array containing flow rates at the nodes and pressure losses for the loops
            '''
            for i, p in enumerate(self.pipes):
                p.Q = Q[i]  # Update the flow rate for each pipe

            qNet = self.getNodeFlowRates()  # Net flow rates at nodes
            lhl = self.getLoopHeadLosses()  # Loop head losses

            return qNet + lhl

        # Using fsolve to find the flow rates
        FR = fsolve(fn, Q0)
        return FR

    def getNodeFlowRates(self):
        '''
        Calculate the net flow rate for each node object.
        :return: a list of net flow rates for the nodes
        '''
        qNet = [n.getNetFlowRate() for n in self.nodes]
        return qNet

    def getLoopHeadLosses(self):
        '''
        Calculate the net head loss for each loop object.
        :return: a list of net head losses for the loops
        '''
        lhl = [l.getLoopHeadLoss() for l in self.loops]
        return lhl

    def getPipe(self, name):
        '''
        Returns a pipe object by its name.
        :param name: the name of the pipe
        :return: the pipe object
        '''
        for p in self.pipes:
            if name == p.Name():
                return p
        return None

    def getNodePipes(self, node):
        '''
        Returns a list of pipe objects that are connected to the node object.
        :param node: the node name
        :return: a list of connected pipe objects
        '''
        return [p for p in self.pipes if p.oContainsNode(node)]

    def nodeBuilt(self, node):
        '''
        Determines if a node object has already been constructed.
        :param node: the node name
        :return: True if the node exists
        '''
        return any(n.name == node for n in self.nodes)

    def getNode(self, name):
        '''
        Returns a node object by name.
        :param name: the name of the node
        :return: the node object
        '''
        for n in self.nodes:
            if n.name == name:
                return n
        return None

    def buildNodes(self):
        '''
        Automatically create the node objects by looking at the pipe ends.
        '''
        for p in self.pipes:
            if not self.nodeBuilt(p.startNode):
                self.nodes.append(Node(p.startNode, self.getNodePipes(p.startNode)))
            if not self.nodeBuilt(p.endNode):
                self.nodes.append(Node(p.endNode, self.getNodePipes(p.endNode)))

    def printPipeFlowRates(self):
        '''
        Prints the flow rate in each pipe.
        '''
        for p in self.pipes:
            p.printPipeFlowRate()

    def printNetNodeFlows(self):
        '''
        Prints the net flow into each node.
        '''
        for n in self.nodes:
            print('net flow into node {} is {:0.2f} L/s'.format(n.name, n.getNetFlowRate()))

    def printLoopHeadLoss(self):
        '''
        Prints the head loss for each loop.
        '''
        for l in self.loops:
            print('head loss for loop {} is {:0.2f}'.format(l.name, l.getLoopHeadLoss()))
    #endregion
# endregion

# region function definitions
# region function definitions
def main():
    '''
    This program analyzes flows in a given pipe network based on the following:
    1. The pipe segments are named by their endpoint node names:  e.g., a-b, b-e, etc.
    2. Flow from the lower letter to the higher letter of a pipe is considered positive.
    3. Pressure decreases in the direction of flow through a pipe.
    4. At each node in the pipe network, mass is conserved.
    5. For any loop in the pipe network, the pressure loss is zero
    Approach to analyzing the pipe network:
    Step 1: build a pipe network object that contains pipe, node, loop and fluid objects
    Step 2: calculate the flow rates in each pipe using fsolve
    Step 3: output results
    Step 4: check results against expected properties of zero head loss around a loop and mass conservation at nodes.
    :return:
    '''
    # instantiate a Fluid object to define the working fluid as water
    water = Fluid(mu=0.00089, rho=1000)  # water properties
    roughness = 0.00025  # in meters

    # instantiate a new PipeNetwork object
    PN = PipeNetwork(fluid=water)  # instantiate PipeNetwork with the water fluid
    # add Pipe objects to the pipe network (see constructor for Pipe class)
    PN.pipes.append(Pipe('a', 'b', 250, 300, roughness, water))
    PN.pipes.append(Pipe('a', 'c', 100, 200, roughness, water))
    PN.pipes.append(Pipe('b', 'e', 100, 200, roughness, water))
    PN.pipes.append(Pipe('c', 'd', 125, 200, roughness, water))
    PN.pipes.append(Pipe('c', 'f', 100, 150, roughness, water))
    PN.pipes.append(Pipe('d', 'e', 125, 200, roughness, water))
    PN.pipes.append(Pipe('d', 'g', 100, 150, roughness, water))
    PN.pipes.append(Pipe('e', 'h', 100, 150, roughness, water))
    PN.pipes.append(Pipe('f', 'g', 125, 250, roughness, water))
    PN.pipes.append(Pipe('g', 'h', 125, 250, roughness, water))
    # add Node objects to the pipe network by calling buildNodes method of PN object
    PN.buildNodes()

    # update the external flow of certain nodes
    PN.getNode('a').extFlow = 60
    PN.getNode('d').extFlow = -30
    PN.getNode('f').extFlow = -15
    PN.getNode('h').extFlow = -15

    # add Loop objects to the pipe network
    PN.loops.append(Loop('A', [PN.getPipe('a-b'), PN.getPipe('b-e'), PN.getPipe('d-e'), PN.getPipe('c-d'), PN.getPipe('a-c')]))
    PN.loops.append(Loop('B', [PN.getPipe('c-d'), PN.getPipe('d-g'), PN.getPipe('f-g'), PN.getPipe('c-f')]))
    PN.loops.append(Loop('C', [PN.getPipe('d-e'), PN.getPipe('e-h'), PN.getPipe('g-h'), PN.getPipe('d-g')]))

    # call the findFlowRates method of the PN (a PipeNetwork object)
    PN.findFlowRates()

    # get output
    PN.printPipeFlowRates()
    print()
    print('Check node flows:')
    PN.printNetNodeFlows()
    print()
    print('Check loop head loss:')
    PN.printLoopHeadLoss()
    #PN.printPipeHeadLosses()
# endregion

# region function calls
if __name__ == "__main__":
    main()
# endregion
