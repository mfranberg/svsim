class OverlappingFeaturesError(Exception):
    """Error for ovelapping features"""
    def __init__(self, contig, position):
        super(OverlappingFeaturesError, self).__init__()
        self.contig = contig
        self.position = position

class SimulatorNotFoundError(Exception):
    """Error for non existing simulator"""
    def __init__(self, simulator):
        super(SimulatorNotFoundError, self).__init__()
        self.name = simulator
