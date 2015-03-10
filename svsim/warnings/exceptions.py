class OverlappingFeaturesError(Exception):
    """Error for ovelapping features"""
    def __init__(self, contig, position):
        super(OverlappingFeaturesError, self).__init__()
        self.contig = contig
        self.position = position