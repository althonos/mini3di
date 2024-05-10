"""A NumPy port of the ``foldseek`` code for encoding structures to 3di.
"""

__all__ = ["Encoder", "FeatureEncoder", "PartnerIndexEncoder", "VirtualCenterEncoder"]
__version__ = "0.1.1"
__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPL-3.0-only"
__credits__ = "Martin Steinegger and his lab for ``foldseek``."

from .encoder import Encoder, FeatureEncoder, PartnerIndexEncoder, VirtualCenterEncoder
