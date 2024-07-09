"""To determine a 2D footprint of the cosmic voids from the Sutter Catalog."""

import pandas as pd
import alphashape # used to determine the outer perimeter of a set of points


voids = pd.read_excel('exported_dataFrames/voids.xlsx')
