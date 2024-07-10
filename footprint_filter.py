"""Determines wether a set of points falls inside a polygon defined by another
set of points. In the context of this research, determines if a celestial
object falls in the void catalog footprint"""

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import pandas as pd


voids = pd.read_excel('exported_dataFrames/voids.xlsx')

cel_obj = pd.read_excel()