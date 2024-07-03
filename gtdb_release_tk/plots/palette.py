from dataclasses import dataclass

@dataclass
class Palette:
    colour1: str
    colour2: str
    colour3: str

# default palette: blue, orange, green
DEFAULT_PALETTE = Palette('#80b1d3', '#fdae6b', '#b3de69')

# color Blind palette from Tableau (https://tableaufriction.blogspot.com/2012/11/finally-you-can-use-tableau-data-colors.html): grey, orange, blue
COLOR_BLIND_PALETTE = Palette('#ABABAB', '#FFBC79', '#5F9ED1') 
