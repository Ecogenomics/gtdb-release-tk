from dataclasses import dataclass

@dataclass

class Palette:
    grey: str
    orange: str
    blue: str


COLOR_BLIND_PALETTE = Palette('#ABABAB', '#FFBC79', '#5F9ED1')