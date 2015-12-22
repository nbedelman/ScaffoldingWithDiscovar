        
class NoSegmentError(Exception):
     def __init__(self, value):
         self.value = value
     def __str__(self):
         return repr(self.value)

class NotInformativeError(Exception):
     def __init__(self, value):
         self.value = value
     def __str__(self):
         return repr(self.value)

class NestedSegsError(Exception):
     def __init__(self, value):
         self.value = value
     def __str__(self):
         return repr(self.value)
