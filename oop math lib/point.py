# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from vector import Vector


# A basic wrapper around the vector class to be able to classify the point.
# This can be made part of the vector class but it more organized this way.
# A point though does not have any size (eg. magnitude) only a location. :P
class Point(object):
    '''Class for a point'''
    # Takes in a vector to define the point.
    # The location is secondary type input (not necessary).
    # The location is used to determine where a point is located relative to a plane.
    def __init__(self, *args):
        self.point = Vector([0.0, 0.0, 0.0])
        self.location = "NONE"
        if isinstance(args[0], Vector):
            # This is a vector. :/ So if you wanna access the contents of this you have to do:
            # your_point.point.vec[0...n] or whatever.
            self.point = args[0]
            # The location can be either FRONT, BEHIND, ON_PLANE.
            if len(args) > 1:
                if args[1] == "FRONT" or args[1] == "BEHIND" or args[1] == "ON_PLANE":
                    self.location = args[1]
                else:
                    print "Oh shi- Location,", self.location, "undefined. Only FRONT, BEHIND and ON_PLANE, allowed. :["
        else:
            print "LOL bro, it's not a vector. Feed me vectors tho. Otherwise there is no point to it. LOL 'point' get it? ROFL"

    def duplicate(self):
        return Point(self.point, self.location)

    def output(self):
        print "Point:"
        print "Values:"
        self.point.output()
        print "Location: ", self.location
