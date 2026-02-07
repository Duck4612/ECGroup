import math
class Point():
    def __init__(self, x: int, y: int, order = None) -> None:
        #O is the point (math.inf, math.inf)
        self.x = x
        self.y = y
        self.order = order
        self.orbit = []

    def __eq__(self, other) -> bool:
        if not(isinstance(other, Point)):
            return False
        return (self.x == other.x) and (self.y == other.y)

    def __hash__(self):
        return hash((self.x, self.y))

    def __str__(self):
        return f"({self.x},{self.y})"

class ECGroup():
    def __init__(self, a: int, b: int, p: int):
        self.a = a
        self.b = b
        self.p = p
        self.validate()
        self.points = self.enumerate_points()
        self.points_set = set(self.points)
        self.generators = self.find_generators()

    def validate(self):
        if not self.is_prime(self.p):
            raise ValueError(f"p ({self.p}) is not prime")
        if (4 * self.a**3 + 27 * self.b**2) % self.p == 0:
            raise ValueError(f"Invalid elliptic curve parameters: a: {self.a} | b: {self.b} | p: {self.p}")
        
    @staticmethod
    def is_prime(n: int) -> bool:
        if n < 2:
            return False
        for i in range(2, int(n**0.5) + 1):
            if n % i == 0:
                return False
        return True

    def multiply(self, p1: Point, k: int) -> Point:
        p2 = Point(math.inf, math.inf)
        for i in range(k):
            p2 = self.add(p1, p2)
        return p2

    def add(self, p1: Point, p2: Point) -> Point:
        if p1 == Point(math.inf,math.inf):
            p3 = p2
        elif p2 == Point(math.inf,math.inf):
            p3 = p1
        elif p1.x != p2.x:
            inv = pow(p2.x - p1.x, -1, self.p)
            m = ((p2.y - p1.y) * inv) % self.p
            x = (m**2 - p1.x - p2.x) % self.p
            y = (m*(p1.x-x)-p1.y) % self.p
            p3 = Point(x,y)
        elif p1.y != p2.y:
            p3 = Point(math.inf, math.inf)
        elif p1.y == 0:
            p3 = Point(math.inf, math.inf)
        else:
            inv = pow(2*p1.y, -1, self.p)
            m = ((3*(p1.x**2)+self.a)* inv) % self.p
            x = (m**2 - 2*p1.x) % self.p
            y = (m*(p1.x-x)-p1.y) % self.p
            p3 = Point(x,y)

        for point in self.points:
            if point == p3:
                return point

        return p3

    def f(self,x: int) -> int:
        return (x**3 + self.a*x + self.b) % self.p

    def enumerate_points(self) -> list:
        points = []
        for x in range(self.p):
            for y in range(self.p):
                if (y*y) % self.p == self.f(x):
                    points.append(Point(x, y))
        points.append(Point(math.inf, math.inf))  # point at infinity
        return points

    def find_generators(self) -> list:
        generators = []
        for start_point in self.points:
            generated = []
            point = start_point
            while not(point in generated):
                generated.append(point)
                point = self.add(point, start_point)
            if set(generated) == self.points_set:
                generators.append(start_point)
        return generators
        
    def isomorphism(self, start_point = None):
        if start_point is None:
            if len(self.generators) == 0:
                return
            start_point = self.generators[0]
        ordered_points = []
        ordered_points.append(Point(math.inf,math.inf, 0))
        point = start_point
        order = 1
        while not(point == Point(math.inf,math.inf)):
            point.order = order
            ordered_points.append(point)
            point = self.add(point,start_point)
            order += 1
        self.points = ordered_points

    def find_orbit(self, p1: Point):
        point = p1
        while not(point in p1.orbit):
            p1.orbit.append(point)
            point = self.add(point, p1)

    def __str__(self):
        string = ""
        for point in self.points:
            string += str(point) +"\n"
        string = string[:-1]
        return string