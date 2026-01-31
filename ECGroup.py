
import math

class Point():
    def __init__(self, x: int, y: int, order = None) -> None:
        #O is the point (math.inf, math.inf)
        self.x = x
        self.y = y
        self.order = order

    def __eq__(self, other) -> bool:
        if not(isinstance(other, Point)):
            return False
        return (self.x == other.x) and (self.y == other.y)

    def __hash__(self):
        return hash((self.x, self.y))

    def __str__(self):
        return f"(x,y): ({self.x},{self.y}) | {self.order}"

class ECGroup():
    def __init__(self, a: int, b: int, p: int):
        self.a = a
        self.b = b
        self.p = p
        self.validate()
        self.points = self.enumerate_points()
        self.points_set = set(self.points)

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

    def __str__(self):
        string = ""
        for point in self.points:
            string += str(point) +"\n"
        string = string[:-1]
        return string


    

def main():
    p = 107

    if False:
        for a in range(0,p):
            for b in range(0,p):
                if ECGroup.is_prime(p) and ((4 * a**3 + 27 * b**2) % p != 0):
                    possible_group = ECGroup(a,b,p)
                    if ECGroup.is_prime(len(possible_group.points)) and ECGroup.is_prime(a) and ECGroup.is_prime(b):
                        print(f"a: {a} | b: {b} | p: {p} | Order: {len(possible_group.points)}")

    a = 11
    b = 13

    ecgroup = ECGroup(11,13,p)

    gen = ecgroup.find_generators()
    print()
    for gener in gen:
        print(gener)
    print(len(gen))
    print(len(ecgroup.points))

    ecgroup.isomorphism(gen[30])
    print(ecgroup)

    

if __name__ == "__main__":
    main()