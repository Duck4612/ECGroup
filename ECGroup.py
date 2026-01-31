
import math

class Point():
    def __init__(self, x, y):
        #O is the point (math.inf, math.inf)
        self.x = x
        self.y = y

    def __eq__(self, other):
        if not(isinstance(other, Point)):
            return False
        return (self.x == other.x) and (self.y == other.y)

    def __str__(self):
        return f"(x,y): ({self.x},{self.y})"

class EllipticCurveGroup():
    def __init__(self,a: int,b: int,p: int):
        self.a = a
        self.b = b
        self.p = p
        self.validate()
        self.points = self.enumerate_points()

    def validate(self):
        if not self.is_prime(self.p):
            raise ValueError(f"p ({self.p}) is not prime")
        if (4 * self.a**3 + 27 * self.b**2) % self.p == 0:
            raise ValueError("Invalid elliptic curve parameters")
        
    @staticmethod
    def is_prime(n: int) -> bool:
        if n < 2:
            return False
        for i in range(2, int(n**0.5) + 1):
            if n % i == 0:
                return False
        return True

    def add(self, p1: Point, p2: Point) -> Point:
        if p1.x != p2.x:
            m = ((p2.y - p1.y)/(p2.x - p1.x)) % self.p
            x = (m**2 - p1.x - p2.x) % self.p
            y = (m*(p1.x-x)-p1.y) % self.p
            p3 = Point(x,y)
        elif p1.y != p2.y:
            p3 = Point(math.inf, math.inf)
        elif p1.y == 0:
            p3 = Point(math.inf, math.inf)
        else:
            m = ((3*(p1.x**2)+self.a)/(2*p1.y)) % self.p
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

    def validate_point(self, p1: Point) -> bool:
        return (p1.y**2 % self.p) == ((p1.x**3 + self.a*p1.x + self.b) % self.p)

def main():
    a = 3
    b = 3
    p = 7
    test_group = EllipticCurveGroup(a,b,p)
    for point in test_group.points:
        print(point)

if __name__ == "__main__":
    main()