class Polynomial:
    def __init__(self, *args):
        if isinstance(args[0], list):
            self.coeffs = args[0]
        elif isinstance(args[0], Polynomial):
            self.coeffs = args[0].coeffs[:]
        elif isinstance(args[0], dict):
            self.coeffs = [args[0][b] if b in args[0] else 0
                           for b in range(max(args[0].keys()) + 1)]
        else:
            self.coeffs = [*args]
        while len(self.coeffs) > 1 and not self.coeffs[-1]:
            self.coeffs.pop()

    def __repr__(self):
        return "Polynomial " + str(self.coeffs)

    def __str__(self):
        s = ''
        l = len(self.coeffs)
        for i in range(l):
            power = l - i - 1
            if self.coeffs[power] != 0:
                k = (abs(self.coeffs[power]) != 1 or power == 0)
                s += ' + ' * (self.coeffs[power] > 0 and power != l - 1)
                s += ' ' * (self.coeffs[power] < 0 and power != l - 1)
                s += '-' * (self.coeffs[power] < 0)
                s += ' ' * (self.coeffs[power] < 0 and power != l - 1)
                s += str(abs(self.coeffs[power])) * k
                s += 'x' * (power >= 1) + ('^' + str(power)) * (power >= 2)
        return s

    def __eq__(self, other):
        other = Polynomial(other)
        return self.coeffs == other.coeffs

    def __add__(self, other):
        other = Polynomial(other)
        m = min(self.degree(), other.degree())
        s = [self.coeffs[i] + other.coeffs[i]
                  for i in range(m + 1)]
        s += self.coeffs[m + 1:]
        s += other.coeffs[m + 1:]
        return Polynomial(s)

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return self * (-1)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -(self - other)

    def __call__(self, val):
        a = 0
        for i in range(len(self.coeffs)):
            a += self.coeffs[i] * (val ** i)
        return a

    def __mul__(self, other):
        if isinstance(other, Polynomial):
            a = [0 for i in range((len(self.coeffs)) *
                                    (len(other.coeffs) - 1) + 1)]
            for i in range(len(self.coeffs)):
                for j in range(len(other.coeffs)):
                    a[i + j] += other.coeffs[j] * self.coeffs[i]
        else:
            a = [x * other for x in self.coeffs]
        return Polynomial(a)

    def __rmul__(self, other):
        return self * other

    def degree(self):
        return len(self.coeffs) - 1

    def der(self, d=1):
        if d > (len(self.coeffs) - 1):
            return Polynomial(0)
        n = Polynomial(self)
        for k in range(d):
            for i in range(len(n.coeffs)):
                n.coeffs[i] *= i
            n.coeffs = n.coeffs[1:]
        if len(n.coeffs) == 0:
            return Polynomial([0])
        return n

    def __iter__(self):
        a = [(x, y) for x, y in enumerate(self.coeffs)]
        return iter(a)

    def __next__(self):
        return next(self)


class NotOddDegreeException(Exception):
    def __init__(self):
        pass


class RealPolynomial(Polynomial):
    def __init__(self, *args):
        super().__init__(*args)
        if self.degree() % 2 == 0:
            raise NotOddDegreeException

    def find_root(self):
        a = 1
        while True:
            if self(-a) * self(a) < 0:
                break
            else:
                a = a * 2
        a, b = -a, a
        while True:
            if abs(self(a)) < 1e-6:
                return a
            if abs(self(b)) < 1e-6:
                return b
            c = (a + b) / 2
            if abs(self(c)) < 1e-6:
                return c
            else:
                if self(a) * self(c) < 0:
                    a, b = a, c
                else:
                    a, b = c, b


class DegreeIsTooBigException(Exception):
    def __init__(self):
        pass


class QuadraticPolynomial(Polynomial):
    def __init__(self, *args):
        super().__init__(*args)
        if self.degree() > 2:
            raise DegreeIsTooBigException

    def solve(self):
        a = []
        if self.degree() == 2:
            d = self.coeffs[1] ** 2 - (4 * self.coeffs[2] * self.coeffs[0])
            if d > 0:
                a1 = (-self.coeffs[1] + d ** 0.5) / (2 * self.coeffs[2])
                a2 = (-self.coeffs[1] - d ** 0.5) / (2 * self.coeffs[2])
                a = [a1, a2]
            elif not d:
                a = [-self.coeffs[1] / (2 * self.coeffs[2])]
            else:
                return []
        elif self.degree() == 1:
            a = [-self.coeffs[0] / self.coeffs[1]]
        return a
