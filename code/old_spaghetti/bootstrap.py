class A(object):
    def __init__(self, a, b):
        self.a = a
        self.b = b
        print(a, b)


class B(A):
    def __init__(self, a, b):
        super().__init__(a, b)
        print(a, b)
        print("test")


B(1, 2)
