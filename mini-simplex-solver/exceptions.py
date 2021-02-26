#!/usr/bin/python3


class NoMinNegativeValue(Exception):
    """No Minimum Negative value is found in minimisation probleme"""

    pass


class NoMinPositiveValue(Exception):
    """No Minimum Positive value is found on ratios"""

    pass


class NoMaxPostiveValue(Exception):
    """No Maximum Positive Value is found in maximisation probleme"""

    pass
class DegeneranceProblem(Exception):
    """It's a degenerance problem"""

    pass