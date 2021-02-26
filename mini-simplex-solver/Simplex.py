#!/usr/bin/python3
from typing import Dict, List, Tuple
from exceptions import *
from decimal import Decimal as dc
from fractions import Fraction

phase1_iterations = list()
phase2_iterations = list()

def get_min_negative(lst: list):
    negative_values: list = list(filter(lambda x: True if x < 0 else False, lst))
    if len(negative_values):
        return lst.index(min(negative_values))
    else:
        raise NoMinNegativeValue


def get_min_positive(lst: list):
    positive_values: list = list(filter(lambda x: True if x > 0 else False, lst))
    if len(positive_values):
        return lst.index(min(positive_values))
    else:
        raise NoMinPositiveValue


def get_max_positive(lst: list):
    positive_values: list = list(filter(lambda x: True if x > 0 else False, lst))
    if len(positive_values):
        return lst.index(max(positive_values))
    else:
        raise NoMaxPostiveValue


def get_pivot(matrix: List[List[float]], costs: List[float], constraints: List[float], nature: bool):
    """Devuelve el índice del elemento pivote como una tupla
        El manejo de errores para esta función se maneja en la clase Simplex: método self.phase1()
        naturaleza: bool
            Verdadero para minimización
            Falso para maximización
    """
    m = len(matrix)  # Number of constraints
    n = len(costs)  # Number of variables

    if nature:
        # Minimisation Problem
        jpivot: int = get_min_negative(costs)  # index of minimum negative value
        ratios: list = list()
        for i in range(m):
            try:
                ratios.append(constraints[i] / matrix[i][jpivot])
            except ZeroDivisionError:
                ratios.append(-1)  # Append an arbitrary negative value so it wont be processed by get_min_positive()
        ipivot: int = get_min_positive(ratios)
        return ipivot, jpivot
    else:
        # Maximisation problem
        jpivot: int = get_max_positive(costs)
        ratios: list = list()
        for i in range(m):
            try:
                ratios.append(constraints[i] / matrix[i][jpivot])
            except ZeroDivisionError:
                ratios.append(-1)

        ipivot = get_min_positive(ratios)
        return ipivot, jpivot


class Simplex:
    def __init__(self, n: int, m: int, a: List[List[float]], b: List[float], constraints: List[str], obj_func: List[
        float], nature: bool = True):
        self.error_message: str = ""
        self.n = n  # number of variables
        self.m = m  # number of constraints
        self.a = a  # matrix  of variables coefficient
        self.b = b  # constraints values
        self.constraints = constraints  # type of constraints
        self.obj_func = obj_func  # objectiv function coefficients
        self.nature: bool = nature  # True for minimisation problem, False for maximisation

        self.unknowns: Dict[str, float] = dict()  # variables: x1, x2, ....
        self.shift_vars: Dict[str, float] = dict()  # shift variables
        self.artificial_vars: Dict[str, float] = dict()  # artificial Variables

        # Variables related to tables constructing of phase 1
        self.vars_names: List[str] = list()  # variables names
        self.base_vars_names: List[str] = list()  # base variables names
        self.table: List[List[float]] = list()  # table corresponding to matrix
        self.table_cost: float = 0  # the cost of table calculated from the objective function Za
        self.table_cost_phase2: float = 0   # we need this value when changing to phase two
        self.Za: List[List[float], float] = [list(), 0]  # The objective function of phase 1 to minimize
        self.table_coef: List[float] = list()  # coefficient of objective function in the table
        self.table_coef_phase2: List[float] = list()    # this value is also needed when changing to phase two
        self.table_constraints: List[float] = list()  # constraints of the table

        # iteration variables
        self.phase1_steps = list()
        self.phase2_steps = list()

        """Read each line of given matrix then decide what variables to add"""
        for i in range(self.m):
            if self.constraints[i] == 'lt':
                self.shift_vars['e' + str(i + 1)] = 1
            elif self.constraints[i] == 'gt':
                self.shift_vars['e' + str(i + 1)] = -1
                if self.b[i] != 0:
                    self.artificial_vars['a' + str(i + 1)] = 1
            else:
                if self.b[i] != 0:
                    self.artificial_vars['a' + str(i + 1)] = 1

        # The main process
        self._construct_table_phase1()
        if len(self.artificial_vars):
            result_phase1 = self.phase1()
            if result_phase1:
                print("CONSTRUYENDO LA TABLA DE LA FASE DOS")
                self._phase1_to_phase2_table()
                self.print_state(True)
                self.phase2()
            else:
                pass
        else:
            print("No se necesita la fase 1")
            self._construct_phase2_table()
            self.print_state(False)
            self.phase2()

    def _construct_table_phase1(self):
        """Construct the first table of phase 1"""
        self.table = self.a
        k: int = 0  # shift vars index
        l: int = len(self.shift_vars)  # artificial vars index
        p: int = len(self.shift_vars) + len(self.artificial_vars)

        for i in range(self.m):
            zeros = [0 for _ in range(p)]
            if 'e' + str(i + 1) in self.shift_vars:
                zeros[k] = self.shift_vars['e' + str(i + 1)]
                k += 1
            if 'a' + str(i + 1) in self.artificial_vars:
                zeros[l] = self.artificial_vars['a' + str(i + 1)]
                l += 1
            self.table[i] = self.table[i] + zeros
            self.table_constraints.append(self.b[i])

        # Init all vars with zeros: unknowns, shift, artificial
        for i in range(self.n):
            self.unknowns['x' + str(i + 1)] = 0
        for var in self.shift_vars:
            self.shift_vars[var] = 0
        for var in self.artificial_vars:
            self.artificial_vars[var] = 0
        self.vars_names = list(self.unknowns.keys()) + list(self.shift_vars.keys()) + list(self.artificial_vars.keys())

        # update base vars
        for i in range(self.m):
            if 'a' + str(i + 1) in self.artificial_vars:
                self.artificial_vars['a' + str(i + 1)] = self.b[i]
                self.base_vars_names.append('a' + str(i + 1))
            else:
                if 'e' + str(i + 1) in self.shift_vars:
                    self.shift_vars['e' + str(i + 1)] = self.b[i]
                    self.base_vars_names.append('e' + str(i + 1))

        # Correct the objective function
        p: int = self.n + len(self.shift_vars) + len(self.artificial_vars)  # number of variables
        za: List[float] = [0 for _ in range(p)]  # list of coefficients multiplied by -1 and summed up
        const_b: float = 0  # constant value on the objective value
        for var in self.artificial_vars:
            i: int = int(var[1]) - 1  # line number on the matrix
            tmp: List[float] = [-1 * self.table[i][j] for j in range(p)]
            for j in range(len(tmp)):
                za[j] += tmp[j]
            const_b += self.b[i]

        self.Za[0] = za[:self.n + len(self.shift_vars)]
        self.Za[1] = const_b
        self.table_coef = self.Za[0] + [0 for _ in range(len(self.artificial_vars))]
        self.table_cost = const_b

        self.table_coef_phase2 = self.obj_func[0] + [0 for _ in range(len(self.shift_vars) + len(self.artificial_vars))]
        self.table_cost_phase2 = self.obj_func[1]

        self.store_iterations(
            list(self.vars_names),
            list(self.base_vars_names),
            self.table,
            self.table_constraints,
            self.table_coef,
            self.table_cost,
            True
        )
        self.print_state(True)

    def _construct_phase2_table(self):
        """Construir tabla de problema simplex
            Este método se llama cuando no se necesitan variables artificiales.
            dado que la operación de configuración de la tabla es redundante, se llama al otro método: _construct_table_phase2
            en el constructor.
            Este método actualizará solo lo necesario
        """
        p: int = self.n + len(self.shift_vars)
        self.table_coef = self.obj_func[0] + [0 for _ in range(p-self.n)]
        self.table_cost = self.obj_func[1]

    def _phase1_to_phase2_table(self):
        """Construya la tabla de la fase 2 eliminando variables artificiales y vuelva a calcular la función objetivo"""
        p:int = self.n + len(self.shift_vars)
        tmp: List[float] = [0 for _ in range(p)]
        const_coef: float = self.obj_func[1]
        self.table = [line[:p] for line in self.table]
        self.vars_names = self.vars_names[:p]

        # use the carried around coefficient from phase one
        self.table_coef = self.table_coef_phase2[:self.n + len(self.shift_vars)]
        self.table_cost = const_coef + self._calculate_table_cost(self.unknowns, self.obj_func)

        self.store_iterations(
            list(self.vars_names),
            list(self.base_vars_names),
            self.table,
            self.table_constraints,
            self.table_coef,
            self.table_cost,
            False         
        )
        self.print_state(False)



    def phase1(self):
        """Perform phase 1 iterations"""
        # TODO: complete this function
        while self.table_cost:
            ipivot: int = 0
            jpivot: int = 0
            pivot: float = 0.0
            p: int = self.n + len(self.shift_vars) + len(self.artificial_vars)

            # Get pivot and devide its line by itself
            try:
                ipivot, jpivot = get_pivot(self.table, self.table_coef, self.table_constraints, True)
                pivot = self.table[ipivot][jpivot]
            except NoMinNegativeValue:
                print("Error: FIN DE FASE UNO")
                if self.table_cost:
                    print("El costo de la tabla no es nulo, no hay solución posible para este problema")
                    self.error_message = "El costo de la mesa no es nulo, no hay solución posible para este problema"
                return False
            except NoMinPositiveValue:
                print("Error: ALGO ESTÁ INCORRECTO CON EL CÁLCULO, PORQUE NO SE ENCUENTRA UN MIN DE RELACIÓN")
                self.error_message = "Error: ALGO ESTÁ INCORRECTO CON EL CÁLCULO, PORQUE NO SE ENCUENTRA UN MIN DE RELACIÓN"
                break
            except DegeneranceProblem:
                print("Problema de degeneración: no se puede resolver")
                self.error_message = "Problema de degeneración: no se puede resolver"
                break
            except Exception as e:
                raise e

            for i in range(p):
                self.table[ipivot][i] /= pivot
            else:
                self.table_constraints[ipivot] /= pivot

            # Update every line in the table according to pivot
            for i in range(self.m):
                if i != ipivot:
                    multiplier: float = self.table[i][jpivot]
                    for j in range(p):
                        self.table[i][j] = self.table[i][j] - self.table[ipivot][j] * multiplier
                    else:
                        self.table_constraints[i] -= self.table_constraints[ipivot] * multiplier
            else:
                # update the coeficient line
                multiplier: float = self.table_coef[jpivot]
                for i in range(p):
                    self.table_coef[i] -= self.table[ipivot][i] * multiplier

                # update the other coefficient line
                multiplier: float = self.table_coef_phase2[jpivot]
                for i in range(p):
                    self.table_coef_phase2[i] -= self.table[ipivot][i] * multiplier

            # Update Variables: leaving and entering one
            leaving: str = self.base_vars_names[ipivot]
            entering: str = self.vars_names[jpivot]
            self.base_vars_names[ipivot] = entering  # Add entering variable to base variables

            # reset all variables and update according to new table
            for var in self.unknowns: self.unknowns[var] = 0
            for var in self.shift_vars: self.shift_vars[var] = 0
            for var in self.artificial_vars: self.artificial_vars[var] = 0
            for i in range(self.m):
                var: str = self.base_vars_names[i]
                if var in self.unknowns:
                    self.unknowns[var] = self.table_constraints[i]
                if var in self.shift_vars:
                    self.shift_vars[var] = self.table_constraints[i]
                if var in self.artificial_vars:
                    self.artificial_vars[var] = self.table_constraints[i]

            # Update the table cost
            self.table_cost = self._calculate_table_cost({**self.unknowns, **self.shift_vars}, self.Za)
            self.table_cost_phase2 = self._calculate_table_cost(self.unknowns, self.obj_func)
        
            self.store_iterations(
                list(self.vars_names),
                list(self.base_vars_names),
                self.table,
                self.table_constraints,
                self.table_coef,
                self.table_cost,
                True
            )
        return True

    def phase2(self):
        """Realizar iteración de la fase 2"""

        while True:
            ipivot: int = 0
            jpivot: int = 0
            pivot: float = 0.0
            p: int = self.n + len(self.shift_vars)

            # Get pivot and devide its line by itself
            if self.nature:
                # Minimisation Problem
                try:
                    ipivot, jpivot = get_pivot(self.table, self.table_coef, self.table_constraints, True)
                except NoMinNegativeValue:
                    print("FIN DEL ALGORITMO (MINIMIZACIÓN)")
                    print("SOLUCIÓN: ", end=" ")
                    for var in self.unknowns:
                        print("{}: {}".format(var, self.unknowns[var]), end=" ")
                    print("")
                    break
                except Exception as e:
                    raise e
            else:
                # Maximisation probleme
                try:
                    ipivot, jpivot = get_pivot(self.table, self.table_coef, self.table_constraints, False)
                except NoMaxPostiveValue:
                    print('FIN DE ALGORITMO (MAXIMIZACIÓN)')
                    print("SOLUCIÓN: ", self.unknowns)
                    break
                except NoMinPositiveValue:
                    print("Error: ALGO ESTÁ INCORRECTO CON EL CÁLCULO, PORQUE NO SE ENCUENTRA UN MIN DE RELACIÓN")
                    self.error_message = "Error: ALGO ESTÁ INCORRECTO CON EL CÁLCULO, PORQUE NO SE ENCUENTRA UN MIN DE RELACIÓN"
                    break
                except Exception as e:
                    raise e
            pivot = self.table[ipivot][jpivot]

            # Devide pivot line by pivot
            for i in range(p):
                self.table[ipivot][i] /= pivot
            else:
                self.table_constraints[ipivot] /= pivot

            # Devide other lines according to pivot
            for i in range(self.m):
                if i is not ipivot:
                    multiplier: float = self.table[i][jpivot]
                    for j in range(p):
                        self.table[i][j] -= multiplier * self.table[ipivot][j]
                    else:
                        self.table_constraints[i] -= multiplier * self.table_constraints[ipivot]
            else:
                # Update table coef
                multiplier = self.table_coef[jpivot]
                for i in range(p):
                    self.table_coef[i] -= multiplier * self.table[ipivot][i]


            # Update Variables: leaving and entering one
            leaving: str = self.base_vars_names[ipivot]
            entering: str = self.vars_names[jpivot]
            self.base_vars_names[ipivot] = entering  # Add entering variable to base variables

            # reset all variables and update according to new table
            for var in self.unknowns: self.unknowns[var] = 0
            for var in self.shift_vars: self.shift_vars[var] = 0
            for i in range(self.m):
                var: str = self.base_vars_names[i]
                if var in self.unknowns:
                    self.unknowns[var] = self.table_constraints[i]
                if var in self.shift_vars:
                    self.shift_vars[var] = self.table_constraints[i]

            # Calculatae table cost
            self.table_cost = self._calculate_table_cost(self.unknowns, self.obj_func)
            self.store_iterations(
                list(self.vars_names),
                list(self.base_vars_names),
                self.table,
                self.table_constraints,
                self.table_coef,
                self.table_cost,
                False
            )
            self.print_state(False)


    ########################################################################
    # General purpose methods
    ########################################################################

    def _calculate_table_cost(self, vars_names: Dict[str, float], Za: List[float]):
        """Calcule el costo de la tabla en la fase 1"""
        res: float = Za[1]
        coef: list = Za[0]
        for key,value in list(zip(list(vars_names.keys()), coef)):
            res += vars_names[key] * value

        return res

    def store_iterations(
        self,
        vars_names,
        base_vars_names,
        table,
        table_constraints,
        table_coef,
        table_cost,
        phase
    ):

        # solve overridden variables issue
        variables_names = [var for var in vars_names]
        base_variables_names = [var for var in base_vars_names]
        matrix_table = [[var for var in line] for line in table]
        constraints = [var for var in table_constraints]
        coefs = [var for var in table_coef]
        cost = table_cost

        if phase:


            self.phase1_steps.append([
                variables_names,
                base_variables_names,
                matrix_table,
                constraints,
                coefs,
                cost
            ])
        else:
            self.phase2_steps.append([
                variables_names,
                base_vars_names,
                matrix_table,
                constraints,
                coefs,
                cost
            ])

    def print_state(self, nature: bool):
        # print("Unknows: ", self.unknowns)
        print("Desconocido: ", end="{")
        for var in self.unknowns:
            print("{}: {} ".format(var, self.unknowns[var]), end="")
        print("}")

        # print("shift Variables:", self.shift_vars)
        print("Variaciones de cambio: ", end="{")
        for var in self.shift_vars:
            print("{}: {} ".format(var, self.shift_vars[var]), end="")
        print("}")

        if nature:
            # print("Artificial vars:", self.artificial_vars)
            print("Variable artificial: ", end="{")
            for var in self.artificial_vars:
                print("{}: {}, ".format(var, self.artificial_vars[var]), end="")
            print("}")

        # print("*", self.vars_names, "constraints")
        print("*", end=" | ")
        for var in self.vars_names:
            print("{}".format(var), end="\t")
        else:
            print(" | Restricciones")

        for i in range(self.m):
            # print(self.base_vars_names[i], self.table[i], self.table_constraints[i])
            print(self.base_vars_names[i], end=" | ")
            for var in self.table[i]:
                print("{}".format(var), end="\t")
            else:
                print(" | {}".format(self.table_constraints[i]))

        # print("costs", self.table_coef, self.table_cost)
        print("Z ", end=" | ")
        for var in self.table_coef:
            print("{}".format(var), end="\t")
        else:
            print(" | {}".format(self.table_cost))

        print("=" * 20)


if __name__ == '__main__':
    # problem with aritificl variables
    # n: int = 2  # nombre de variables
    # m: int = 4  # nombre de contrainte
    # a = [
    #     [Fraction(10), Fraction(5)],
    #     [Fraction(2), Fraction(3)],
    #     [Fraction(1), Fraction(0)],
    #     [Fraction(0), Fraction(1)]
    # ]
    # b = [Fraction(200), Fraction(60), Fraction(12), Fraction(6)]
    # const = ['lt', 'eq', 'lt', 'gt']
    # objective_function = [Fraction(1000), Fraction(2000)]
    # simplex = Simplex(n, m, a, b, const, [objective_function, 0], False)
    #
    # print("+"*100)
    # print("ANOTHER EXAMPLE")
    # print("+"*100)
    # probleme with aritficial variables
    # n: int = 2
    # m: int = 4
    # a = [
    #     [2, 1],
    #     [1, 1],
    #     [5, 4],
    #     [1, 2]
    # ]
    # b = [600, 225, 1000, 150]
    # const = ['lt', 'lt', 'lt', 'lt']
    # objective_function = [3, 4]

    # probelem without artificial variables
    # n: int = 2
    # m: int = 3
    # a = [
    #     [1, 0],
    #     [0, 2],
    #     [3, 2]
    # ]
    # b = [4, 12, 18]
    # const = ['lt', 'lt', 'lt']
    # objective_function = [3, 5]
    n:int = 2
    m: int = 4
    a = [
        [Fraction(1, 10), Fraction(0)],
        [Fraction(0), Fraction(1, 10)],
        [Fraction(1, 10), Fraction(2, 10)],
        [Fraction(2, 10), Fraction(1, 10)]
    ]
    b = [Fraction(4, 10), Fraction(6, 10), Fraction(2), Fraction(17, 10)]
    const = ['gt', 'gt', 'gt', 'gt']
    objective_function = [Fraction(100), Fraction(40)]
    simplex = Simplex(n, m, a, b, const, [objective_function, 0], True)
    print(len(simplex.phase1_steps))
    for iteration in simplex.phase1_steps:
        print(iteration)
        print("="*100)

    # print("pivot Minimisation: ", get_pivot([
    #     [10, 5, 1, 0, 0, 0, 0, 0],
    #     [2, 3, 0, 1, 0, 0, 0, 0],
    #     [1, 0, 0, 0, 1, 0, 0],
    #     [0, 1, 0, 0, 0, -1, -1]
    # ], [-2, -4, 0, 0, 0, 1, 0], [200, 60,12, 6], True))
    #
    # print("pivot maximisation: ", get_pivot([
    #     [1, 2, 3, 1, 0, 0],
    #     [15, 21, 30, 0, 1, 0],
    #     [1, 1, 1, 0, 0, 1]
    # ], [87, 147, 258, 0, 0, 0], [90, 1260, 84, 0], False))
    print("Fin de la ejecución")
