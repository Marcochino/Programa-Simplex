# Simplex Algorithm

## Description
Solve LP with Simplex algorithm (Maximization & Minimization)

* Simplex.py: Simplex Algorithm
* exceptions.py: User defined exceptions for specific errors.
* main.py: GUI with PyQt5

Source code is intented to work with a GUI.
<br>
To run on Terminal use Simplex.py

## Simplex Algorithm
The solution of a linear program is accomplished in two steps. In the first step, known as Phase I, a starting extreme point is found. Depending on the nature of the program this may be trivial, but in general it can be solved by applying the simplex algorithm to a modified version of the original program. The possible results of Phase I are either that a basic feasible solution is found or that the feasible region is empty. In the latter case the linear program is called infeasible. In the second step, Phase II, the simplex algorithm is applied using the basic feasible solution found in Phase I as a starting point. The possible results from Phase II are either an optimum basic feasible solution or an infinite edge on which the objective function is unbounded above

## Requirement
Install PyQt5

```
> pip3 install pyqt5 pyqt5-tools
```
## Installation

```
> python3 main.py
```
## Samples
![mainwindow](.img/mainWindow.png)
![inputwidget](.img/inputWindow.png)
