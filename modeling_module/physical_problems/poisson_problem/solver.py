from fenics import *
import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.tri as mtri
import sympy as sym


# from mpl_toolkits.mplot3d import *
# from matplotlib import cm
# import matplotlib.animation as animation
############################
#####УРАВНЕНИЕ ПУАССОНА#####
############################
path = 'modeling_module/physical_problems/math_problems/'

class Task_maker():

    def __init__(self, config):
        self.config = config

        self.name_vector = []
        self.name_scalar = []
        self.description_vector = []
        self.description_scalar = []
        self.type_vector = []
        self.type_scalar = []
        self.vector_condition = []
        self.scalar_condition = []
        self.notation_scalar = []
        self.source = "0"
        self.kappa = "1"

        if self.config.mesh_name:
             self.mesh = path + 'mesh/' + self.config.mesh_name

        if self.config.source:
            self.source = self.config.source

        # if self.config.kappa:
        #     self.kappa = self.config.kappa

        # if self.config.vector_functions:
        #     self.conditions = self.config.vector_functions
        #     for boundaries in self.conditions:
        #         if boundaries.boundary_name:
        #             self.name_vector.append(boundaries.boundary_name)
        #         if boundaries.description:
        #             self.description_vector.append(boundaries.description)
        #         if boundaries.condition_type:
        #             self.type_vector.append(boundaries.condition_type)
        #         if boundaries.condition:
        #             self.vector_condition.append(boundaries.condition)

        if self.config.scalar_functions:
            self.conditions = self.config.scalar_functions
            for boundaries in self.conditions:
                if boundaries.name:
                    self.name_scalar.append(boundaries.name)
                if boundaries.description:
                    self.description_scalar.append(boundaries.description)
                if boundaries.condition_type:
                    self.type_scalar.append(boundaries.condition_type)
                if boundaries.condition:
                    self.scalar_condition.append(boundaries.condition)
                if boundaries.notation:
                    self.notation_scalar.append(boundaries.notation)
                else:
                    self.notation_scalar.append("EXPR")


class BVP_solver():

    def __init__(self, task):

        self.mesh = Mesh(task.mesh)

        self.V = FunctionSpace(self.mesh, 'P', 1)
        #объявление искомых функций и пробных функций. Они являются частью V.
        self.u = Function(self.V)
        self.v = TestFunction(self.V)

        self.Dc = []
        self.Nc = []
        self.bx = []

        self.boundary_parts = MeshFunction('size_t', self.mesh, self.mesh.topology().dim() - 1)

        ds = Measure('ds', domain = self.mesh, subdomain_data = self.boundary_parts)

        k=0
        m=0
        for i in task.name_scalar:
            if task.type_scalar[k]=="Dirichlet":
                if task.notation_scalar[k]=="SYM":
                    u_D = sym.sympify(task.scalar_condition[k])
                    x, y, z, t = sym.symbols('x[0], x[1], x[2], t')
                    sub = [('x',x),('y',y),('z',z),('t',t)]
                    u_D = u_D.subs(sub)
                    u_code = sym.printing.ccode(u_D)
                    u_D = Expression(u_code, degree=2)
                else:
                    u_D = Expression(task.scalar_condition[k], degree=2)
                DC = DirichletBC(self.V, u_D, task.description_scalar[k])
                self.Dc.append(DC)
            if task.type_scalar[k]=="Neumann":
                if task.notation_scalar[k]=="SYM":
                    g_N = sym.sympify(task.scalar_condition[k])
                    x, y, z, t = sym.symbols('x[0], x[1], x[2], t')
                    sub = [('x',x),('y',y),('z',z),('t',t)]
                    g_N = g_N.subs(sub)
                    g_code = sym.printing.ccode(g_N)
                    g = Expression(g_code, degree=2)
                else:
                    g = Expression(task.scalar_condition[k], degree=2)
                bx1 = CompiledSubDomain(task.description_scalar[k])
                self.bx.append(bx1)
                self.bx[m].mark(self.boundary_parts, m)
                self.Nc.append(g*self.v*ds(m))
                m = m+1
            k = k+1

            #источник в правой части уравнения Пуассона
        self.f = Expression(task.source, degree=2, u=self.u)

        self.kappa1 = Expression(task.kappa, degree=2,u=self.u)


    def Solving_eq(self):
        #постановка вариационной задачи и ее решение с граничными условиями
        F = self.kappa1*dot(grad(self.u), grad(self.v))*dx - self.f*self.v*dx - sum(self.Nc)
        solve(F == 0, self.u, self.Dc)

        # Save solution to file in VTK format
        vtkfile = File(path+'results/'+'solution.pvd')
        vtkfile << self.u


#пример задания граничных условий через привычные переменные x,y на языке Python
# x, y = sym.symbols('x[0], x[1]')
# u_D = 1 + x + 2*y
#
# #пример сложного выражения для правой части уравнения Пуассона в терминах sympy
# f = - sym.diff(q(u_D)*sym.diff(u_D, x), x) - sym.diff(q(u_D)*sym.diff(u_D, y), y)
# f = sym.simplify(f)
#
# #конвертация в код C++
# u_code = sym.printing.ccode(u_D)
# f_code = sym.printing.ccode(f)
# print('u =', u_code)
# print('f =', f_code)
#
# #граничные условия Дирихле
# u_D = Expression(u_code, degree=2)
# def boundary(x, on_boundary):
#     return on_boundary
# bc = DirichletBC(V, u_D, boundary)
#
# f = Expression(f_code, degree=1)
# F = q(u)*dot(grad(u), grad(v))*dx - f*v*dx
#
# solve(F == 0, u, bc)
#
# # Plot solution and mesh
# plot(u)
# plot(mesh)
#
# # Save solution to file in VTK format
# vtkfile = File('ex5/solution.pvd')
# vtkfile << u
